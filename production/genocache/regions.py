"""Bounded, strand-aware region proposals from uncertain embedding hits.

This module proposes reference windows. It does not align bases, assign MAPQ,
prove genome-wide uniqueness, or emit SAM. An ANN reference window is an interval
of uncertainty, not an exact point anchor. All query coordinates supplied to the
API refer to the original forward read, including negative-strand hits.

For an oriented query interval beginning at q and spanning l bases, a reference
window [r0, r1) supports read-start diagonals in [r0-q, r1-l-q], inclusive. The
explicit positional_slop and indel_slack expand those bounds. A sweep computes
simultaneous support from distinct query intervals and their union coverage.
Several ANN neighbors or seed-ID aliases for one query interval get one vote.

Runtime is O(M log M + K*M), memory O(M), where inspected hits M and output
candidates K are bounded. K cannot exceed 64. Input and output budget truncation
are reported; the function never claims candidate-set completeness.
"""

from __future__ import annotations

import heapq
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import islice
from numbers import Integral, Real
from typing import Any, Iterable, Mapping


MAX_HITS = 100_000
MAX_CANDIDATES = 64
MAX_COORDINATE = (1 << 63) - 1


@dataclass(frozen=True, slots=True)
class Anchor:
    read_id: str
    read_length: int
    query_start: int
    query_end: int
    seed_id: str | int
    chrom: str
    ref_start: int
    ref_end: int
    strand: str
    similarity: float


@dataclass(frozen=True, slots=True)
class CandidateRegion:
    """A reference proposal, with inclusive uncertain read-start bounds.

    start/end are a 0-based, half-open reference window. diag_low/diag_high
    describe possible starts under the selected support hypothesis; negative
    values can describe reads extending beyond a contig boundary. score is the
    union of supported query bases, not an alignment score or probability.
    seed_ids contains one representative ID per unique query interval; reused
    IDs can therefore appear more than once. query_intervals disambiguates them.
    """

    chrom: str
    strand: str
    start: int
    end: int
    distinct_seeds: int
    query_coverage_bases: int
    query_coverage_fraction: float
    score: float
    seed_ids: tuple[str, ...]
    diag_low: int
    diag_high: int
    mean_similarity: float
    query_intervals: tuple[tuple[int, int], ...]
    clipped_to_contig: bool

    def to_dict(self) -> dict[str, Any]:
        return {
            "chrom": self.chrom,
            "strand": self.strand,
            "start": self.start,
            "end": self.end,
            "distinct_seeds": self.distinct_seeds,
            "query_coverage_bases": self.query_coverage_bases,
            "query_coverage_fraction": self.query_coverage_fraction,
            "score": self.score,
            "seed_ids": list(self.seed_ids),
            "diag_low": self.diag_low,
            "diag_high": self.diag_high,
            "mean_similarity": self.mean_similarity,
            "query_intervals": [list(interval) for interval in self.query_intervals],
            "clipped_to_contig": self.clipped_to_contig,
        }


@dataclass(frozen=True, slots=True)
class ProposalResult:
    proposals: tuple[CandidateRegion, ...]
    diagnostics: dict[str, Any]

    def to_dict(self) -> dict[str, Any]:
        return {
            "proposals": [proposal.to_dict() for proposal in self.proposals],
            "diagnostics": dict(self.diagnostics),
        }

    as_dict = to_dict


@dataclass(frozen=True, slots=True)
class _Hit:
    anchor: Anchor
    low: int
    high: int

    @property
    def query_interval(self) -> tuple[int, int]:
        return self.anchor.query_start, self.anchor.query_end


@dataclass(frozen=True, slots=True)
class _Cell:
    chrom: str
    strand: str
    low: int
    high: int
    distinct_seeds: int
    coverage: int
    mean_similarity: float

    @property
    def point(self) -> int:
        return self.low + (self.high - self.low) // 2


class _InvalidAnchor(ValueError):
    pass


def _integer(value: Any, label: str, *, minimum: int = 0, maximum: int = MAX_COORDINATE) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise ValueError(f"{label} must be an integer")
    result = int(value)
    if not minimum <= result <= maximum:
        raise ValueError(f"{label} must be between {minimum} and {maximum}")
    return result


def _validate_anchor(raw: Anchor | Mapping[str, Any], contigs: Mapping[str, int]) -> Anchor:
    if isinstance(raw, Anchor):
        values = {field: getattr(raw, field) for field in Anchor.__dataclass_fields__}
    elif isinstance(raw, Mapping):
        try:
            values = {field: raw[field] for field in Anchor.__dataclass_fields__}
        except KeyError as exc:
            raise _InvalidAnchor("missing_fields") from exc
    else:
        raise _InvalidAnchor("not_an_anchor")

    read_id = values["read_id"]
    if not isinstance(read_id, str) or not read_id.strip():
        raise _InvalidAnchor("invalid_read_id")
    chrom = values["chrom"]
    if not isinstance(chrom, str) or chrom not in contigs:
        raise _InvalidAnchor("unknown_contig")
    if values["strand"] not in ("+", "-"):
        raise _InvalidAnchor("invalid_strand")

    seed_id = values["seed_id"]
    if isinstance(seed_id, bool):
        raise _InvalidAnchor("invalid_seed_id")
    if isinstance(seed_id, Integral):
        if seed_id < 0:
            raise _InvalidAnchor("invalid_seed_id")
        seed_id = str(int(seed_id))
    elif not isinstance(seed_id, str) or not seed_id.strip():
        raise _InvalidAnchor("invalid_seed_id")

    try:
        length = _integer(values["read_length"], "read_length", minimum=1)
        qs = _integer(values["query_start"], "query_start")
        qe = _integer(values["query_end"], "query_end")
        rs = _integer(values["ref_start"], "ref_start")
        re = _integer(values["ref_end"], "ref_end")
    except ValueError as exc:
        raise _InvalidAnchor("invalid_integer_coordinate") from exc
    if not 0 <= qs < qe <= length:
        raise _InvalidAnchor("invalid_query_interval")
    if not 0 <= rs < re <= contigs[chrom]:
        raise _InvalidAnchor("invalid_reference_interval")

    similarity = values["similarity"]
    if isinstance(similarity, bool) or not isinstance(similarity, Real):
        raise _InvalidAnchor("invalid_similarity")
    similarity = float(similarity)
    if not math.isfinite(similarity):
        raise _InvalidAnchor("nonfinite_similarity")

    return Anchor(
        read_id, length, qs, qe, seed_id, chrom, rs, re, values["strand"], similarity
    )


class _QueryCoverage:
    """Dynamic union length of query intervals in O(log M) per update."""

    def __init__(self, query_intervals: Iterable[tuple[int, int]]):
        self.coordinates = sorted({coordinate for pair in query_intervals for coordinate in pair})
        self.indices = {coordinate: i for i, coordinate in enumerate(self.coordinates)}
        self.segments = len(self.coordinates) - 1
        self.count = [0] * (4 * self.segments + 4)
        self.length = [0] * (4 * self.segments + 4)

    @property
    def covered(self) -> int:
        return self.length[1]

    def update(self, interval: tuple[int, int], delta: int) -> None:
        first, last = (self.indices[coordinate] for coordinate in interval)
        self._update(1, 0, self.segments, first, last, delta)

    def _update(self, node: int, left: int, right: int, first: int, last: int, delta: int) -> None:
        if first <= left and right <= last:
            self.count[node] += delta
        else:
            middle = (left + right) // 2
            if first < middle:
                self._update(node * 2, left, middle, first, last, delta)
            if last > middle:
                self._update(node * 2 + 1, middle, right, first, last, delta)
        if self.count[node] > 0:
            self.length[node] = self.coordinates[right] - self.coordinates[left]
        elif right - left == 1:
            self.length[node] = 0
        else:
            self.length[node] = self.length[node * 2] + self.length[node * 2 + 1]


def _sweep(
    chrom: str,
    strand: str,
    hits: list[_Hit],
    read_length: int,
    min_seeds: int,
    min_coverage: float,
) -> list[_Cell]:
    # Inclusive integer start intervals become half-open sweep intervals.
    events: list[tuple[int, int, int]] = []
    for index, hit in enumerate(hits):
        events.append((hit.low, 1, index))
        events.append((hit.high + 1, -1, index))
    events.sort()

    coverage = _QueryCoverage(hit.query_interval for hit in hits)
    counts: Counter[tuple[int, int]] = Counter()
    heaps: dict[tuple[int, int], list[tuple[float, str, int]]] = defaultdict(list)
    active = [False] * len(hits)
    distinct = 0
    # Scaling prevents overflow even for finite but unusually large similarities.
    scale = max(1.0, max(abs(hit.anchor.similarity) for hit in hits))
    scaled_similarity_sum = 0.0
    cells: list[_Cell] = []

    def best_similarity(query_interval: tuple[int, int]) -> float:
        heap = heaps[query_interval]
        while heap and not active[heap[0][2]]:
            heapq.heappop(heap)
        return -heap[0][0] if heap else 0.0

    cursor = 0
    while cursor < len(events):
        coordinate = events[cursor][0]
        while cursor < len(events) and events[cursor][0] == coordinate:
            _, delta, index = events[cursor]
            hit = hits[index]
            query_interval = hit.query_interval
            previous_best = best_similarity(query_interval)
            if delta == 1:
                if counts[query_interval] == 0:
                    distinct += 1
                    coverage.update(query_interval, 1)
                counts[query_interval] += 1
                active[index] = True
                heapq.heappush(
                    heaps[query_interval],
                    (-hit.anchor.similarity, str(hit.anchor.seed_id), index),
                )
            else:
                active[index] = False
                counts[query_interval] -= 1
                if counts[query_interval] == 0:
                    distinct -= 1
                    coverage.update(query_interval, -1)
            current_best = best_similarity(query_interval)
            scaled_similarity_sum += current_best / scale - previous_best / scale
            cursor += 1

        if cursor == len(events) or distinct < min_seeds:
            continue
        covered_bases = coverage.covered
        if covered_bases / read_length < min_coverage:
            continue
        next_coordinate = events[cursor][0]
        mean_similarity = max(-1.0, min(1.0, scaled_similarity_sum / distinct)) * scale
        cells.append(
            _Cell(
                chrom, strand, coordinate, next_coordinate - 1,
                distinct, covered_bases, mean_similarity,
            )
        )
    return cells


def _union_length(sorted_intervals: Iterable[tuple[int, int]]) -> int:
    total = 0
    previous_end = 0
    for start, end in sorted_intervals:
        total += max(0, end - max(start, previous_end))
        previous_end = max(previous_end, end)
    return total


def _build_region(
    cell: _Cell, hits: list[_Hit], contig_length: int, read_length: int, indel_slack: int
) -> CandidateRegion:
    # hits are sorted by original query interval. Assignment preserves that order.
    best: dict[tuple[int, int], _Hit] = {}
    for hit in hits:
        if not hit.low <= cell.point <= hit.high:
            continue
        query_interval = hit.query_interval
        previous = best.get(query_interval)
        rank = (-hit.anchor.similarity, hit.anchor.ref_start, hit.anchor.ref_end, str(hit.anchor.seed_id))
        if previous is None:
            best[query_interval] = hit
        else:
            previous_rank = (
                -previous.anchor.similarity, previous.anchor.ref_start,
                previous.anchor.ref_end, str(previous.anchor.seed_id),
            )
            if rank < previous_rank:
                best[query_interval] = hit

    query_intervals = tuple(best)
    coverage = _union_length(query_intervals)
    if len(best) != cell.distinct_seeds or coverage != cell.coverage:
        raise RuntimeError("Region sweep support invariant failed")

    low = max(hit.low for hit in best.values())
    high = min(hit.high for hit in best.values())
    if not low <= cell.point <= high:
        raise RuntimeError("Region sweep uncertainty invariant failed")
    raw_end = high + read_length + indel_slack
    start, end = max(0, low), min(contig_length, raw_end)
    if start >= end:
        raise RuntimeError("Region proposal has no overlap with reference")

    scale = max(1.0, max(abs(hit.anchor.similarity) for hit in best.values()))
    scaled_mean = math.fsum(hit.anchor.similarity / scale for hit in best.values()) / len(best)
    mean_similarity = max(-1.0, min(1.0, scaled_mean)) * scale
    return CandidateRegion(
        chrom=cell.chrom,
        strand=cell.strand,
        start=start,
        end=end,
        distinct_seeds=len(best),
        query_coverage_bases=coverage,
        query_coverage_fraction=coverage / read_length,
        score=float(coverage),
        seed_ids=tuple(str(hit.anchor.seed_id) for hit in best.values()),
        diag_low=low,
        diag_high=high,
        mean_similarity=mean_similarity,
        query_intervals=query_intervals,
        clipped_to_contig=low < 0 or raw_end > contig_length,
    )


def propose_regions(
    anchors: Iterable[Anchor | Mapping[str, Any]],
    contig_lengths: Mapping[str, int],
    *,
    positional_slop: int = 64,
    indel_slack: int = 0,
    max_candidates: int = 8,
    min_seeds: int = 2,
    min_coverage: float = 0.1,
    max_hits: int = 10_000,
) -> ProposalResult:
    """Propose candidate regions for exactly one read, with explicit limitations.

    Invalid hits are rejected with reason counts. Mixed valid read identities or
    read lengths raise ValueError so evidence from separate reads cannot mix.
    All coordinates are integer, 0-based half-open intervals. Positional slop
    expands each start-uncertainty endpoint. Indel slack is an additional absolute
    bound in bases, cannot exceed the read length, and extends the proposal end.
    Reference windows shorter than a query seed are rejected unless the declared
    slop/slack makes a nonempty feasible start interval.

    The endpoint sweep requires actual common overlap, so a transitive chain of
    nearby but incompatible intervals does not become one strong candidate.
    Rankings prioritize union query coverage, then distinct query intervals,
    then the mean of the best active similarity per query interval. Neither the
    rank nor score is a calibrated correctness probability.

    At most max_hits inputs plus one lookahead are consumed. If truncated, only
    that submitted prefix participates: permutation invariance is guaranteed
    only when the input budget is not exceeded. At most max_candidates support
    sets are reconstructed. Covered sweep cells can be suppressed by a selected
    proposal's diagonal interval; genuinely omitted cells trigger diagnostics.
    """
    positional_slop = _integer(positional_slop, "positional_slop")
    indel_slack = _integer(indel_slack, "indel_slack")
    max_hits = _integer(max_hits, "max_hits", minimum=1, maximum=MAX_HITS)
    max_candidates = _integer(max_candidates, "max_candidates", minimum=1, maximum=MAX_CANDIDATES)
    min_seeds = _integer(min_seeds, "min_seeds", minimum=1, maximum=MAX_HITS)
    if isinstance(min_coverage, bool) or not isinstance(min_coverage, Real):
        raise ValueError("min_coverage must be a finite fraction between 0 and 1")
    min_coverage = float(min_coverage)
    if not math.isfinite(min_coverage) or not 0 <= min_coverage <= 1:
        raise ValueError("min_coverage must be a finite fraction between 0 and 1")
    if not isinstance(contig_lengths, Mapping):
        raise ValueError("contig_lengths must be a mapping")
    contigs: dict[str, int] = {}
    for name, length in contig_lengths.items():
        if not isinstance(name, str) or not name:
            raise ValueError("Contig names must be nonempty strings")
        contigs[name] = _integer(length, "contig length")

    prefix = list(islice(iter(anchors), max_hits + 1))
    input_truncated = len(prefix) > max_hits
    inspected = prefix[:max_hits]
    rejected: Counter[str] = Counter()
    identity: tuple[str, int] | None = None
    duplicate_hits = 0
    accepted_raw = 0
    unique: dict[tuple[Any, ...], _Hit] = {}
    total_slop = positional_slop + indel_slack

    for raw in inspected:
        try:
            anchor = _validate_anchor(raw, contigs)
        except _InvalidAnchor as exc:
            rejected[str(exc)] += 1
            continue
        current_identity = anchor.read_id, anchor.read_length
        if identity is not None and current_identity != identity:
            raise ValueError("propose_regions accepts exactly one read_id and read_length per call")
        identity = current_identity
        if indel_slack > anchor.read_length:
            raise ValueError("indel_slack cannot exceed the read length")
        oriented_query = (
            anchor.query_start if anchor.strand == "+" else anchor.read_length - anchor.query_end
        )
        span = anchor.query_end - anchor.query_start
        low = anchor.ref_start - oriented_query - total_slop
        high = anchor.ref_end - span - oriented_query + total_slop
        if high < low:
            rejected["reference_window_too_short"] += 1
            continue
        # Allow boundary-spanning reads, but every proposed read must overlap the
        # reference by at least one base under the declared indel allowance.
        low = max(low, 1 - anchor.read_length - indel_slack)
        high = min(high, contigs[anchor.chrom] - 1)
        if high < low:
            rejected["no_possible_reference_overlap"] += 1
            continue
        accepted_raw += 1
        key = (
            anchor.chrom, anchor.strand, anchor.query_start, anchor.query_end,
            anchor.ref_start, anchor.ref_end,
        )
        previous = unique.get(key)
        hit = _Hit(anchor, low, high)
        if previous is not None:
            duplicate_hits += 1
            if (-anchor.similarity, str(anchor.seed_id)) >= (
                -previous.anchor.similarity, str(previous.anchor.seed_id)
            ):
                continue
        unique[key] = hit

    groups: dict[tuple[str, str], list[_Hit]] = defaultdict(list)
    for key in sorted(unique):
        hit = unique[key]
        groups[(hit.anchor.chrom, hit.anchor.strand)].append(hit)
    read_length = identity[1] if identity is not None else None
    cells: list[_Cell] = []
    for (chrom, strand), hits in groups.items():
        cells.extend(_sweep(chrom, strand, hits, read_length, min_seeds, min_coverage))
    cells.sort(
        key=lambda cell: (
            -cell.coverage, -cell.distinct_seeds, -cell.mean_similarity,
            cell.chrom, cell.strand, cell.low, cell.high,
        )
    )

    proposals: list[CandidateRegion] = []
    redundant_cells = 0
    budget_omitted_cells = 0
    for cell in cells:
        if any(
            proposal.chrom == cell.chrom and proposal.strand == cell.strand
            and proposal.diag_low <= cell.low and cell.high <= proposal.diag_high
            for proposal in proposals
        ):
            redundant_cells += 1
            continue
        if len(proposals) >= max_candidates:
            budget_omitted_cells += 1
            continue
        proposals.append(
            _build_region(
                cell, groups[(cell.chrom, cell.strand)], contigs[cell.chrom],
                read_length, indel_slack,
            )
        )

    truncation_reasons = []
    if input_truncated:
        truncation_reasons.append("max_hits")
    if budget_omitted_cells:
        truncation_reasons.append("max_candidates")
    diagnostics = {
        "authoritative": False,
        "purpose": "candidate_region_proposals_only",
        "score_definition": "union_query_coverage_bases",
        "read_id": identity[0] if identity else None,
        "read_length": read_length,
        "anchors_inspected": len(inspected),
        "anchors_seen_at_least": len(prefix),
        "anchors_accepted_before_deduplication": accepted_raw,
        "anchors_accepted": len(unique),
        "anchors_rejected": sum(rejected.values()),
        "rejected_by_reason": dict(sorted(rejected.items())),
        "duplicate_hits_removed": duplicate_hits,
        "unique_query_intervals": len({hit.query_interval for hit in unique.values()}),
        "groups": len(groups),
        "qualifying_sweep_intervals": len(cells),
        "redundant_sweep_intervals": redundant_cells,
        "budget_omitted_sweep_intervals": budget_omitted_cells,
        "truncated": bool(truncation_reasons),
        "truncation_reasons": truncation_reasons,
        "max_hits": max_hits,
        "max_candidates": max_candidates,
        "positional_slop": positional_slop,
        "indel_slack": indel_slack,
        "min_seeds": min_seeds,
        "min_coverage": min_coverage,
    }
    return ProposalResult(tuple(proposals), diagnostics)
