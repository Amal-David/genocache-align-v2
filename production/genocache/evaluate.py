"""Truth-aware mapping and retrieval evaluation with explicit denominators.

Truth is supplied independently: minimap2 agreement alone is not biological truth.
All coordinate fields in truth/proposals are zero-based. Unmapped/missing reads
stay in the denominator. An ambiguous read may have multiple valid origins.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from collections.abc import Mapping
import json
from pathlib import Path
import tempfile
import os

from .model import wilson_interval
from .retrieval import file_sha256, json_lines


def _truth_integer(value, field: str, *, minimum: int = 1) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value < (1 << 63):
        raise ValueError(f"{field} must be an integer between {minimum} and {(1 << 63) - 1}")
    return value


def load_truth(path: Path) -> dict:
    truth = {}
    for record in json_lines(path):
        name = record.get("read_id")
        if not isinstance(name, str) or not name.strip() or name in truth:
            raise ValueError("truth read_id must be a unique nonempty string")
        if not isinstance(record.get("origins"), list):
            raise ValueError("truth origins must be a list; [] means deliberately unmappable")
        stratum = record.get("stratum", "unclassified")
        if not isinstance(stratum, str) or not stratum.strip():
            raise ValueError("truth stratum must be a nonempty string")
        for field in ("read_length", "reference_span"):
            if field in record:
                _truth_integer(record[field], field)
        for origin in record["origins"]:
            if (
                not isinstance(origin, dict)
                or not isinstance(origin.get("chrom"), str)
                or not origin["chrom"] or any(c.isspace() for c in origin["chrom"])
                or origin.get("strand") not in ("+", "-")
            ):
                raise ValueError("invalid truth origin")
            _truth_integer(origin.get("start"), "truth origin start", minimum=0)
        truth[name] = record
    if not truth:
        raise ValueError("empty truth set")
    return truth


def _summarize(rows: list[dict], field: str) -> dict:
    strata = defaultdict(list)
    for row in rows:
        strata[row["stratum"]].append(row)
    def metrics(items):
        eligible = [r for r in items if r["has_origin"]]
        correct = sum(r[field] for r in eligible)
        count = len(eligible)
        return {
            "reads": len(items),
            "reads_with_known_origins": count,
            "correct": correct,
            "recall": correct / count if count else None,
            "recall_wilson_95": wilson_interval(correct, count) if count else None,
        }
    return {"overall": metrics(rows), "strata": {k: metrics(v) for k, v in sorted(strata.items())}}


def evaluate_alignments(alignment_file: Path, truth_file: Path, *, tolerance: int = 100) -> dict:
    import pysam

    if not isinstance(tolerance, int) or isinstance(tolerance, bool) or tolerance < 0:
        raise ValueError("tolerance must be a nonnegative integer")
    truth = load_truth(truth_file)
    primary = {}
    unknown = set()
    secondary = supplementary = 0
    with pysam.AlignmentFile(str(alignment_file)) as stream:
        for alignment in stream.fetch(until_eof=True):
            if alignment.query_name not in truth:
                unknown.add(alignment.query_name)
            if alignment.is_secondary:
                secondary += 1
                continue
            if alignment.is_supplementary:
                supplementary += 1
                continue
            if alignment.query_name in primary:
                raise ValueError("duplicate primary read names: use unique QNAMEs in evaluation fixtures")
            primary[alignment.query_name] = alignment
    if unknown:
        raise ValueError(f"{len(unknown)} output read IDs are absent from truth")
    rows = []
    for name, item in truth.items():
        hit = primary.get(name)
        mapped = hit is not None and not hit.is_unmapped
        origins = item["origins"]
        correct = False
        if mapped:
            strand = "-" if hit.is_reverse else "+"
            correct = any(
                hit.reference_name == origin["chrom"]
                and strand == origin["strand"]
                and abs(hit.reference_start - origin["start"]) <= tolerance
                for origin in origins
            )
        rows.append({
            "read_id": name,
            "stratum": item.get("stratum", "unclassified"),
            "has_origin": bool(origins),
            "mapped": mapped,
            "primary_record_present": hit is not None,
            "correct_origin": correct,
            "mapq": hit.mapping_quality if mapped else None,
            "chrom": hit.reference_name if mapped else None,
            "start": hit.reference_start if mapped else None,
            "strand": ("-" if hit.is_reverse else "+") if mapped else None,
        })
    report = _summarize(rows, "correct_origin")
    report["overall"].update({
        "mapped": sum(r["mapped"] for r in rows),
        "mapping_rate": sum(r["mapped"] for r in rows) / len(rows),
        "missing_primary_records": sum(not r["primary_record_present"] for r in rows),
        "false_mapped_known_unmappable": sum(r["mapped"] and not r["has_origin"] for r in rows),
        "secondary_records": secondary,
        "supplementary_records": supplementary,
        "mapped_with_unavailable_mapq": sum(r["mapped"] and r["mapq"] == 255 for r in rows),
    })
    # SAM 255 means mapping quality is unavailable; it is not high confidence.
    high_mapq = [r for r in rows if r["mapped"] and 20 <= r["mapq"] < 255]
    report["mapq_ge_20"] = {
        "reads": len(high_mapq),
        "wrong": sum(not r["correct_origin"] for r in high_mapq),
        "empirical_error_rate": (
            sum(not r["correct_origin"] for r in high_mapq) / len(high_mapq)
            if high_mapq else None
        ),
    }
    report.update({
        "schema": "genocache.mapping-evaluation.v1",
        "truth_sha256": file_sha256(truth_file),
        "alignment_sha256": file_sha256(alignment_file),
        "position_tolerance_bp": tolerance,
        "metric": "primary reference start and strand agree with a supplied valid origin",
        "cigar_or_variant_accuracy_evaluated": False,
        "per_read": rows,
    })
    return report


def evaluate_retrieval(
    anchors_file: Path, truth_file: Path, contig_lengths: dict, *,
    positional_slop: int = 256, max_candidates: int = 8, min_seeds: int = 2,
) -> dict:
    from .regions import propose_regions

    truth = load_truth(truth_file)
    if not isinstance(contig_lengths, Mapping) or not contig_lengths:
        raise ValueError("contig_lengths must be a nonempty mapping")
    for chrom, length in contig_lengths.items():
        if not isinstance(chrom, str) or not chrom or any(c.isspace() for c in chrom):
            raise ValueError("contig names must be nonempty strings without whitespace")
        _truth_integer(length, "contig length")
    for item in truth.values():
        _truth_integer(item.get("read_length"), "read_length")
        if item["origins"]:
            span = _truth_integer(item.get("reference_span"), "reference_span")
            for origin in item["origins"]:
                if (
                    origin["chrom"] not in contig_lengths
                    or origin["start"] + span > contig_lengths[origin["chrom"]]
                ):
                    raise ValueError("truth origin reference span does not fit the declared contig")
    grouped = defaultdict(list)
    # This diagnostic is for a declared held-out panel, not unbounded service ingestion.
    for record in json_lines(anchors_file):
        name = record.get("read_id")
        if not isinstance(name, str) or name not in truth:
            raise ValueError("retrieval output contains a read absent from truth")
        if _truth_integer(record.get("read_length"), "anchor read_length") != truth[name]["read_length"]:
            raise ValueError("anchor read_length differs from independent truth read_length")
        if len(grouped[name]) >= 10000:
            raise ValueError("more than 10,000 hits for one read; evaluate with bounded retrieval K")
        grouped[name].append(record)
    rows = []
    for name, item in truth.items():
        anchors = grouped.get(name, [])
        result = propose_regions(
            anchors, contig_lengths,
            positional_slop=positional_slop,
            max_candidates=max_candidates,
            min_seeds=min_seeds,
        )
        proposal_dict = result.to_dict()
        rejected = proposal_dict["diagnostics"]["rejected_by_reason"]
        # Geometrically incompatible windows are legitimate retrieval misses.
        # Malformed coordinates or missing fields must not become silent misses.
        malformed = set(rejected) - {"reference_window_too_short", "no_possible_reference_overlap"}
        if malformed:
            raise ValueError(f"invalid anchor input for {name}: {', '.join(sorted(malformed))}")
        candidates = proposal_dict["proposals"]
        correct = any(
            candidate["chrom"] == origin["chrom"]
            and candidate["strand"] == origin["strand"]
            and candidate["start"] <= origin["start"]
            and candidate["end"] >= origin["start"] + item["reference_span"]
            for candidate in candidates for origin in item["origins"]
        )
        rows.append({
            "read_id": name,
            "stratum": item.get("stratum", "unclassified"),
            "has_origin": bool(item["origins"]),
            "candidate_contains_origin": correct,
            "candidate_count": len(candidates),
            "diagnostics": proposal_dict["diagnostics"],
        })
    report = _summarize(rows, "candidate_contains_origin")
    report.update({
        "schema": "genocache.retrieval-evaluation.v1",
        "truth_sha256": file_sha256(truth_file),
        "anchors_sha256": file_sha256(anchors_file),
        "authoritative_alignment": False,
        "metric": "a same-strand proposed interval contains the supplied true reference span",
        "config": {"positional_slop": positional_slop, "max_candidates": max_candidates, "min_seeds": min_seeds},
        "per_read": rows,
    })
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    mapped = sub.add_parser("alignments")
    mapped.add_argument("--alignment", type=Path, required=True)
    mapped.add_argument("--tolerance", type=int, default=100)
    retrieval = sub.add_parser("retrieval")
    retrieval.add_argument("--anchors", type=Path, required=True)
    retrieval.add_argument("--contigs", type=Path, required=True, help="JSON name:contig_length mapping")
    retrieval.add_argument("--positional-slop", type=int, default=256)
    retrieval.add_argument("--max-candidates", type=int, default=8)
    retrieval.add_argument("--min-seeds", type=int, default=2)
    for command in (mapped, retrieval):
        command.add_argument("--truth", type=Path, required=True)
        command.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == "alignments":
            report = evaluate_alignments(args.alignment, args.truth, tolerance=args.tolerance)
        else:
            report = evaluate_retrieval(
                args.anchors, args.truth, json.loads(args.contigs.read_text()),
                positional_slop=args.positional_slop,
                max_candidates=args.max_candidates, min_seeds=args.min_seeds,
            )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        fd, temporary = tempfile.mkstemp(prefix=".evaluation-", dir=args.output.parent)
        try:
            with os.fdopen(fd, "w") as stream:
                json.dump(report, stream, indent=2, allow_nan=False)
                stream.write("\n")
            os.replace(temporary, args.output)
        finally:
            Path(temporary).unlink(missing_ok=True)
        print(json.dumps({k: v for k, v in report.items() if k != "per_read"}, indent=2, allow_nan=False))
        return 0
    except (OSError, ValueError, RuntimeError) as exc:
        parser.exit(2, f"genocache-evaluate: {exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
