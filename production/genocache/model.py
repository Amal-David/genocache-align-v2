"""Transparent, dependency-free cost and error models for genome mapping.

Inputs are supplied measurements or explicit assumptions, never measurements made
by this module. Missing quantities remain ``None``. All durations must use the
same workload and timing boundary. See ``docs/RESEARCH_MODEL.md`` for definitions.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, fields
from statistics import NormalDist
from typing import Any


_DTYPE_BYTES = {"float64": 8, "float32": 4, "float16": 2, "bfloat16": 2,
                "int8": 1, "uint8": 1}


@dataclass(frozen=True)
class IndexInputs:
    genome_bases: int | None = None
    contig_lengths: Sequence[int] | None = None
    stride_bases: int | None = None
    window_bases: int | None = None
    max_query_reference_span: int | None = None
    orientations: int | None = None
    embedding_dim: int | None = None
    dtype: str | None = None
    id_bytes: int | None = None
    metadata_bytes_per_row: int | None = None
    hnsw_m: int | None = None
    graph_link_bytes: int | None = None
    pq_subquantizers: int | None = None
    pq_bits: int | None = None
    pq_codebook_dtype: str | None = None
    pq_retain_full_vectors: bool | None = None
    ivf_lists: int | None = None
    ivf_nprobe: int | None = None
    ivf_centroid_dtype: str | None = None
    ivf_list_header_bytes: int | None = None
    query_vectors_per_read: int | None = None
    batch_query_vectors: int | None = None
    score_dtype: str | None = None
    ann_distance_evaluations_per_query: int | None = None
    encoder_flops_per_vector: float | None = None


@dataclass(frozen=True)
class StageTiming:
    name: str
    baseline_seconds: float | None = None
    speedup: float | None = None
    eliminated: bool = False


@dataclass(frozen=True)
class RoutingInputs:
    baseline_seconds: float | None = None
    front_end_seconds: float | None = None
    front_end_fraction: float | None = None
    fallback_baseline_seconds: float | None = None
    fallback_time_fraction: float | None = None
    fallback_read_fraction: float | None = None
    target_speedup: float | None = None


@dataclass(frozen=True)
class ErrorBudgetInputs:
    candidate_recall: float | None = None
    chain_survival_given_candidate: float | None = None
    verification_success_given_chain: float | None = None
    reporting_correct_given_verified: float | None = None
    accepted_fraction: float | None = None
    accepted_error_probability: float | None = None
    fallback_error_probability: float | None = None
    maximum_total_error_probability: float | None = None
    acceptance_given_candidate_miss: float | None = None
    maximum_miss_error_probability: float | None = None


@dataclass(frozen=True)
class UnionRecallInputs:
    recall_a: float | None = None
    recall_b: float | None = None
    observed_union_recall: float | None = None


@dataclass(frozen=True)
class CacheInputs:
    baseline_seconds_per_request: float | None = None
    lookup_seconds_per_request: float | None = None
    hit_seconds_per_request: float | None = None
    miss_seconds_per_request: float | None = None
    hit_rate: float | None = None
    build_seconds: float | None = None
    amortization_requests: int | None = None
    observed_requests: int | None = None
    distinct_keys: int | None = None
    capacity_entries: int | None = None
    uniform_key_space: int | None = None
    request_counts_by_key: Sequence[int] | None = None


@dataclass(frozen=True)
class Observation:
    successes: int | None = None
    trials: int | None = None
    kind: str = "recall"
    confidence: float = 0.95


@dataclass(frozen=True)
class ModelInputs:
    index: IndexInputs | Mapping[str, Any] | None = None
    stages: Sequence[StageTiming | Mapping[str, Any]] | None = None
    new_overhead_seconds: float | None = None
    routing: RoutingInputs | Mapping[str, Any] | None = None
    error_budget: ErrorBudgetInputs | Mapping[str, Any] | None = None
    union_recall: UnionRecallInputs | Mapping[str, Any] | None = None
    cache: CacheInputs | Mapping[str, Any] | None = None
    observations: Mapping[str, Observation | Mapping[str, Any]] | None = None


def _coerce(cls: type, value: Any, name: str) -> Any:
    if value is None:
        try:
            return cls()
        except TypeError as exc:
            raise ValueError(f"{name} is missing required fields") from exc
    if isinstance(value, cls):
        return value
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a {cls.__name__} or mapping")
    unknown = set(value) - {f.name for f in fields(cls)}
    if unknown:
        raise ValueError(f"unknown {name} fields: {', '.join(sorted(map(str, unknown)))}")
    try:
        return cls(**value)
    except TypeError as exc:
        raise ValueError(f"invalid {name}: {exc}") from exc


def _number(value: Any, name: str, *, positive: bool = False,
            probability: bool = False) -> float | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite number")
    try:
        result = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be finite") from exc
    if not math.isfinite(result) or result < 0 or (positive and result == 0):
        raise ValueError(f"{name} must be finite and {'positive' if positive else 'nonnegative'}")
    if probability and result > 1:
        raise ValueError(f"{name} must be between zero and one")
    return result


def _integer(value: Any, name: str, *, positive: bool = False) -> int | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{name} must be an integer")
    _number(value, name, positive=positive)
    return value


def _dtype(value: Any, name: str) -> int | None:
    if value is None:
        return None
    if not isinstance(value, str) or value not in _DTYPE_BYTES:
        raise ValueError(f"{name} must be one of {', '.join(_DTYPE_BYTES)}")
    return _DTYPE_BYTES[value]


def _product(*values: Any) -> Any:
    if any(v is None for v in values):
        return None
    return math.prod(values)


def _sum_known(values: Sequence[Any]) -> Any:
    return None if any(v is None for v in values) else sum(values)


def _ratio(numerator: Any, denominator: Any) -> float | None:
    if numerator is None or denominator is None or denominator == 0:
        return None
    return numerator / denominator


def _validate_finite_tree(value: Any, path: str = "result") -> None:
    """Reject arithmetic overflow rather than returning nonstandard JSON NaN/Inf."""
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError(f"{path} overflowed; choose smaller, consistently scaled inputs")
    if isinstance(value, Mapping):
        for key, child in value.items():
            _validate_finite_tree(child, f"{path}.{key}")
    elif isinstance(value, (tuple, list)):
        for pos, child in enumerate(value):
            _validate_finite_tree(child, f"{path}[{pos}]")


def _index_model(raw: Any) -> dict[str, Any]:
    x = _coerce(IndexInputs, raw, "index")
    positive = ("stride_bases", "window_bases", "max_query_reference_span",
                "orientations", "embedding_dim", "hnsw_m", "graph_link_bytes",
                "pq_subquantizers", "pq_bits", "ivf_lists", "ivf_nprobe",
                "query_vectors_per_read", "batch_query_vectors")
    integers = ("genome_bases", "id_bytes", "metadata_bytes_per_row",
                "ivf_list_header_bytes", "ann_distance_evaluations_per_query")
    v = {name: _integer(getattr(x, name), f"index.{name}", positive=True)
         for name in positive}
    v.update({name: _integer(getattr(x, name), f"index.{name}") for name in integers})
    if v["orientations"] not in (None, 1, 2):
        raise ValueError("index.orientations must be one or two")
    lengths = None
    if x.contig_lengths is not None:
        if not isinstance(x.contig_lengths, Sequence) or isinstance(x.contig_lengths, (str, bytes)):
            raise ValueError("index.contig_lengths must be a sequence of integers")
        lengths = [_integer(n, "index.contig_lengths item") for n in x.contig_lengths]
        if any(n is None for n in lengths):
            raise ValueError("index.contig_lengths cannot contain unknown entries")
        if v["genome_bases"] is not None and sum(lengths) != v["genome_bases"]:
            raise ValueError("index.genome_bases must equal the sum of contig_lengths")
        v["genome_bases"] = sum(lengths)
    stride, orientations, dim = v["stride_bases"], v["orientations"], v["embedding_dim"]
    rows = None
    if stride is not None and orientations is not None:
        if lengths is not None:
            rows = orientations * sum((n + stride - 1) // stride for n in lengths)
        elif v["genome_bases"] is not None:
            rows = orientations * ((v["genome_bases"] + stride - 1) // stride)
    scalar_bytes = _dtype(x.dtype, "index.dtype")
    score_bytes = _dtype(x.score_dtype, "index.score_dtype")
    codebook_scalar_bytes = _dtype(x.pq_codebook_dtype, "index.pq_codebook_dtype")
    centroid_scalar_bytes = _dtype(x.ivf_centroid_dtype, "index.ivf_centroid_dtype")
    if x.pq_retain_full_vectors is not None and type(x.pq_retain_full_vectors) is not bool:
        raise ValueError("index.pq_retain_full_vectors must be a boolean")
    m, bits = v["pq_subquantizers"], v["pq_bits"]
    if bits is not None and bits > 16:
        raise ValueError("index.pq_bits must be at most 16")
    if m is not None and dim is not None and dim % m:
        raise ValueError("index.embedding_dim must be divisible by pq_subquantizers")
    nlist, nprobe = v["ivf_lists"], v["ivf_nprobe"]
    if nlist is not None and nprobe is not None and nprobe > nlist:
        raise ValueError("index.ivf_nprobe cannot exceed ivf_lists")
    payload = _product(rows, dim, scalar_bytes)
    ids = _product(rows, v["id_bytes"])
    metadata = _product(rows, v["metadata_bytes_per_row"])
    graph = _product(rows, 2, v["hnsw_m"], v["graph_link_bytes"])
    code_bytes = (m * bits + 7) // 8 if m is not None and bits is not None else None
    codes = _product(rows, code_bytes)
    codebook = _product(dim, (2 ** bits if bits is not None else None), codebook_scalar_bytes)
    centroids = _product(nlist, dim, centroid_scalar_bytes)
    headers = _product(nlist, v["ivf_list_header_bytes"])
    pq_retained = (payload if x.pq_retain_full_vectors else 0) if x.pq_retain_full_vectors is not None else None
    pq_subtotal = _sum_known([codes, codebook, ids, metadata, pq_retained])
    ivf_pq_subtotal = _sum_known([pq_subtotal, centroids, headers])
    flat_flops = _product(2, rows, dim)
    ann_flops = _product(2, v["ann_distance_evaluations_per_query"], dim)
    visited = rows * nprobe / nlist if rows is not None and nlist is not None and nprobe is not None else None
    ivf_scan_flops = _product(2, visited, dim)
    ivf_centroid_flops = _product(2, nlist, dim)
    encoder_flops = _number(x.encoder_flops_per_vector, "index.encoder_flops_per_vector")
    coverage = None
    if stride is not None and v["window_bases"] is not None and v["max_query_reference_span"] is not None:
        # Inclusive integer starting positions permit the +1 term.
        coverage = stride <= v["window_bases"] - v["max_query_reference_span"] + 1
    return {
        "status": "calculated" if rows is not None else "unknown",
        "genome_bases": v["genome_bases"], "rows": rows,
        "row_count_scope": "per_contig" if lengths is not None else "single_contig_approximation",
        "row_rule": "One right-clipped window at starts 0, stride, ... < contig_length, per orientation.",
        "full_chunk_containment_guaranteed_by_tiling": coverage,
        "memory_bytes": {
            "vector_payload": payload, "ids": ids, "metadata": metadata,
            "raw_payload_ids_metadata_subtotal": _sum_known([payload, ids, metadata]),
            "hnsw_base_layer_links": graph,
            "raw_with_hnsw_base_links_subtotal": _sum_known([payload, ids, metadata, graph]),
            "pq_code_bytes_per_row": code_bytes, "pq_codes": codes,
            "pq_codebooks": codebook, "pq_retained_full_vectors": pq_retained,
            "pq_payload_bookkeeping_subtotal": pq_subtotal,
            "ivf_centroids": centroids, "ivf_list_headers": headers,
            "ivf_pq_payload_bookkeeping_subtotal": ivf_pq_subtotal,
            "dense_batch_score_matrix": _product(v["batch_query_vectors"], rows, score_bytes),
        },
        "flops": {
            "dense_exact_per_query_vector": flat_flops,
            "dense_exact_per_read": _product(flat_flops, v["query_vectors_per_read"]),
            "ann_full_precision_distances_per_query": ann_flops,
            "ivf_balanced_list_rows_per_query": visited,
            "ivf_flat_scan_per_query": ivf_scan_flops,
            "ivf_flat_centroid_search_per_query": ivf_centroid_flops,
            "ivf_flat_scan_and_centroids_per_query": _sum_known([ivf_scan_flops, ivf_centroid_flops]),
            "reference_encoding": _product(rows, encoder_flops),
            "query_encoding_per_read": _product(v["query_vectors_per_read"], encoder_flops),
        },
        "limitations": [
            "Memory subtotals exclude allocator alignment, upper graph layers, model weights, build/search scratch, replication and framework overhead.",
            "Dense dot products count approximately 2*rows*dimension FLOPs; top-K, normalization, memory traffic and encoder overhead are separate.",
            "IVF row estimates assume balanced lists and are not recall guarantees; PQ lookup costs are not modeled as full-precision distances.",
            "Containment assumes the supplied span already includes indels and windows never cross contig boundaries; it is not a biological recall guarantee.",
        ],
    }


def _stage_model(raw: Any, overhead_raw: Any) -> dict[str, Any]:
    overhead = _number(overhead_raw, "new_overhead_seconds")
    if raw is not None and (not isinstance(raw, Sequence) or isinstance(raw, (str, bytes))):
        raise ValueError("stages must be a sequence")
    stages, names = [], set()
    for item in raw or []:
        s = _coerce(StageTiming, item, "stage")
        if not isinstance(s.name, str) or not s.name.strip() or s.name in names:
            raise ValueError("stage names must be nonempty and unique")
        names.add(s.name)
        baseline = _number(s.baseline_seconds, f"stage.{s.name}.baseline_seconds")
        factor = _number(s.speedup, f"stage.{s.name}.speedup", positive=True)
        if type(s.eliminated) is not bool or (s.eliminated and factor is not None):
            raise ValueError("stage.eliminated must be boolean and cannot accompany speedup")
        projected = (0.0 if s.eliminated else _ratio(baseline, factor)) if baseline is not None else None
        stages.append({"name": s.name, "baseline_seconds": baseline, "speedup": factor,
                       "eliminated": s.eliminated, "projected_seconds": projected,
                       "selected_for_acceleration": s.eliminated or (factor is not None and factor > 1)})
    baseline_total = _sum_known([s["baseline_seconds"] for s in stages]) if stages else None
    projected_total = _sum_known([s["projected_seconds"] for s in stages]) if stages else None
    total_with_overhead = _sum_known([projected_total, overhead])
    fraction = None
    if baseline_total is not None and baseline_total > 0:
        fraction = sum(s["baseline_seconds"] for s in stages if s["selected_for_acceleration"]) / baseline_total
    for s in stages:
        s["baseline_fraction"] = _ratio(s["baseline_seconds"], baseline_total)
    return {
        "status": "calculated" if projected_total is not None else "unknown",
        "stages": stages, "baseline_seconds": baseline_total,
        "projected_stage_seconds": projected_total, "new_overhead_seconds": overhead,
        "projected_seconds_including_overhead": total_with_overhead,
        "stage_only_speedup": _ratio(baseline_total, projected_total),
        "speedup_including_overhead": _ratio(baseline_total, total_with_overhead),
        "selected_baseline_fraction": fraction,
        "idealized_selected_stage_speedup_limit": _ratio(1.0, 1 - fraction) if fraction is not None else None,
        "idealized_limit_unbounded": fraction == 1 if fraction is not None else None,
        "limitations": ["Amdahl model assumes sequential additive stage times; overlapping pipelines need a critical-path/throughput model.",
                        "The idealized limit eliminates only the selected accelerated stages and holds other baseline work fixed; zero denominators have no finite ratio."],
    }


def _routing_model(raw: Any) -> dict[str, Any]:
    x = _coerce(RoutingInputs, raw, "routing")
    base = _number(x.baseline_seconds, "routing.baseline_seconds")
    front = _number(x.front_end_seconds, "routing.front_end_seconds")
    fallback = _number(x.fallback_baseline_seconds, "routing.fallback_baseline_seconds")
    alpha = _number(x.front_end_fraction, "routing.front_end_fraction")
    frac = _number(x.fallback_time_fraction, "routing.fallback_time_fraction", probability=True)
    read_frac = _number(x.fallback_read_fraction, "routing.fallback_read_fraction", probability=True)
    target = _number(x.target_speedup, "routing.target_speedup", positive=True)
    if base is not None and fallback is not None and fallback > base:
        raise ValueError("fallback_baseline_seconds is a subset of baseline_seconds and cannot exceed it")
    for name, explicit, duration in (("front_end_fraction", alpha, front), ("fallback_time_fraction", frac, fallback)):
        inferred = _ratio(duration, base)
        if explicit is not None and inferred is not None and not math.isclose(explicit, inferred, rel_tol=1e-9, abs_tol=1e-12):
            raise ValueError(f"routing.{name} disagrees with supplied durations")
    if alpha is None:
        alpha = _ratio(front, base)
    if frac is None:
        frac = _ratio(fallback, base)
    if front is None:
        front = _product(alpha, base)
    if fallback is None:
        fallback = _product(frac, base)
    normalized = _sum_known([alpha, frac])
    projected = _sum_known([front, fallback])
    limit = 1 / target - alpha if target is not None and alpha is not None else None
    return {
        "status": "calculated" if normalized is not None else "unknown",
        "baseline_seconds": base, "front_end_fraction": alpha,
        "fallback_time_fraction": frac, "fallback_read_fraction": read_frac,
        "normalized_runtime": normalized, "projected_seconds": projected,
        "speedup": _ratio(1.0, normalized),
        "zero_runtime_in_idealized_model": normalized == 0 if normalized is not None else None,
        "target_speedup": target, "target_maximum_fallback_time_fraction": limit,
        "target_possible_even_with_no_fallback": limit >= 0 if limit is not None else None,
        "target_met_in_model": normalized <= 1 / target if normalized is not None and target is not None else None,
        "limitations": [
            "front_end cost must include all new non-fallback work: embedding, ANN, candidate verification, routing, transfers and relevant output/merge overhead.",
            "fallback_time_fraction is sum of baseline work for fallback reads divided by total baseline work, not the fraction of reads.",
            "This additive counterfactual assumes the fallback engine and per-read baseline costs remain comparable; batching and concurrency can change them.",
        ],
    }


def _error_model(raw: Any) -> dict[str, Any]:
    x = _coerce(ErrorBudgetInputs, raw, "error_budget")
    v = {f.name: _number(getattr(x, f.name), f"error_budget.{f.name}", probability=True) for f in fields(x)}
    r, a = v["candidate_recall"], v["accepted_fraction"]
    chain = _product(r, v["chain_survival_given_candidate"],
                     v["verification_success_given_chain"], v["reporting_correct_given_verified"])
    fast_term = (0.0 if a == 0 else _product(a, v["accepted_error_probability"])) if a is not None else None
    fallback_term = (0.0 if a == 1 else _product(1 - a, v["fallback_error_probability"])) if a is not None else None
    error = _sum_known([fast_term, fallback_term])
    epsilon = v["maximum_total_error_probability"]
    remaining = epsilon - fallback_term if epsilon is not None and fallback_term is not None else None
    accepted_limit = min(1.0, remaining / a) if remaining is not None and remaining >= 0 and a is not None and a > 0 else None
    miss = 1 - r if r is not None else None
    miss_error = _product(miss, v["acceptance_given_candidate_miss"])
    miss_budget, miss_accept = v["maximum_miss_error_probability"], v["acceptance_given_candidate_miss"]
    recall_required = max(0.0, 1 - miss_budget / miss_accept) if miss_budget is not None and miss_accept is not None and miss_accept > 0 else None
    return {
        "status": "calculated" if any(z is not None for z in (chain, error, miss_error)) else "unknown",
        "candidate_miss_probability": miss,
        "pipeline_correct_probability": chain,
        "pipeline_failure_probability": 1 - chain if chain is not None else None,
        "accepted_error_contribution": fast_term, "fallback_error_contribution": fallback_term,
        "total_error_probability": error,
        "total_error_budget_met": error <= epsilon if error is not None and epsilon is not None else None,
        "budget_after_fallback_errors": remaining,
        "maximum_accepted_error_probability": accepted_limit,
        "fallback_alone_exceeds_budget": remaining < 0 if remaining is not None else None,
        "candidate_miss_error_contribution": miss_error,
        "candidate_miss_error_budget_met": miss_error <= miss_budget if miss_error is not None and miss_budget is not None else None,
        "minimum_candidate_recall_for_miss_budget": recall_required,
        "miss_acceptance_is_zero": miss_accept == 0 if miss_accept is not None else None,
        "limitations": [
            "The pipeline product uses conditional probabilities conditioned on all preceding successful stages; it assumes no independence.",
            "Candidate misses produce wrong accepted mappings if no acceptable locus exists in the returned set. DP verification alone cannot establish that another locus was not missed.",
            "The accepted/fallback error mixture uses read fractions and errors conditional on each route. The fallback subset can be harder than the overall population.",
            "Candidate-miss error is one contribution, not the total error. Do not add it to total_error_probability because events may overlap.",
        ],
    }


def _union_model(raw: Any) -> dict[str, Any]:
    x = _coerce(UnionRecallInputs, raw, "union_recall")
    a = _number(x.recall_a, "union_recall.recall_a", probability=True)
    b = _number(x.recall_b, "union_recall.recall_b", probability=True)
    observed = _number(x.observed_union_recall, "union_recall.observed_union_recall", probability=True)
    low, high = (max(a, b), min(1.0, a + b)) if a is not None and b is not None else (None, None)
    if observed is not None and low is not None and not (low - 1e-12 <= observed <= high + 1e-12):
        raise ValueError("observed_union_recall is inconsistent with recall_a and recall_b on the same population")
    return {
        "status": "calculated" if low is not None or observed is not None else "unknown",
        "recall_lower_bound": low, "recall_upper_bound": high,
        "miss_lower_bound": 1 - high if high is not None else None,
        "miss_upper_bound": 1 - low if low is not None else None,
        "observed_union_recall": observed,
        "observed_joint_miss_probability": 1 - observed if observed is not None else None,
        "limitations": ["Bounds require the same reads, reference and success definition. No independence of the two candidate sources is assumed.",
                        "Measure paired misses directly; combining two 99% recalls does not by itself imply 99.99% union recall."],
    }


def _cache_model(raw: Any) -> dict[str, Any]:
    x = _coerce(CacheInputs, raw, "cache")
    base = _number(x.baseline_seconds_per_request, "cache.baseline_seconds_per_request")
    lookup = _number(x.lookup_seconds_per_request, "cache.lookup_seconds_per_request")
    hit = _number(x.hit_seconds_per_request, "cache.hit_seconds_per_request")
    miss = _number(x.miss_seconds_per_request, "cache.miss_seconds_per_request")
    rate = _number(x.hit_rate, "cache.hit_rate", probability=True)
    build = _number(x.build_seconds, "cache.build_seconds")
    amort_n = _integer(x.amortization_requests, "cache.amortization_requests", positive=True)
    requests = _integer(x.observed_requests, "cache.observed_requests")
    distinct = _integer(x.distinct_keys, "cache.distinct_keys")
    capacity = _integer(x.capacity_entries, "cache.capacity_entries")
    universe = _integer(x.uniform_key_space, "cache.uniform_key_space", positive=True)
    if requests is not None and distinct is not None:
        if distinct > requests or (requests > 0 and distinct == 0):
            raise ValueError("cache.distinct_keys must be between one and observed_requests for a nonempty trace")
    cold_bound = (requests - distinct) / requests if requests and distinct is not None else None
    uniform_bound = min(1.0, capacity / universe) if capacity is not None and universe is not None else None
    optimal_bound = None
    if x.request_counts_by_key is not None:
        counts_raw = x.request_counts_by_key
        if not isinstance(counts_raw, Sequence) or isinstance(counts_raw, (str, bytes)):
            raise ValueError("cache.request_counts_by_key must be a sequence of integers")
        counts = [_integer(n, "cache.request_counts_by_key item") for n in counts_raw]
        if any(n is None for n in counts):
            raise ValueError("cache.request_counts_by_key cannot contain unknown entries")
        if requests is not None and sum(counts) != requests:
            raise ValueError("cache.request_counts_by_key must sum to observed_requests")
        if distinct is not None and sum(n > 0 for n in counts) != distinct:
            raise ValueError("cache.request_counts_by_key must agree with distinct_keys")
        if capacity is not None and sum(counts) > 0:
            optimal_bound = sum(sorted(counts, reverse=True)[:capacity]) / sum(counts)
    amortized_build = 0.0 if build == 0 else _ratio(build, amort_n)
    warm_cost = None
    if lookup is not None and rate is not None:
        hit_term = 0.0 if rate == 0 else _product(rate, hit)
        miss_term = 0.0 if rate == 1 else _product(1 - rate, miss)
        warm_cost = _sum_known([lookup, hit_term, miss_term])
    all_cost = _sum_known([warm_cost, amortized_build])
    boundary, operator, can_win, always_wins = None, None, None, None
    if all(n is not None for n in (base, lookup, hit, miss, amortized_build)):
        at_zero = lookup + miss + amortized_build
        at_one = lookup + hit + amortized_build
        saving = miss - hit
        can_win = min(at_zero, at_one) < base
        always_wins = max(at_zero, at_one) < base
        if saving:
            boundary = (at_zero - base) / saving
            operator = ">" if saving > 0 else "<"
        else:
            operator = "all" if at_zero < base else ("none" if at_zero > base else "equal_for_all")
    return {
        "status": "calculated" if any(z is not None for z in (warm_cost, cold_bound, uniform_bound, optimal_bound)) else "unknown",
        "warm_seconds_per_request": warm_cost, "amortized_build_seconds_per_request": amortized_build,
        "seconds_per_request_including_build": all_cost,
        "warm_speedup": _ratio(base, warm_cost), "speedup_including_build": _ratio(base, all_cost),
        "strict_speedup_hit_rate_boundary": boundary, "strict_speedup_requires_hit_rate": operator,
        "some_feasible_hit_rate_is_faster": can_win, "all_feasible_hit_rates_are_faster": always_wins,
        "cold_trace_hit_rate_upper_bound": cold_bound,
        "uniform_stationary_hit_rate_upper_bound": uniform_bound,
        "optimal_static_histogram_hit_rate_upper_bound": optimal_bound,
        "limitations": [
            "Lookup is charged for every request; hit and miss durations exclude lookup. Build cost is additional to the comparable baseline and needs an amortization horizon.",
            "The cold-trace bound (requests-distinct_keys)/requests assumes no initial entries; the uniform bound assumes equiprobable keys; the histogram bound is for an oracle static cache under the supplied distribution, not a universal bound for adaptive caches.",
            "Result-cache keys must include every input affecting alignment. Reuse of a reference index is distinct from reuse of an alignment result.",
        ],
    }


def wilson_interval(successes: int, trials: int, confidence: float = 0.95) -> dict[str, Any]:
    """Two-sided Wilson score interval; ``successes`` may count errors instead.

    A zero-trial observation has unknown probability, not a zero error rate.
    These binomial intervals do not account for correlated reads or sample drift.
    """
    successes = _integer(successes, "successes")
    trials = _integer(trials, "trials")
    if successes is None or trials is None or successes > trials:
        raise ValueError("successes and trials must be known, with successes <= trials")
    conf = _number(confidence, "confidence", positive=True, probability=True)
    if conf is None or conf == 1:
        raise ValueError("confidence must be strictly between zero and one")
    if trials == 0:
        return {"successes": successes, "trials": trials, "confidence": conf,
                "estimate": None, "lower": None, "upper": None}
    # Use the lower tail to avoid rounding (1 + confidence) / 2 to exactly one.
    z = -NormalDist().inv_cdf((1 - conf) / 2)
    p = successes / trials
    z2 = z * z
    denominator = 1 + z2 / trials
    center = (p + z2 / (2 * trials)) / denominator
    half_width = z * math.sqrt((p * (1 - p) + z2 / (4 * trials)) / trials) / denominator
    return {"successes": successes, "trials": trials, "confidence": conf,
            "estimate": p, "lower": max(0.0, center - half_width),
            "upper": min(1.0, center + half_width)}


def _observation_model(raw: Any) -> dict[str, Any]:
    if raw is None:
        return {}
    if not isinstance(raw, Mapping):
        raise ValueError("observations must be a mapping of names to observations")
    result = {}
    for name, item in raw.items():
        if not isinstance(name, str) or not name.strip():
            raise ValueError("observation names must be nonempty strings")
        x = _coerce(Observation, item, f"observations.{name}")
        if x.kind not in ("recall", "error"):
            raise ValueError("observation.kind must be 'recall' or 'error'")
        successes = _integer(x.successes, f"observations.{name}.successes")
        trials = _integer(x.trials, f"observations.{name}.trials")
        confidence = _number(x.confidence, f"observations.{name}.confidence", positive=True, probability=True)
        if confidence is None or confidence == 1:
            raise ValueError("observation.confidence must be strictly between zero and one")
        if successes is None or trials is None:
            interval = {"successes": successes, "trials": trials, "confidence": confidence,
                        "estimate": None, "lower": None, "upper": None}
        else:
            interval = wilson_interval(successes, trials, confidence)
        result[name] = {"kind": x.kind, **interval,
                        "interpretation": "Interval for the supplied counted event; zero observed errors do not prove zero population error."}
    return result


def analyze_model(inputs: ModelInputs | Mapping[str, Any]) -> dict[str, Any]:
    """Return JSON-ready stage, search, routing, reliability and cache estimates.

    Unknown keys, nonfinite/negative inputs and contradictory measurements raise
    ``ValueError``. Timing speedups are mathematical predictions, not benchmarks.
    """
    if inputs is None:
        raise ValueError("inputs must be a ModelInputs or mapping; use {} for unknown inputs")
    x = _coerce(ModelInputs, inputs, "model")
    result = {
        "schema_version": 1,
        "evidence": "Calculated from supplied measurements or assumptions; this function does not run a benchmark.",
        "unknown_representation": None,
        "index": _index_model(x.index),
        "staged_amdahl": _stage_model(x.stages, x.new_overhead_seconds),
        "routing": _routing_model(x.routing),
        "error_budget": _error_model(x.error_budget),
        "union_recall": _union_model(x.union_recall),
        "cache": _cache_model(x.cache),
        "observations": _observation_model(x.observations),
    }
    _validate_finite_tree(result)
    return result
