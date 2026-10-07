"""Numerical contract tests for the analytical model, not speed benchmarks."""

import json
import math
import unittest

from genocache.model import (
    CacheInputs, IndexInputs, ModelInputs, Observation, RoutingInputs,
    StageTiming, analyze_model, wilson_interval,
)


class ModelTests(unittest.TestCase):
    def test_unknown_is_not_zero(self):
        result = analyze_model({})
        self.assertIsNone(result["index"]["rows"])
        self.assertIsNone(result["routing"]["speedup"])
        self.assertIsNone(result["error_budget"]["total_error_probability"])
        self.assertIsNone(result["cache"]["speedup_including_build"])
        self.assertEqual(result["observations"], {})
        json.dumps(result, allow_nan=False)

    def test_dataclass_and_mapping_have_identical_results(self):
        actual = analyze_model(ModelInputs(
            routing=RoutingInputs(front_end_fraction=0.2, fallback_time_fraction=0.1),
            observations={"errors": Observation(successes=0, trials=100, kind="error")},
        ))
        expected = analyze_model({
            "routing": {"front_end_fraction": 0.2, "fallback_time_fraction": 0.1},
            "observations": {"errors": {"successes": 0, "trials": 100, "kind": "error"}},
        })
        self.assertEqual(actual, expected)

    def test_human_genome_payload_and_dense_cost(self):
        result = analyze_model({"index": {
            "genome_bases": 3_100_000_000, "stride_bases": 1000,
            "orientations": 1, "embedding_dim": 256, "dtype": "float16",
            "id_bytes": 8, "metadata_bytes_per_row": 16,
            "query_vectors_per_read": 20, "batch_query_vectors": 32,
            "score_dtype": "float32", "hnsw_m": 32, "graph_link_bytes": 4,
        }})["index"]
        self.assertEqual(result["rows"], 3_100_000)
        self.assertEqual(result["memory_bytes"]["vector_payload"], 1_587_200_000)
        self.assertEqual(result["memory_bytes"]["hnsw_base_layer_links"], 793_600_000)
        self.assertEqual(result["memory_bytes"]["dense_batch_score_matrix"], 396_800_000)
        self.assertEqual(result["flops"]["dense_exact_per_query_vector"], 1_587_200_000)
        self.assertEqual(result["flops"]["dense_exact_per_read"], 31_744_000_000)

    def test_contig_boundaries_orientations_and_clipped_tail(self):
        result = analyze_model(ModelInputs(index=IndexInputs(
            contig_lengths=[1001, 1001, 0], stride_bases=1000,
            orientations=2, embedding_dim=4, dtype="int8",
        )))["index"]
        self.assertEqual(result["genome_bases"], 2002)
        self.assertEqual(result["rows"], 8)
        self.assertEqual(result["memory_bytes"]["vector_payload"], 32)
        self.assertIsNone(result["memory_bytes"]["raw_payload_ids_metadata_subtotal"])

    def test_containment_includes_integer_boundary(self):
        config = {"window_bases": 1250, "max_query_reference_span": 250, "stride_bases": 1001}
        self.assertTrue(analyze_model({"index": config})["index"]["full_chunk_containment_guaranteed_by_tiling"])
        config["stride_bases"] = 1002
        self.assertFalse(analyze_model({"index": config})["index"]["full_chunk_containment_guaranteed_by_tiling"])

    def test_pq_bookkeeping_and_ivf_cost(self):
        result = analyze_model({"index": {
            "genome_bases": 3_100_000_000, "stride_bases": 1000,
            "orientations": 1, "embedding_dim": 256, "dtype": "float16",
            "id_bytes": 8, "metadata_bytes_per_row": 16,
            "pq_subquantizers": 32, "pq_bits": 8, "pq_codebook_dtype": "float32",
            "pq_retain_full_vectors": False, "ivf_lists": 4096, "ivf_nprobe": 32,
            "ivf_centroid_dtype": "float32", "ivf_list_header_bytes": 64,
        }})["index"]
        self.assertEqual(result["memory_bytes"]["pq_code_bytes_per_row"], 32)
        self.assertEqual(result["memory_bytes"]["pq_codebooks"], 262_144)
        self.assertEqual(result["memory_bytes"]["pq_payload_bookkeeping_subtotal"], 173_862_144)
        self.assertEqual(result["memory_bytes"]["ivf_pq_payload_bookkeeping_subtotal"], 178_318_592)
        self.assertEqual(result["flops"]["ivf_balanced_list_rows_per_query"], 24_218.75)
        self.assertEqual(result["flops"]["ivf_flat_scan_and_centroids_per_query"], 14_497_152)

    def test_pq_bit_packing_rounds_up_per_row(self):
        result = analyze_model({"index": {
            "genome_bases": 100, "stride_bases": 10, "orientations": 1,
            "embedding_dim": 12, "pq_subquantizers": 3, "pq_bits": 6,
        }})["index"]
        self.assertEqual(result["memory_bytes"]["pq_code_bytes_per_row"], 3)
        self.assertEqual(result["memory_bytes"]["pq_codes"], 30)
        self.assertIsNone(result["memory_bytes"]["pq_payload_bookkeeping_subtotal"])

    def test_staged_amdahl_and_overhead(self):
        result = analyze_model(ModelInputs(stages=[
            StageTiming("seed", 20, 1), StageTiming("chain", 40, 4),
            StageTiming("align", 40, 4),
        ], new_overhead_seconds=10))["staged_amdahl"]
        self.assertEqual(result["projected_stage_seconds"], 40)
        self.assertEqual(result["stage_only_speedup"], 2.5)
        self.assertEqual(result["speedup_including_overhead"], 2)
        self.assertAlmostEqual(result["idealized_selected_stage_speedup_limit"], 5)

    def test_missing_stage_acceleration_is_unknown(self):
        result = analyze_model({"stages": [{"name": "chain", "baseline_seconds": 10}]})["staged_amdahl"]
        self.assertIsNone(result["projected_stage_seconds"])
        self.assertIsNone(result["stage_only_speedup"])

    def test_missing_overhead_does_not_imply_zero(self):
        result = analyze_model({"stages": [{"name": "chain", "baseline_seconds": 10, "speedup": 2}]})["staged_amdahl"]
        self.assertEqual(result["stage_only_speedup"], 2)
        self.assertIsNone(result["speedup_including_overhead"])

    def test_eliminating_all_work_is_not_infinity_json(self):
        result = analyze_model({"stages": [{"name": "work", "baseline_seconds": 10, "eliminated": True}], "new_overhead_seconds": 0})
        self.assertTrue(result["staged_amdahl"]["idealized_limit_unbounded"])
        self.assertIsNone(result["staged_amdahl"]["stage_only_speedup"])
        json.dumps(result, allow_nan=False)

    def test_routing_uses_time_not_read_fraction(self):
        result = analyze_model({"routing": {
            "baseline_seconds": 100, "front_end_seconds": 10,
            "fallback_baseline_seconds": 75, "fallback_read_fraction": 0.05,
            "target_speedup": 2,
        }})["routing"]
        self.assertAlmostEqual(result["speedup"], 1 / 0.85)
        self.assertEqual(result["fallback_time_fraction"], 0.75)
        self.assertFalse(result["target_met_in_model"])
        self.assertAlmostEqual(result["target_maximum_fallback_time_fraction"], 0.4)

    def test_routing_fraction_only_and_unknown_read_time(self):
        result = analyze_model({"routing": {"front_end_fraction": 0.2, "fallback_time_fraction": 0.1}})["routing"]
        self.assertAlmostEqual(result["speedup"], 10 / 3)
        self.assertIsNone(result["projected_seconds"])
        incomplete = analyze_model({"routing": {"front_end_fraction": 0.2, "fallback_read_fraction": 0.1}})["routing"]
        self.assertIsNone(incomplete["speedup"])

    def test_conditional_pipeline_and_error_budget(self):
        result = analyze_model({"error_budget": {
            "candidate_recall": 0.99, "chain_survival_given_candidate": 0.98,
            "verification_success_given_chain": 0.97, "reporting_correct_given_verified": 0.96,
            "accepted_fraction": 0.8, "accepted_error_probability": 0.001,
            "fallback_error_probability": 0.002, "maximum_total_error_probability": 0.001,
            "acceptance_given_candidate_miss": 0.05, "maximum_miss_error_probability": 0.0001,
        }})["error_budget"]
        self.assertAlmostEqual(result["pipeline_correct_probability"], 0.90345024)
        self.assertAlmostEqual(result["total_error_probability"], 0.0012)
        self.assertFalse(result["total_error_budget_met"])
        self.assertAlmostEqual(result["maximum_accepted_error_probability"], 0.00075)
        self.assertAlmostEqual(result["candidate_miss_error_contribution"], 0.0005)
        self.assertAlmostEqual(result["minimum_candidate_recall_for_miss_budget"], 0.998)

    def test_no_acceptance_requires_only_fallback_error(self):
        result = analyze_model({"error_budget": {"accepted_fraction": 0, "fallback_error_probability": 0.02}})["error_budget"]
        self.assertEqual(result["total_error_probability"], 0.02)
        self.assertIsNone(result["maximum_accepted_error_probability"])

    def test_impossible_error_budget_is_explicit(self):
        result = analyze_model({"error_budget": {
            "accepted_fraction": 0.5, "fallback_error_probability": 0.1,
            "maximum_total_error_probability": 0.01,
        }})["error_budget"]
        self.assertTrue(result["fallback_alone_exceeds_budget"])
        self.assertIsNone(result["maximum_accepted_error_probability"])

    def test_union_recall_does_not_assume_independence(self):
        result = analyze_model({"union_recall": {"recall_a": 0.99, "recall_b": 0.99}})["union_recall"]
        self.assertEqual(result["recall_lower_bound"], 0.99)
        self.assertEqual(result["recall_upper_bound"], 1)
        self.assertIsNone(result["observed_union_recall"])
        measured = analyze_model({"union_recall": {"recall_a": 0.9, "recall_b": 0.8, "observed_union_recall": 0.95}})["union_recall"]
        self.assertAlmostEqual(measured["observed_joint_miss_probability"], 0.05)

    def test_cache_break_even_and_amortization(self):
        result = analyze_model(ModelInputs(cache=CacheInputs(
            baseline_seconds_per_request=10, lookup_seconds_per_request=0.1,
            hit_seconds_per_request=1, miss_seconds_per_request=11, hit_rate=0.8,
            build_seconds=100, amortization_requests=100,
        )))["cache"]
        self.assertAlmostEqual(result["warm_seconds_per_request"], 3.1)
        self.assertAlmostEqual(result["seconds_per_request_including_build"], 4.1)
        self.assertAlmostEqual(result["speedup_including_build"], 10 / 4.1)
        self.assertAlmostEqual(result["strict_speedup_hit_rate_boundary"], 0.21)
        self.assertEqual(result["strict_speedup_requires_hit_rate"], ">")

    def test_cache_hit_bounds_have_distinct_assumptions(self):
        result = analyze_model({"cache": {
            "observed_requests": 100, "distinct_keys": 4, "capacity_entries": 2,
            "uniform_key_space": 100, "request_counts_by_key": [60, 20, 15, 5],
        }})["cache"]
        self.assertEqual(result["cold_trace_hit_rate_upper_bound"], 0.96)
        self.assertEqual(result["uniform_stationary_hit_rate_upper_bound"], 0.02)
        self.assertEqual(result["optimal_static_histogram_hit_rate_upper_bound"], 0.8)

    def test_cache_warm_cost_does_not_hide_unknown_build(self):
        result = analyze_model({"cache": {
            "baseline_seconds_per_request": 1, "lookup_seconds_per_request": 0.01,
            "hit_seconds_per_request": 0.1, "miss_seconds_per_request": 1,
            "hit_rate": 0.5,
        }})["cache"]
        self.assertAlmostEqual(result["warm_seconds_per_request"], 0.56)
        self.assertIsNone(result["seconds_per_request_including_build"])
        self.assertIsNone(result["strict_speedup_hit_rate_boundary"])

    def test_expensive_cache_hit_reverses_threshold(self):
        result = analyze_model({"cache": {
            "baseline_seconds_per_request": 5, "lookup_seconds_per_request": 0,
            "hit_seconds_per_request": 10, "miss_seconds_per_request": 1, "build_seconds": 0,
        }})["cache"]
        self.assertEqual(result["strict_speedup_requires_hit_rate"], "<")
        self.assertAlmostEqual(result["strict_speedup_hit_rate_boundary"], 4 / 9)

    def test_wilson_known_reference_values(self):
        zero_errors = wilson_interval(0, 100)
        self.assertEqual(zero_errors["estimate"], 0)
        self.assertAlmostEqual(zero_errors["lower"], 0)
        self.assertAlmostEqual(zero_errors["upper"], 0.0369934982, places=9)
        half = wilson_interval(50, 100)
        self.assertAlmostEqual(half["lower"], 0.4038315304, places=9)
        self.assertAlmostEqual(half["upper"], 0.5961684696, places=9)
        perfect_recall = wilson_interval(100, 100)
        self.assertAlmostEqual(perfect_recall["lower"], 1 - zero_errors["upper"])

    def test_empty_observation_has_unknown_rate(self):
        self.assertIsNone(wilson_interval(0, 0)["estimate"])
        result = analyze_model({"observations": {"errors": {"successes": 0, "kind": "error"}}})
        self.assertIsNone(result["observations"]["errors"]["upper"])

    def test_invalid_and_contradictory_inputs_fail(self):
        invalid = [
            {"unexpected": 1}, {"index": {"stride_bases": 0}},
            {"index": {"embedding_dim": True}}, {"index": {"embedding_dim": 1.5}},
            {"index": {"dtype": "fp17"}}, {"index": {"orientations": 3}},
            {"index": {"genome_bases": 10, "contig_lengths": [9]}},
            {"index": {"embedding_dim": 10, "pq_subquantizers": 3}},
            {"index": {"ivf_lists": 4, "ivf_nprobe": 5}},
            {"routing": {"front_end_fraction": -1}},
            {"routing": {"front_end_fraction": math.nan}},
            {"routing": {"front_end_fraction": math.inf}},
            {"routing": {"baseline_seconds": 100, "front_end_seconds": 20, "front_end_fraction": 0.1}},
            {"routing": {"baseline_seconds": 100, "fallback_baseline_seconds": 101}},
            {"union_recall": {"recall_a": 0.9, "recall_b": 0.8, "observed_union_recall": 0.7}},
            {"error_budget": {"candidate_recall": 1.1}},
            {"cache": {"observed_requests": 100, "distinct_keys": 101}},
            {"cache": {"request_counts_by_key": [1, 2], "observed_requests": 4}},
            {"observations": {"bad": {"successes": 101, "trials": 100}}},
            {"observations": {"bad": {"successes": 1, "trials": 1, "confidence": 1}}},
            {"stages": [{"name": "chain", "speedup": 0}]},
            {"stages": [{"name": "same"}, {"name": "same"}]},
            {"stages": [None]}, {"stages": [{}]},
        ]
        for config in invalid:
            with self.subTest(config=config), self.assertRaises(ValueError):
                analyze_model(config)


if __name__ == "__main__":
    unittest.main()
