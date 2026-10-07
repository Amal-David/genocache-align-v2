"""Correctness contracts for non-authoritative embedding region proposals."""

import json
import math
import random
import unittest
from dataclasses import asdict, replace

from genocache.regions import Anchor, MAX_CANDIDATES, MAX_HITS, propose_regions


CONTIGS = {"chr1": 20_000, "chr2": 20_000}


def hit(qs=0, qe=100, *, origin=1000, strand="+", seed_id=None,
        read_id="read1", read_length=1000, chrom="chr1", similarity=0.9,
        ref_start=None, ref_end=None):
    oriented_query = qs if strand == "+" else read_length - qe
    start = origin + oriented_query if ref_start is None else ref_start
    end = start + qe - qs if ref_end is None else ref_end
    return Anchor(
        read_id=read_id, read_length=read_length,
        query_start=qs, query_end=qe,
        seed_id=seed_id if seed_id is not None else f"seed_{qs}_{qe}",
        chrom=chrom, ref_start=start, ref_end=end,
        strand=strand, similarity=similarity,
    )


def propose(anchors, **kwargs):
    options = {"positional_slop": 0, "min_coverage": 0}
    options.update(kwargs)
    return propose_regions(anchors, CONTIGS, **options)


class TestRegionProposals(unittest.TestCase):
    def test_forward_coordinates_are_derived_from_query_offsets(self):
        result = propose([hit(100, 200), hit(600, 700)])
        self.assertEqual(len(result.proposals), 1)
        region = result.proposals[0]
        self.assertEqual((region.chrom, region.strand, region.start, region.end),
                         ("chr1", "+", 1000, 2000))
        self.assertEqual((region.diag_low, region.diag_high), (1000, 1000))
        self.assertEqual(region.distinct_seeds, 2)
        self.assertEqual(region.query_coverage_bases, 200)
        self.assertAlmostEqual(region.query_coverage_fraction, 0.2)
        self.assertFalse(result.diagnostics["authoritative"])

    def test_negative_strand_uses_original_forward_read_offsets_once(self):
        result = propose([hit(100, 200, strand="-"), hit(600, 700, strand="-")])
        self.assertEqual(len(result.proposals), 1)
        region = result.proposals[0]
        self.assertEqual((region.strand, region.start, region.end), ("-", 1000, 2000))
        self.assertEqual(region.query_intervals, ((100, 200), (600, 700)))

    def test_positive_hits_do_not_create_negative_evidence(self):
        result = propose([hit(100, 200), hit(600, 700)])
        self.assertEqual({region.strand for region in result.proposals}, {"+"})

    def test_opposite_strands_and_contigs_never_pool_support(self):
        anchors = [
            hit(0, 100, strand="+"),
            hit(300, 400, strand="-"),
            hit(500, 600, chrom="chr2"),
        ]
        self.assertEqual(propose(anchors).proposals, ())

    def test_thirty_two_neighbors_from_one_seed_are_one_vote(self):
        anchors = [hit(origin=1000 + i * 32) for i in range(32)]
        self.assertEqual(propose(anchors, positional_slop=600).proposals, ())
        result = propose(anchors, positional_slop=600, min_seeds=1)
        self.assertTrue(result.proposals)
        for region in result.proposals:
            self.assertEqual(region.distinct_seeds, 1)
            self.assertEqual(region.query_coverage_bases, 100)
            self.assertEqual(region.score, 100)

    def test_seed_id_aliases_do_not_increase_support(self):
        anchors = [hit(seed_id=f"alias_{i}") for i in range(32)]
        result = propose(anchors)
        self.assertEqual(result.proposals, ())
        self.assertEqual(result.diagnostics["duplicate_hits_removed"], 31)
        self.assertEqual(result.diagnostics["anchors_accepted"], 1)

    def test_aliases_with_different_windows_still_count_one_query_interval(self):
        anchors = [hit(seed_id="a"), hit(seed_id="b", origin=1001)]
        result = propose(anchors, positional_slop=1, min_seeds=1)
        self.assertTrue(result.proposals)
        self.assertTrue(all(region.distinct_seeds == 1 for region in result.proposals))

    def test_overlapping_query_seeds_use_union_coverage(self):
        anchors = [hit(0, 512), hit(32, 544)]
        region = propose(anchors).proposals[0]
        self.assertEqual(region.distinct_seeds, 2)
        self.assertEqual(region.query_coverage_bases, 544)
        self.assertEqual(region.score, 544)
        self.assertAlmostEqual(region.query_coverage_fraction, 0.544)
        self.assertEqual(propose(anchors, min_coverage=0.75).proposals, ())

    def test_adjacent_query_intervals_can_cover_the_whole_read(self):
        region = propose([hit(0, 500), hit(500, 1000)], min_coverage=1).proposals[0]
        self.assertEqual(region.query_coverage_bases, 1000)
        self.assertEqual(region.query_coverage_fraction, 1)

    def test_separate_loci_on_same_contig_remain_separate(self):
        anchors = [hit(qs, qs + 100, origin=origin)
                   for origin in (1000, 6000) for qs in (0, 500)]
        result = propose(anchors)
        self.assertEqual([(region.start, region.end) for region in result.proposals],
                         [(1000, 2000), (6000, 7000)])
        self.assertFalse(result.diagnostics["truncated"])

    def test_transitive_proximity_does_not_create_three_way_support(self):
        anchors = [hit(i * 200, i * 200 + 100, origin=1000 + i * 900)
                   for i in range(3)]
        self.assertEqual(propose(anchors, positional_slop=500, min_seeds=3).proposals, ())
        result = propose(anchors, positional_slop=500, min_seeds=2)
        self.assertEqual(len(result.proposals), 2)
        self.assertTrue(all(region.distinct_seeds == 2 for region in result.proposals))
        self.assertTrue(all(region.diag_high - region.diag_low <= 100
                            for region in result.proposals))

    def test_ann_window_uncertainty_is_preserved(self):
        anchors = [
            hit(0, 100, ref_start=1000, ref_end=1600),
            hit(300, 400, ref_start=1450, ref_end=2150),
        ]
        region = propose(anchors).proposals[0]
        self.assertEqual((region.diag_low, region.diag_high), (1150, 1500))
        self.assertEqual((region.start, region.end), (1150, 2500))
        self.assertLess(region.diag_low, 1280)
        self.assertGreater(region.diag_high, 1280)

    def test_a_short_reference_window_requires_explicit_slop(self):
        anchors = [
            hit(0, 100, ref_start=1000, ref_end=1080),
            hit(300, 400, ref_start=1300, ref_end=1380),
        ]
        rejected = propose(anchors)
        self.assertEqual(rejected.proposals, ())
        self.assertEqual(rejected.diagnostics["rejected_by_reason"],
                         {"reference_window_too_short": 2})
        region = propose(anchors, positional_slop=10).proposals[0]
        self.assertEqual((region.diag_low, region.diag_high), (990, 990))

    def test_indel_slack_is_explicit_and_extends_the_candidate_window(self):
        region = propose([hit(0, 100), hit(300, 400)], indel_slack=20).proposals[0]
        self.assertEqual((region.diag_low, region.diag_high), (980, 1020))
        self.assertEqual((region.start, region.end), (980, 2040))
        with self.assertRaisesRegex(ValueError, "indel_slack"):
            propose([hit()], indel_slack=1001)

    def test_inclusive_diagonal_endpoints_share_a_single_possible_start(self):
        anchors = [
            hit(0, 100, ref_start=1000, ref_end=1101),
            hit(300, 400, ref_start=1301, ref_end=1402),
        ]
        region = propose(anchors).proposals[0]
        self.assertEqual((region.diag_low, region.diag_high), (1001, 1001))

    def test_adjacent_nonoverlapping_diagonals_do_not_share_support(self):
        anchors = [hit(0, 100, origin=1000), hit(300, 400, origin=1001)]
        self.assertEqual(propose(anchors).proposals, ())

    def test_invalid_hits_are_rejected_with_reason_counts(self):
        good = [hit(0, 100), hit(300, 400)]
        invalid = [
            replace(good[0], ref_start=-1),
            replace(good[0], ref_end=20_001),
            replace(good[0], query_end=1001),
            replace(good[0], similarity=math.nan),
            replace(good[0], similarity=math.inf),
            replace(good[0], chrom="absent"),
            replace(good[0], strand="?"),
            replace(good[0], seed_id=-1),
            replace(good[0], query_start=True),
            {},
            None,
        ]
        result = propose(good + invalid)
        self.assertEqual(len(result.proposals), 1)
        self.assertEqual(result.diagnostics["anchors_rejected"], len(invalid))
        self.assertEqual(result.diagnostics["rejected_by_reason"]["nonfinite_similarity"], 2)
        self.assertFalse(result.diagnostics["authoritative"])

    def test_input_budget_reports_truncation_and_consumes_only_one_lookahead(self):
        def stream():
            yield hit(0, 100)
            yield hit(300, 400)
            yield hit(600, 700)
            raise AssertionError("Consumed beyond max_hits plus one lookahead")

        result = propose(stream(), max_hits=2)
        self.assertEqual(len(result.proposals), 1)
        self.assertTrue(result.diagnostics["truncated"])
        self.assertEqual(result.diagnostics["truncation_reasons"], ["max_hits"])
        self.assertEqual(result.diagnostics["anchors_inspected"], 2)
        self.assertEqual(result.diagnostics["anchors_seen_at_least"], 3)

    def test_output_budget_reports_omitted_loci(self):
        anchors = [hit(qs, qs + 100, origin=origin)
                   for origin in (1000, 6000) for qs in (0, 500)]
        result = propose(anchors, max_candidates=1)
        self.assertEqual(len(result.proposals), 1)
        self.assertTrue(result.diagnostics["truncated"])
        self.assertEqual(result.diagnostics["truncation_reasons"], ["max_candidates"])
        self.assertEqual(result.diagnostics["budget_omitted_sweep_intervals"], 1)

    def test_shuffle_is_deterministic_within_budget(self):
        anchors = [hit(qs, qs + 100, origin=origin, similarity=sim)
                   for origin, sim in ((1000, 0.9), (6000, 0.95))
                   for qs in (0, 300, 600)]
        anchors += [replace(anchors[0], seed_id="alias"), replace(anchors[3], seed_id="alias2")]
        expected = propose(anchors, positional_slop=15).to_dict()
        rng = random.Random(42)
        for _ in range(10):
            rng.shuffle(anchors)
            self.assertEqual(propose(anchors, positional_slop=15).to_dict(), expected)

    def test_reference_windows_are_clipped_but_diagonal_uncertainty_remains(self):
        anchors = [hit(0, 100, origin=0), hit(200, 300, origin=0)]
        result = propose_regions(anchors, {"chr1": 600}, positional_slop=10, min_coverage=0)
        region = result.proposals[0]
        self.assertEqual((region.start, region.end), (0, 600))
        self.assertEqual((region.diag_low, region.diag_high), (-10, 10))
        self.assertTrue(region.clipped_to_contig)

    def test_mixed_read_identity_is_an_error(self):
        with self.assertRaisesRegex(ValueError, "exactly one"):
            propose([hit(), hit(300, 400, read_id="other")])
        with self.assertRaisesRegex(ValueError, "exactly one"):
            propose([hit(), hit(300, 400, read_length=2000)])

    def test_plain_mappings_are_accepted_and_output_is_json_serializable(self):
        result = propose([asdict(hit(0, 100)), asdict(hit(300, 400))])
        decoded = json.loads(json.dumps(result.to_dict(), allow_nan=False))
        self.assertEqual(decoded["proposals"][0]["query_intervals"], [[0, 100], [300, 400]])
        self.assertFalse(decoded["diagnostics"]["authoritative"])

    def test_empty_input_is_non_authoritative_and_not_truncated(self):
        result = propose([])
        self.assertEqual(result.proposals, ())
        self.assertEqual(result.diagnostics["anchors_inspected"], 0)
        self.assertFalse(result.diagnostics["truncated"])
        self.assertFalse(result.diagnostics["authoritative"])

    def test_finite_extreme_similarities_cannot_overflow_rank_metadata(self):
        anchors = [hit(0, 100, similarity=1.7e308), hit(300, 400, similarity=1.7e308)]
        region = propose(anchors).proposals[0]
        self.assertTrue(math.isfinite(region.mean_similarity))
        self.assertEqual(region.score, 200)
        json.dumps(region.to_dict(), allow_nan=False)

    def test_sweep_best_coverage_matches_an_independent_small_reference_oracle(self):
        # Enumerate possible integer starts and actual query bases. This does not
        # reproduce the endpoint sweep, heaps, compressed-coordinate union tree,
        # or candidate reconstruction used by the implementation.
        rng = random.Random(719)
        for trial in range(80):
            anchors = []
            for i in range(12):
                qs = rng.randrange(0, 70)
                span = rng.randrange(5, 25)
                origin = rng.randrange(30, 90)
                strand = rng.choice(("+", "-"))
                oriented_query = qs if strand == "+" else 100 - qs - span
                ref_start = origin + oriented_query
                anchors.append(hit(
                    qs, qs + span, read_length=100, seed_id=str(i), strand=strand,
                    ref_start=ref_start, ref_end=ref_start + span + rng.randrange(0, 15),
                    similarity=rng.uniform(-1, 1),
                ))
            anchors.append(replace(anchors[0], seed_id="same_interval_alias"))
            slop = rng.randrange(0, 6)
            oracle_best_coverage = 0
            for strand in ("+", "-"):
                for origin in range(0, 130):
                    supported_intervals = set()
                    for anchor in anchors:
                        if anchor.strand != strand:
                            continue
                        q = (anchor.query_start if strand == "+"
                             else anchor.read_length - anchor.query_end)
                        span = anchor.query_end - anchor.query_start
                        if (anchor.ref_start - slop <= origin + q
                                and origin + q + span <= anchor.ref_end + slop):
                            supported_intervals.add((anchor.query_start, anchor.query_end))
                    if len(supported_intervals) < 2:
                        continue
                    query_bases = {base for qs, qe in supported_intervals for base in range(qs, qe)}
                    oracle_best_coverage = max(oracle_best_coverage, len(query_bases))
            result = propose(anchors, positional_slop=slop)
            observed = result.proposals[0].query_coverage_bases if result.proposals else 0
            with self.subTest(trial=trial):
                self.assertEqual(observed, oracle_best_coverage)

    def test_invalid_parameters_fail_explicitly(self):
        options = [
            {"positional_slop": -1}, {"indel_slack": -1},
            {"max_hits": MAX_HITS + 1}, {"max_candidates": MAX_CANDIDATES + 1},
            {"max_hits": 0}, {"max_candidates": 0}, {"min_seeds": 0},
            {"min_coverage": -0.1}, {"min_coverage": 1.1},
            {"min_coverage": math.nan}, {"min_coverage": True},
            {"positional_slop": 1.5},
        ]
        for option in options:
            with self.subTest(option=option), self.assertRaises(ValueError):
                propose([], **option)


if __name__ == "__main__":
    unittest.main()
