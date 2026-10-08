"""Executable evidence of historical defects; these are NOT acceptance tests.

The assertions deliberately describe known broken behavior in the archived code.
A green result means that the audit reproduced the defect, not that the legacy
aligner is safe. Run with pytest, or with the Python standard library:

    python -m unittest discover -s production/tests -p test_legacy_contracts.py -v

The repository backup does not include model weights, indexes, or its old runtime.
To keep this audit independent of PyTorch, FAISS, BioPython and native WFA, the
loader executes unchanged class/function AST definitions directly from archived
source in an isolated namespace. Imports and CLI entrypoints are not executed.
Only the WFA dependency is substituted, with an explicit negative-cost test
double. No entries are added to sys.modules and no historical file is modified.
These are component-contract reproductions, not an end-to-end alignment run.
"""

import ast
import contextlib
import csv
import io
import re
import types
import unittest
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple


ROOT = Path(__file__).resolve().parents[2]
V4 = Path("genocache-v4.1-production")
CORE = V4 / "genocache_core"
NAL = V4 / "development/training/nal_aligned"
CHR22 = V4 / "development/training/nal_chr22_PRODUCTION_BACKUP/alignment"


def archived_definitions(relative_path, *names, **bindings):
    """Load actual archived definitions without importing the old dependency tree."""
    path = ROOT / relative_path
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    nodes = [
        node for node in tree.body
        if isinstance(node, (ast.ClassDef, ast.FunctionDef)) and node.name in names
    ]
    found = {node.name for node in nodes}
    if found != set(names):
        raise AssertionError(f"Missing archived definitions in {path}: {set(names) - found}")
    namespace = {
        "__name__": "isolated_legacy_contract_audit",
        "__file__": str(path),
        "Dict": Dict,
        "List": List,
        "Optional": Optional,
        "Tuple": Tuple,
        "Path": Path,
        "defaultdict": defaultdict,
        "re": re,
    }
    namespace.update(bindings)
    selected = ast.Module(body=nodes, type_ignores=[])
    exec(compile(selected, str(path), "exec"), namespace)
    return namespace


class NegativeCostWfaDouble:
    """Only a score-contract fixture; this does not implement sequence alignment."""

    def __init__(self, *args, **kwargs):
        self.constructor_args = args
        self.constructor_kwargs = kwargs
        self.score = 0
        self.cigarstring = ""

    def wavefront_align(self, query, reference):
        self.score = 0 if query == reference else -200
        self.cigarstring = f"{len(query)}M"
        return self.score


def archived_fast_aligner(genome, mode="semi-global"):
    pywfa_double = types.ModuleType("isolated_pywfa_double")
    pywfa_double.WavefrontAligner = NegativeCostWfaDouble
    namespace = archived_definitions(
        CORE / "fast_alignment.py", "FastAligner", pywfa=pywfa_double
    )
    with contextlib.redirect_stdout(io.StringIO()):
        return namespace["FastAligner"](genome, mode=mode)


def archived_extend(aligner):
    namespace = archived_definitions(
        CORE / "extend_phase.py", "ExtendPhase", FastAligner=type(aligner)
    )
    return namespace["ExtendPhase"](aligner)


def nearby_hits_from_one_seed():
    return [
        {
            "seed_idx": 0,
            "seed_pos": 0,
            "ref_chr": "chr1",
            "ref_pos": 10000 + i * 32,
            "strand": "+",
            "similarity": 0.9,
            "rank": i,
        }
        for i in range(32)
    ]


class TestLegacyContractEvidence(unittest.TestCase):
    """Passing tests reproduce archived defects; never use as a release gate."""

    def test_evidence_extend_selects_larger_wfa_cost_over_perfect_match(self):
        aligner = archived_fast_aligner({"correct": "A" * 100, "wrong": "C" * 100})
        candidates = [
            {"chr": name, "start": 0, "end": 100, "score": 3, "num_seeds": 3}
            for name in ("correct", "wrong")
        ]
        result = archived_extend(aligner).extend_and_score("A" * 100, candidates)
        self.assertEqual(result["chr"], "wrong")
        self.assertEqual(result["alignment_score"], 200)

    def test_evidence_extend_rejects_the_only_perfect_zero_cost_alignment(self):
        aligner = archived_fast_aligner({"correct": "A" * 100})
        candidate = {"chr": "correct", "start": 0, "end": 100, "score": 3, "num_seeds": 3}
        result = archived_extend(aligner).extend_and_score("A" * 100, [candidate])
        self.assertIsNone(result)

    def test_evidence_requested_semiglobal_mode_never_reaches_wfa_backend(self):
        for mode in ("global", "semi-global"):
            aligner = archived_fast_aligner({"chr1": "A" * 100}, mode=mode)
            self.assertEqual(aligner.mode, mode)
            self.assertEqual(aligner.wfa_aligner.constructor_args, ())
            self.assertEqual(aligner.wfa_aligner.constructor_kwargs, {})

    def test_evidence_reference_coordinates_are_always_the_full_search_window(self):
        aligner = archived_fast_aligner({"chr1": "C" * 10 + "A" * 100 + "C" * 10})
        result = aligner.align("A" * 100, "C" * 10 + "A" * 100 + "C" * 10)
        self.assertEqual(result["start_ref"], 0)
        self.assertEqual(result["end_ref"], 119)
        # The only exact full-query occurrence starts at 10, not at the returned 0.
        self.assertEqual(("C" * 10 + "A" * 100 + "C" * 10).find("A" * 100), 10)

    def test_evidence_thirty_two_neighbors_count_as_thirty_two_supporting_seeds(self):
        chainer_type = archived_definitions(CHR22 / "chaining_nal.py", "NALChaining")["NALChaining"]
        chains = chainer_type().chain_anchors(nearby_hits_from_one_seed(), read_len=2000)
        self.assertEqual(chains[0]["score"], 32)
        self.assertEqual(chains[0]["num_seeds"], 1)

    def test_evidence_positive_only_anchors_create_chains_on_both_strands(self):
        chainer_type = archived_definitions(CHR22 / "chaining_nal.py", "NALChaining")["NALChaining"]
        chains = chainer_type().chain_anchors(nearby_hits_from_one_seed(), read_len=2000)
        self.assertEqual({chain["strand"] for chain in chains}, {"+", "-"})
        self.assertTrue(all(anchor["strand"] == "+" for anchor in nearby_hits_from_one_seed()))

    def test_evidence_missing_wfa_still_produces_a_mapped_record_with_mapq_sixty(self):
        nal_type = archived_definitions(NAL / "align_nal.py", "NALAligner", WFA_AVAILABLE=False)["NALAligner"]
        aligner = nal_type.__new__(nal_type)
        chain = {"ref_chr": "chr1", "ref_pos": 25, "strand": "+", "score": 6}
        aligner.seeder = types.SimpleNamespace(get_anchors=lambda *args, **kwargs: [])
        aligner.chainer = types.SimpleNamespace(chain_with_rescue_check=lambda *args, **kwargs: ([chain], False))
        aligner.seed_len = 512
        aligner.K = 32
        aligner.wfa_aligner = None
        aligner.load_reference_region = lambda *args: "T" * 100
        result = aligner.align_read("read1", "A" * 100)
        self.assertTrue(result["mapped"])
        self.assertEqual(result["cigar"], "100M")
        self.assertEqual(result["alignment_score"], 0)
        self.assertEqual(result["mapq"], 60)

    def test_evidence_missing_reference_returns_fabricated_full_length_cigar(self):
        nal_type = archived_definitions(NAL / "align_nal.py", "NALAligner", WFA_AVAILABLE=False)["NALAligner"]
        aligner = nal_type.__new__(nal_type)
        aligner.wfa_aligner = None

        def unavailable_reference(*args):
            raise FileNotFoundError("controlled audit fixture: reference missing")

        aligner.load_reference_region = unavailable_reference
        chain = {"ref_chr": "chr1", "ref_pos": 25, "strand": "+", "score": 6}
        with contextlib.redirect_stdout(io.StringIO()):
            result = aligner._align_with_wfa("A" * 100, chain, 100)
        self.assertEqual(result, {"cigar": "100M", "score": 0})

    def test_evidence_nm_ignores_all_mismatches_inside_m_cigar_operations(self):
        parse_nm = archived_definitions(CORE / "sam_output.py", "parse_cigar_for_nm")["parse_cigar_for_nm"]
        self.assertEqual(parse_nm("4M", "AAAA", "TTTT"), 0)
        self.assertEqual(sum(a != b for a, b in zip("AAAA", "TTTT")), 4)

    def test_evidence_sam_formatter_invents_high_quality_scores(self):
        namespace = archived_definitions(
            CORE / "sam_output.py",
            "parse_cigar_for_nm", "calculate_mapq", "calculate_alignment_identity",
            "format_sam_tags", "format_sam_line",
        )
        record = namespace["format_sam_line"](
            "read1", "ACGT",
            {"chr": "chr1", "start": 0, "cigar": "4M", "alignment_score": 4, "num_seeds": 6},
            {"chr1": 100},
        )
        self.assertEqual(record.split("\t")[10], "IIII")

    def test_evidence_tier_summary_omits_cost_of_earlier_rescue_attempts(self):
        tier_type = archived_definitions(NAL / "test/align_adaptive_3tier.py", "ThreeTierAligner")["ThreeTierAligner"]
        aligner = tier_type.__new__(tier_type)
        result = {
            "mapped": True,
            "tier_used": 2,
            "timing": {"tier1": 0.1, "tier2": 0.2, "tier3": 0},
        }
        aligner.align_read = lambda *args: result
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            aligner.align_reads([("read1", "ACGT")])
        self.assertIn("Weighted average: 200.0ms/read", output.getvalue())
        self.assertAlmostEqual(sum(result["timing"].values()) * 1000, 300)

    def test_evidence_checked_in_improvement_did_not_increase_chromosome_correctness(self):
        observed = {}
        for name in ("baseline_results.txt", "improved_results.txt"):
            with (ROOT / V4 / "validation/results" / name).open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            observed[name] = (
                len(rows),
                sum(row["correct"] == "True" for row in rows),
                sum(row["predicted_chr"] != "unmapped" for row in rows),
            )
        self.assertEqual(observed["baseline_results.txt"], (100, 26, 78))
        self.assertEqual(observed["improved_results.txt"], (100, 26, 74))


if __name__ == "__main__":
    unittest.main()
