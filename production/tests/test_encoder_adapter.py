"""Coordinate, streaming, and artifact tests; real Torch execution is optional."""

from __future__ import annotations

import argparse
import gzip
import importlib.util
import json
from pathlib import Path
import tempfile
import types
import unittest
from unittest import mock


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "encode_torchscript.py"
SPEC = importlib.util.spec_from_file_location("genocache_encode_torchscript", SCRIPT)
encoder = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(encoder)
HAS_NUMPY = importlib.util.find_spec("numpy") is not None
HAS_TORCH = importlib.util.find_spec("torch") is not None


class EncoderHelperTests(unittest.TestCase):
    def test_uniform_distinct_full_seed_windows(self):
        self.assertEqual(encoder.uniform_seed_starts(511), ())
        self.assertEqual(encoder.uniform_seed_starts(512), (0,))
        self.assertEqual(encoder.uniform_seed_starts(513), (0, 1))
        self.assertEqual(encoder.uniform_seed_starts(15, 4, 4), (0, 3, 7, 11))
        self.assertEqual(encoder.uniform_seed_starts(15, 4, 1), (5,))
        self.assertEqual(encoder.full_window_count(9, 4, 3), 2)
        self.assertEqual(encoder.full_window_count(2, 4, 3), 0)
        for invalid in (0, -1, True, 1.2):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                encoder.uniform_seed_starts(20, window=invalid)

    def test_reverse_complement_iupac(self):
        original = b"ACGTRYSWKMBDHVN"
        self.assertEqual(encoder.reverse_complement(original), b"NBDHVKMWSRYACGT")
        self.assertEqual(encoder.reverse_complement(encoder.reverse_complement(original)), original)
        self.assertEqual(encoder.reverse_complement(b"aaccg"), b"CGGTT")
        with self.assertRaises(ValueError):
            encoder.reverse_complement(b"ACU")

    def test_both_orientations_keep_original_forward_coordinates(self):
        rows = list(encoder.seed_windows("read-1", b"AACCGTTAG", window=4, seeds=2))
        self.assertEqual([sequence for sequence, _ in rows], [b"AACC", b"GGTT", b"TTAG", b"CTAA"])
        plus, minus = rows[2][1], rows[3][1]
        self.assertEqual((plus["query_start"], plus["query_end"]), (5, 9))
        self.assertEqual((minus["query_start"], minus["query_end"]), (5, 9))
        self.assertEqual(minus["read_length"], 9)
        self.assertEqual(plus["seed_id"], minus["seed_id"])
        self.assertEqual((plus["query_strand"], minus["query_strand"]), ("+", "-"))

    def test_full_reference_windows_are_chunk_partition_invariant(self):
        contigs = [("chrA", b"AACCGTTAGCG"), ("tiny", b"A"), ("chrB", b"TTACCGGA")]
        for window in (1, 3, 8):
            for stride in (1, 3, 12):
                expected = [(sequence[start:start + window], {"chrom": name, "start": start,
                             "end": start + window, "strand": "+"})
                            for name, sequence in contigs
                            for start in range(0, len(sequence) - window + 1, stride)]
                for chunk in (1, 2, 8, 30):
                    events = []
                    for name, sequence in contigs:
                        events.append(("start", name))
                        events.extend(("bases", sequence[i:i + chunk]) for i in range(0, len(sequence), chunk))
                        events.append(("end", len(sequence)))
                    with self.subTest(window=window, stride=stride, chunk=chunk):
                        self.assertEqual(list(encoder.reference_windows(events, window, stride)), expected)

    def test_batching_preserves_order_and_flushes_final_partial_batch(self):
        rows = [(bytes([65 + i]), {"row": i}) for i in range(7)]
        groups = list(encoder.batches(iter(rows), 3))
        self.assertEqual([len(group) for group in groups], [3, 3, 1])
        self.assertEqual([item for group in groups for item in group], rows)
        self.assertEqual(list(encoder.batches([], 3)), [])
        with self.assertRaises(ValueError):
            list(encoder.batches(rows, 0))

    def test_short_read_skip_error_and_limit(self):
        events = [("start", "r"), ("bases", b"AAC"), ("end", 3)]
        self.assertEqual(list(encoder.read_windows(events, window=4, seeds=8, max_read_length=9)), [])
        with self.assertRaisesRegex(encoder.ValidationError, "shorter"):
            list(encoder.read_windows(events, window=4, seeds=8, max_read_length=9, short_reads="error"))
        with self.assertRaisesRegex(encoder.ValidationError, "max-read-length"):
            list(encoder.read_windows(events, window=4, seeds=8, max_read_length=2))

    def test_missing_torch_is_explicit(self):
        with mock.patch.object(encoder.importlib, "import_module", side_effect=ImportError):
            with self.assertRaisesRegex(RuntimeError, "Torch is required"):
                encoder.TorchEncoder(Path("unused.pt"), "cpu", 1, 17)

    def test_requested_cuda_does_not_silently_fall_back(self):
        fake_torch = types.SimpleNamespace(cuda=types.SimpleNamespace(is_available=lambda: False))
        with mock.patch.object(encoder.importlib, "import_module", return_value=fake_torch):
            with self.assertRaisesRegex(RuntimeError, "CUDA was requested"):
                encoder.TorchEncoder(Path("unused.pt"), "cuda", 1, 17)


class EncoderParserTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.addCleanup(self.temporary.cleanup)

    def put(self, data: bytes, *, compressed: bool = False) -> Path:
        path = self.root / ("input.gz" if compressed else "input.fa")
        path.write_bytes(gzip.compress(data) if compressed else data)
        return path

    def test_wrapped_crlf_fasta_and_gzip(self):
        data = b">chr1 description\r\naAcc\r\ngtN\r\n>chr2\r\nTA\r\n"
        expected = [("start", "chr1"), ("bases", b"AACC"), ("bases", b"GTN"),
                    ("end", 7), ("start", "chr2"), ("bases", b"TA"), ("end", 2)]
        for compressed in (False, True):
            self.assertEqual(list(encoder.sequence_events(self.put(data, compressed=compressed), reference=True)), expected)

    def test_wrapped_fastq_quality_header_characters_and_no_final_newline(self):
        events = list(encoder.sequence_events(self.put(b"@r\nAC\nGT\n+r\n!!\n@@\n@s\nAA\n+\n##")))
        self.assertEqual(events, [("start", "r"), ("bases", b"AC"), ("bases", b"GT"), ("end", 4),
                                  ("start", "s"), ("bases", b"AA"), ("end", 2)])

    def test_malformed_sequences_fail(self):
        for data in (b"", b">r\n", b">r\nAC U\n", b"@r\nAC\n+\n!\n", b"@r\nAC\n+other\n!!\n",
                     b"@r\nAC\n+\n!!!\n", b"@r\nAC\n+\n \x7f\n"):
            with self.subTest(data=data), self.assertRaises(encoder.ValidationError):
                list(encoder.sequence_events(self.put(data)))
        with self.assertRaisesRegex(encoder.ValidationError, "Reference must be FASTA"):
            list(encoder.sequence_events(self.put(b"@r\nAC\n+\n!!\n"), reference=True))

    def test_read_limit_checked_during_sequence_stream(self):
        events = encoder.sequence_events(self.put(b">r\nAACCGG\n"), max_length=4)
        self.assertEqual(next(events), ("start", "r"))
        with self.assertRaisesRegex(encoder.ValidationError, "max-read-length"):
            next(events)

    def test_unwrapped_reference_is_streamed_in_bounded_pieces(self):
        from genocache.io import CHUNK_BYTES
        length = CHUNK_BYTES * 2 + 31
        path = self.put(b">r\n" + b"A" * length)
        pieces = (value for kind, value in encoder.sequence_events(path, reference=True) if kind == "bases")
        sizes = [len(piece) for piece in pieces]
        self.assertGreater(len(sizes), 1)
        self.assertLessEqual(max(sizes), CHUNK_BYTES + 1)
        self.assertEqual(sum(sizes), length)

    def config(self, mode="reads"):
        return {"mode": mode, "window": 4, "stride": 3 if mode == "reference" else None,
                "seeds": 8 if mode == "reads" else None, "max_read_length": 100,
                "short_reads": "skip", "ambiguous": "zero"}

    def test_census_counts_short_sequences_and_orientations(self):
        stats = encoder.census(self.put(b">r1\nACGTN\n>short\nAA\n"), self.config(), self.root / "names.sqlite")
        self.assertEqual(stats["sequence_count"], 2)
        self.assertEqual(stats["base_count"], 7)
        self.assertEqual(stats["skipped_short_sequences"], 1)
        self.assertEqual(stats["seed_intervals"], 2)
        self.assertEqual(stats["vector_count"], 4)
        self.assertEqual(stats["ambiguous_bases"], 1)

    def test_census_rejects_duplicate_first_token_names(self):
        with self.assertRaisesRegex(encoder.ValidationError, "Duplicate sequence identifier"):
            encoder.census(self.put(b">same a\nACGT\n>same b\nACGT\n"), self.config(), self.root / "names.sqlite")


@unittest.skipUnless(HAS_NUMPY, "NumPy is not installed; array/artifact tests skipped")
class EncoderArrayTests(unittest.TestCase):
    setUp = EncoderParserTests.setUp
    put = EncoderParserTests.put

    def test_one_hot_channel_order_and_ambiguity(self):
        import numpy as np
        values = encoder.one_hot_batch([b"aCgT", b"NNRY"], 4)
        self.assertEqual(values.dtype, np.float32)
        self.assertEqual(values.shape, (2, 4, 4))
        np.testing.assert_array_equal(values[0], np.eye(4, dtype=np.float32))
        np.testing.assert_array_equal(values[1], np.zeros((4, 4), dtype=np.float32))
        with self.assertRaisesRegex(ValueError, "ambiguous"):
            encoder.one_hot_batch([b"ACGN"], 4, "error")
        with self.assertRaisesRegex(ValueError, "exactly"):
            encoder.one_hot_batch([b"ACG"], 4)

    def test_output_numeric_contract_rejects_invalid_batches(self):
        import numpy as np
        valid = encoder.validate_vectors(np.array([[1.0, 2.0]], dtype=np.float64), 1)
        self.assertEqual(valid.dtype, np.float32)
        for values in (np.array([[0.0, 0.0]]), np.array([[1e-20, 0.0]]),
                       np.array([[1e30, 0.0]]), np.array([[1e300, 0.0]]),
                       np.array([[float("nan"), 1.0]]), np.array([[float("inf"), 1.0]]),
                       np.array([[1, 2]]), np.ones((2, 4)), np.ones((1, 0)), np.ones(4)):
            with self.subTest(shape=values.shape, dtype=values.dtype), self.assertRaises(ValueError):
                encoder.validate_vectors(values, 1)

    def arguments(self, output: str) -> argparse.Namespace:
        checkpoint = self.root / "encoder.fixture"
        checkpoint.write_bytes(b"NOT A REAL CHECKPOINT - runtime is explicitly stubbed in this artifact test")
        return argparse.Namespace(mode="reads", input=self.put(b">r\nAACCGTTA\n>short\nAA\n"),
            encoder=checkpoint, output=self.root / output, window=4, seeds=3, max_read_length=100,
            short_reads="skip", ambiguous="zero", batch_size=4, device="cpu", threads=1, seed=17)

    def test_stub_runtime_artifacts_are_batched_ordered_and_content_identified(self):
        import numpy as np
        seen_batches = []

        class StubRuntime:
            def __init__(self, *_args):
                self.runtime = {"device": "stub", "note": "unit-test fixture; no Torch/model execution"}

            def encode(self, sequences, window, ambiguous):
                seen_batches.append(len(sequences))
                return encoder.one_hot_batch(sequences, window, ambiguous).mean(axis=2) + np.float32(1)

        args = self.arguments("encoded-a")
        with mock.patch.object(encoder, "TorchEncoder", StubRuntime):
            first = encoder.encode_file(args)
            args.output = self.root / "encoded-b"
            second = encoder.encode_file(args)
        self.assertEqual(seen_batches, [4, 2, 4, 2])
        self.assertEqual(first["output_id"], second["output_id"])
        self.assertEqual(first["recipe_id"], second["recipe_id"])
        self.assertEqual(first["counts"]["vector_count"], 6)
        self.assertEqual(first["counts"]["skipped_short_sequences"], 1)
        self.assertEqual(first["reads_sha256"], encoder.hash_file(args.input)["sha256"])
        self.assertEqual(first["encoder_sha256"], encoder.hash_file(args.encoder)["sha256"])
        self.assertGreaterEqual(first["timing"]["full_elapsed_seconds"], first["timing"]["encoding_seconds"])
        path = self.root / "encoded-a"
        values = np.load(path / "vectors.npy", mmap_mode="r", allow_pickle=False)
        metadata = [json.loads(line) for line in (path / "queries.jsonl").read_text().splitlines()]
        self.assertEqual(values.shape, (len(metadata), 4))
        self.assertEqual([row["query_start"] for row in metadata], [0, 0, 2, 2, 4, 4])
        self.assertEqual([row["query_strand"] for row in metadata], ["+", "-"] * 3)
        self.assertEqual(set(item.name for item in path.iterdir()), {"vectors.npy", "queries.jsonl", "manifest.json"})

    def test_runtime_failure_leaves_no_published_or_temporary_artifacts(self):
        class FailingRuntime:
            def __init__(self, *_args):
                self.runtime = {"device": "stub"}

            def encode(self, *_args):
                raise ValueError("deliberate model-contract failure")

        args = self.arguments("failed")
        with mock.patch.object(encoder, "TorchEncoder", FailingRuntime):
            with self.assertRaisesRegex(ValueError, "model-contract"):
                encoder.encode_file(args)
        self.assertFalse(args.output.exists())
        self.assertFalse(list(self.root.glob(".encode-*")))
        self.assertFalse(list(self.root.glob(".encoder-inputs-*")))

    def test_existing_output_is_rejected(self):
        args = self.arguments("existing")
        args.output.mkdir()
        with self.assertRaises(FileExistsError):
            encoder.encode_file(args)


@unittest.skipUnless(HAS_TORCH and HAS_NUMPY, "Torch is not installed; real TorchScript execution is unverified")
class OptionalTorchExecutionTests(unittest.TestCase):
    def test_cpu_torchscript_contract_on_test_only_pooling_module(self):
        import numpy as np
        import torch

        # This is a contract fixture, not trained weights or a genomic accuracy test.
        with tempfile.TemporaryDirectory() as work:
            artifact = Path(work) / "test-contract.pt"
            model = torch.nn.Sequential(torch.nn.AdaptiveAvgPool1d(1), torch.nn.Flatten(1)).eval()
            torch.jit.trace(model, torch.ones(2, 4, 4)).save(str(artifact))
            runtime = encoder.TorchEncoder(artifact, "cpu", 1, 17)
            values = runtime.encode([b"AACC", b"GGTT", b"ACGT"], 4, "zero")
            np.testing.assert_allclose(values, [[.5, .5, 0, 0], [0, 0, .5, .5], [.25] * 4])
            with self.assertRaisesRegex(ValueError, "zero"):
                runtime.encode([b"NNNN"], 4, "zero")


if __name__ == "__main__":
    unittest.main()
