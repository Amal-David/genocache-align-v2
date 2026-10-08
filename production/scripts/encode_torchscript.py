#!/usr/bin/env python3
"""Stream a trusted, trained TorchScript DNA encoder into GenoCache vector files.

Contract: float32 [batch, 4, window] in A,C,G,T channel order -> finite,
nonzero floating [batch, dimension]. No model downloads or training occur.
Torch and NumPy are imported only by execution/array helpers, never for --help.
"""

from __future__ import annotations

import argparse
import importlib
import os
from pathlib import Path
import re
import shutil
import sqlite3
import sys
import tempfile
import time
from typing import Any, Iterable, Iterator

# Permit running the checked-out script without an editable installation.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from genocache.io import (  # noqa: E402
    DNA_ALPHABET,
    QUALITY_ALPHABET,
    EngineError,
    ValidationError,
    _header_name,
    _line_pieces,
    _sequence_stream,
    canonical_bytes,
    content_id,
    fsync_directory,
    hash_file,
    snapshot_file,
    write_json,
)


SCHEMA = "genocache.encoder-output.v1"
COMPLEMENT = bytes.maketrans(b"ACGTRYSWKMBDHVN", b"TGCAYRSWMKVHDBN")
CHANNELS = bytes("ACGT".find(chr(i).upper()) if chr(i).upper() in "ACGT" else 4
                 for i in range(256))


def positive_int(value: Any, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return value


def full_window_count(length: int, window: int, stride: int) -> int:
    positive_int(window, "window")
    positive_int(stride, "stride")
    if isinstance(length, bool) or not isinstance(length, int) or length < 0:
        raise ValueError("length must be a nonnegative integer")
    return max(0, (length - window) // stride + 1)


def uniform_seed_starts(length: int, window: int = 512, seeds: int = 8) -> tuple[int, ...]:
    """Distinct, deterministic full windows; one requested seed uses the midpoint."""
    positive_int(window, "window")
    positive_int(seeds, "seeds")
    if isinstance(length, bool) or not isinstance(length, int) or length < 0:
        raise ValueError("length must be a nonnegative integer")
    if length < window:
        return ()
    last = length - window
    count = min(seeds, last + 1)
    if count == 1:
        return (last // 2,)
    return tuple(i * last // (count - 1) for i in range(count))


def reverse_complement(sequence: bytes) -> bytes:
    if sequence.translate(None, DNA_ALPHABET):
        raise ValueError("sequence contains a non-IUPAC DNA base")
    return sequence.upper().translate(COMPLEMENT)[::-1]


def sequence_events(
    path: Path, *, reference: bool = False, max_length: int | None = None
) -> Iterator[tuple[str, Any]]:
    """Validated start/name, bases/bytes, end/length events with bounded line reads.

    Qualities are checked and discarded. Reference chromosomes are never loaded
    into memory. Consumers must exhaust the iterator to validate the full input.
    """
    if max_length is not None:
        positive_int(max_length, "max_length")
    sequence_format, state, current_name = None, "header", None
    length = quality_length = sequence_count = 0

    def check_length() -> None:
        if not length:
            raise ValidationError("Empty sequences are not supported")
        if length >= 2**31:
            raise ValidationError("Individual sequences must be shorter than 2^31 bases")
        if max_length is not None and length > max_length:
            raise ValidationError(f"Read exceeds --max-read-length ({max_length} bases)")

    with _sequence_stream(Path(path)) as (stream, _compressed):
        for value, at_start, at_end, line in _line_pieces(stream, None):
            if sequence_format is None:
                if not value:
                    continue
                if at_start and value.startswith(b">"):
                    sequence_format = "fasta"
                elif at_start and value.startswith(b"@") and not reference:
                    sequence_format = "fastq"
                else:
                    raise ValidationError("Reference must be FASTA" if reference else "Reads must be FASTA or FASTQ")
            if sequence_format == "fasta":
                if at_start and value.startswith(b">"):
                    if current_name is not None:
                        check_length()
                        yield "end", length
                    current_name = _header_name(value, at_end, line, reference)
                    sequence_count += 1
                    length = 0
                    yield "start", current_name
                else:
                    if value.translate(None, DNA_ALPHABET):
                        raise ValidationError(f"Non-IUPAC DNA character at line {line}")
                    length += len(value)
                    if value:
                        check_length()
                        yield "bases", value.upper()
                continue
            if state == "header":
                if not value:
                    continue
                if not at_start or not value.startswith(b"@"):
                    raise ValidationError(f"Expected a FASTQ header at line {line}")
                current_name = _header_name(value, at_end, line, False)
                sequence_count += 1
                length = quality_length = 0
                state = "sequence"
                yield "start", current_name
            elif state == "sequence":
                if at_start and value.startswith(b"+"):
                    if not at_end:
                        raise ValidationError(f"FASTQ separator is too long at line {line}")
                    check_length()
                    repeated = value[1:].split(None, 1)
                    if repeated and repeated[0] != current_name.encode("ascii"):
                        raise ValidationError(f"FASTQ separator identifier mismatch at line {line}")
                    state = "quality"
                else:
                    if value.translate(None, DNA_ALPHABET):
                        raise ValidationError(f"Non-IUPAC DNA character at line {line}")
                    length += len(value)
                    if value:
                        check_length()
                        yield "bases", value.upper()
            else:
                if value.translate(None, QUALITY_ALPHABET):
                    raise ValidationError(f"Invalid FASTQ quality character at line {line}")
                quality_length += len(value)
                if quality_length > length:
                    raise ValidationError(f"FASTQ sequence/quality length mismatch at line {line}")
                if quality_length == length and at_end:
                    check_length()
                    yield "end", length
                    state = "header"
    if not sequence_count:
        raise ValidationError("Input contains no sequences")
    if sequence_format == "fasta":
        check_length()
        yield "end", length
    elif state == "quality" and quality_length == length:
        check_length()
        yield "end", length
    elif state != "header":
        raise ValidationError("Truncated FASTQ record or sequence/quality length mismatch")


def reference_windows(
    events: Iterable[tuple[str, Any]], window: int, stride: int
) -> Iterator[tuple[bytes, dict]]:
    """Full windows on a contig-local grid; RAM is one line chunk plus one window."""
    positive_int(window, "window")
    positive_int(stride, "stride")
    buffer = bytearray()
    chrom = None
    consumed = next_start = buffer_start = 0
    for kind, value in events:
        if kind == "start":
            chrom = value
            buffer.clear()
            consumed = next_start = buffer_start = 0
        elif kind == "bases":
            previous = consumed
            consumed += len(value)
            skip = min(len(value), max(0, next_start - previous)) if not buffer else 0
            if not buffer:
                buffer_start = previous + skip
            buffer.extend(value[skip:])
            while next_start + window <= consumed:
                offset = next_start - buffer_start
                seed = bytes(buffer[offset : offset + window])
                if len(seed) != window:
                    raise ValueError("reference window buffer inconsistency")
                yield seed, {"chrom": chrom, "start": next_start, "end": next_start + window, "strand": "+"}
                next_start += stride
            discard = min(len(buffer), max(0, next_start - buffer_start))
            del buffer[:discard]
            buffer_start += discard
        elif kind == "end":
            buffer.clear()


def seed_windows(
    read_id: str, sequence: bytes | bytearray, *, window: int = 512, seeds: int = 8,
    read_ordinal: int = 0,
) -> Iterator[tuple[bytes, dict]]:
    """Both orientations share the same original forward coordinates and seed_id."""
    for start in uniform_seed_starts(len(sequence), window, seeds):
        end = start + window
        seed = bytes(sequence[start:end]).upper()
        metadata = {
            "read_id": read_id, "read_length": len(sequence), "read_ordinal": read_ordinal,
            "query_start": start, "query_end": end, "seed_id": f"{start}:{end}",
        }
        yield seed, {**metadata, "query_strand": "+"}
        yield reverse_complement(seed), {**metadata, "query_strand": "-"}


def read_windows(
    events: Iterable[tuple[str, Any]], *, window: int, seeds: int, max_read_length: int,
    short_reads: str = "skip",
) -> Iterator[tuple[bytes, dict]]:
    positive_int(max_read_length, "max_read_length")
    if short_reads not in {"skip", "error"}:
        raise ValueError("short_reads must be skip or error")
    buffer = bytearray()
    name, ordinal = None, -1
    for kind, value in events:
        if kind == "start":
            buffer.clear()
            name = value
            ordinal += 1
        elif kind == "bases":
            if len(buffer) + len(value) > max_read_length:
                raise ValidationError(f"Read exceeds --max-read-length ({max_read_length} bases)")
            buffer.extend(value)
        elif kind == "end":
            if len(buffer) != value:
                raise ValueError("read length differs from parser length")
            if value < window and short_reads == "error":
                raise ValidationError("Read is shorter than --window; no padding is defined")
            yield from seed_windows(name, buffer, window=window, seeds=seeds, read_ordinal=ordinal)
            buffer.clear()


def batches(rows: Iterable[tuple[bytes, dict]], batch_size: int) -> Iterator[list[tuple[bytes, dict]]]:
    positive_int(batch_size, "batch_size")
    batch: list[tuple[bytes, dict]] = []
    for row in rows:
        batch.append(row)
        if len(batch) == batch_size:
            yield batch
            batch = []
    if batch:
        yield batch


def one_hot_batch(sequences: list[bytes], window: int, ambiguous: str = "zero"):
    """Batched float32 A,C,G,T, with explicit all-zero IUPAC ambiguity handling."""
    import numpy as np

    positive_int(window, "window")
    if ambiguous not in {"zero", "error"}:
        raise ValueError("ambiguous must be zero or error")
    if not sequences or any(len(sequence) != window for sequence in sequences):
        raise ValueError("every input sequence must have exactly --window bases")
    joined = b"".join(sequences).upper()
    if joined.translate(None, DNA_ALPHABET):
        raise ValueError("sequence contains a non-IUPAC DNA base")
    codes = np.frombuffer(CHANNELS, dtype=np.uint8)[np.frombuffer(joined, dtype=np.uint8)]
    if ambiguous == "error" and (codes == 4).any():
        raise ValueError("ambiguous bases are forbidden by --ambiguous error")
    codes = codes.reshape(len(sequences), window)
    return (codes[:, None, :] == np.arange(4, dtype=np.uint8)[None, :, None]).astype(np.float32)


def validate_vectors(values, expected_rows: int):
    """Apply the numeric contract to a whole CPU output batch."""
    import numpy as np

    values = np.asarray(values)
    if values.dtype.kind != "f":
        raise ValueError("Encoder output must be floating point")
    if values.ndim != 2 or values.shape[0] != expected_rows or values.shape[1] < 1:
        raise ValueError("Encoder output must have shape [batch, positive_dimension]")
    with np.errstate(over="ignore", invalid="ignore"):
        values = np.ascontiguousarray(values, dtype=np.float32)
        norms = np.linalg.norm(values, axis=1)
    if not np.isfinite(values).all():
        raise ValueError("Encoder produced nonfinite float32 vectors")
    if not np.isfinite(norms).all() or (norms <= 1e-12).any():
        raise ValueError("Encoder produced zero, near-zero, or overflowing vectors")
    return values


def census(path: Path, config: dict, database: Path) -> dict:
    """Count output rows and reject duplicate names using a bounded-memory disk DB."""
    reference = config["mode"] == "reference"
    counts = {"sequence_count": 0, "base_count": 0, "max_sequence_length": 0,
              "skipped_short_sequences": 0, "seed_intervals": 0, "vector_count": 0,
              "ambiguous_bases": 0}
    with sqlite3.connect(database) as connection:
        connection.execute("PRAGMA cache_size=-4096")
        connection.execute("PRAGMA temp_store=FILE")
        connection.execute("CREATE TABLE names (name TEXT PRIMARY KEY) WITHOUT ROWID")
        for kind, value in sequence_events(path, reference=reference,
                                          max_length=None if reference else config["max_read_length"]):
            if kind == "start":
                try:
                    connection.execute("INSERT INTO names VALUES (?)", (value,))
                except sqlite3.IntegrityError as exc:
                    raise ValidationError(f"Duplicate sequence identifier: {value}") from exc
                counts["sequence_count"] += 1
            elif kind == "bases":
                ambiguous_count = len(value.translate(None, b"ACGT"))
                counts["ambiguous_bases"] += ambiguous_count
                if ambiguous_count and config["ambiguous"] == "error":
                    raise ValidationError("Input has ambiguous bases forbidden by --ambiguous error")
            elif kind == "end":
                counts["base_count"] += value
                counts["max_sequence_length"] = max(counts["max_sequence_length"], value)
                if value < config["window"]:
                    counts["skipped_short_sequences"] += 1
                    if not reference and config["short_reads"] == "error":
                        raise ValidationError("Read is shorter than --window; no padding is defined")
                intervals = (full_window_count(value, config["window"], config["stride"])
                             if reference else len(uniform_seed_starts(value, config["window"], config["seeds"])))
                counts["seed_intervals"] += intervals
                counts["vector_count"] += intervals * (1 if reference else 2)
    return counts


class TorchEncoder:
    """A single loaded module; each call is a whole batch and one host transfer."""

    def __init__(self, artifact: Path, device: str, threads: int, seed: int):
        if not re.fullmatch(r"cpu|cuda(?::[0-9]+)?", device):
            raise ValueError("device must be cpu, cuda, or cuda:N")
        positive_int(threads, "threads")
        os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
        try:
            self.torch = importlib.import_module("torch")
        except ImportError as exc:
            raise RuntimeError("Torch is required for encoding; install a suitable trusted PyTorch build or the production[torch] extra") from exc
        torch = self.torch
        if device.startswith("cuda"):
            if not torch.cuda.is_available():
                raise RuntimeError("CUDA was requested but PyTorch reports no available CUDA device")
            index = int(device.split(":", 1)[1]) if ":" in device else 0
            if index >= torch.cuda.device_count():
                raise RuntimeError(f"CUDA device index {index} is unavailable")
            device = f"cuda:{index}"
        self.device = torch.device(device)
        torch.set_num_threads(threads)
        torch.manual_seed(seed)
        torch.use_deterministic_algorithms(True)
        torch.backends.cudnn.benchmark = False
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.allow_tf32 = False
        torch.backends.cuda.matmul.allow_tf32 = False
        self.model = torch.jit.load(str(artifact), map_location="cpu").eval().to(self.device)
        if self.device.type == "cuda":
            torch.cuda.synchronize(self.device)
        self.runtime = {
            "torch_version": str(torch.__version__), "device": str(self.device), "threads": threads,
            "seed": seed, "deterministic_algorithms": True, "allow_tf32": False,
            "cuda_version": torch.version.cuda,
            "device_name": torch.cuda.get_device_name(self.device) if self.device.type == "cuda" else "CPU",
            "cublas_workspace_config": os.environ.get("CUBLAS_WORKSPACE_CONFIG"),
        }

    def encode(self, sequences: list[bytes], window: int, ambiguous: str):
        torch = self.torch
        tensor = torch.from_numpy(one_hot_batch(sequences, window, ambiguous)).to(self.device)
        with torch.inference_mode():
            output = self.model(tensor)
        if not isinstance(output, torch.Tensor) or not output.is_floating_point():
            raise ValueError("Encoder must return one floating Tensor, not a tuple, dict, or integer tensor")
        if output.ndim != 2 or output.shape[0] != len(sequences) or output.shape[1] < 1:
            raise ValueError("Encoder output must have shape [batch, positive_dimension]")
        # Blocking transfer once per batch includes completion of CUDA work.
        values = output.detach().to(device="cpu", dtype=torch.float32).contiguous().numpy()
        return validate_vectors(values, len(sequences))


def encode_file(args: argparse.Namespace) -> dict:
    """Snapshot, validate/count, encode by batch, and atomically publish artifacts."""
    started = time.perf_counter()
    config = {
        "mode": args.mode, "window": positive_int(args.window, "window"),
        "stride": positive_int(args.stride, "stride") if args.mode == "reference" else None,
        "seeds": positive_int(args.seeds, "seeds") if args.mode == "reads" else None,
        "max_read_length": positive_int(args.max_read_length, "max_read_length") if args.mode == "reads" else None,
        "short_reads": args.short_reads if args.mode == "reads" else None,
        "orientations": ["+"] if args.mode == "reference" else ["+", "-"],
        "batch_size": positive_int(args.batch_size, "batch_size"), "ambiguous": args.ambiguous,
        "channel_order": "ACGT", "input_dtype": "float32", "output_dtype": "float32",
        "window_policy": "full_windows_only", "query_coordinates": "original_forward_half_open",
    }
    output = Path(args.output).expanduser().absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Output already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=".encode-", dir=output.parent))
    published = False
    try:
        with tempfile.TemporaryDirectory(prefix=".encoder-inputs-", dir=output.parent) as work:
            working = Path(work)
            source = snapshot_file(Path(args.input), working / "input.snapshot")
            encoder_source = snapshot_file(Path(args.encoder), working / "encoder.snapshot")
            counts = census(working / "input.snapshot", config, working / "names.sqlite")
            if not counts["vector_count"]:
                raise ValueError("No full encoder windows: every sequence is shorter than --window")
            load_started = time.perf_counter()
            runtime = TorchEncoder(working / "encoder.snapshot", args.device, args.threads, args.seed)
            model_load_seconds = time.perf_counter() - load_started
            import numpy as np

            events = sequence_events(working / "input.snapshot", reference=args.mode == "reference",
                                     max_length=config["max_read_length"])
            rows = (reference_windows(events, config["window"], config["stride"])
                    if args.mode == "reference" else read_windows(events, window=config["window"],
                        seeds=config["seeds"], max_read_length=config["max_read_length"],
                        short_reads=config["short_reads"]))
            metadata_name = "windows.jsonl" if args.mode == "reference" else "queries.jsonl"
            vectors, offset, dimension, encoding_seconds, batch_count = None, 0, None, 0.0, 0
            try:
                with (temporary / metadata_name).open("xb") as metadata:
                    for batch in batches(rows, config["batch_size"]):
                        encoding_started = time.perf_counter()
                        values = runtime.encode([sequence for sequence, _ in batch], config["window"], config["ambiguous"])
                        encoding_seconds += time.perf_counter() - encoding_started
                        if vectors is None:
                            dimension = int(values.shape[1])
                            vectors = np.lib.format.open_memmap(temporary / "vectors.npy", mode="w+",
                                dtype="<f4", shape=(counts["vector_count"], dimension))
                        if values.shape != (len(batch), dimension):
                            raise ValueError("Encoder dimension changed between batches")
                        if offset + len(batch) > counts["vector_count"]:
                            raise ValueError("More generated windows than the input census")
                        vectors[offset : offset + len(batch)] = values
                        for _sequence, row in batch:
                            metadata.write(canonical_bytes(row) + b"\n")
                        offset += len(batch)
                        batch_count += 1
                    if offset != counts["vector_count"]:
                        raise ValueError("Generated window count differs from the input census")
                    metadata.flush()
                    os.fsync(metadata.fileno())
            finally:
                if vectors is not None:
                    vectors.flush()
                    # Explicitly release the mmap before directory rename/cleanup.
                    vectors._mmap.close()
            with (temporary / "vectors.npy").open("rb") as stored:
                os.fsync(stored.fileno())
            files = {name: hash_file(temporary / name) for name in ("vectors.npy", metadata_name)}
            input_kind = "reference" if args.mode == "reference" else "reads"
            recipe = {"schema": SCHEMA, f"{input_kind}_sha256": source["sha256"],
                      "encoder_sha256": encoder_source["sha256"], "config": config,
                      "runtime": runtime.runtime, "numpy_version": np.__version__}
            manifest = {
                **recipe, "recipe_id": content_id(recipe), "source": source,
                "encoder_source": encoder_source, "counts": {**counts, "batch_count": batch_count},
                "dimension": dimension, "metadata_file": metadata_name, "vectors_file": "vectors.npy",
                "files": files, "output_id": content_id({"recipe": recipe, "files": files}),
                "timing": {"encoding_seconds": encoding_seconds, "model_load_seconds": model_load_seconds,
                           "full_elapsed_seconds": None},
            }
        manifest["timing"]["full_elapsed_seconds"] = time.perf_counter() - started
        write_json(temporary / "manifest.json", manifest)
        fsync_directory(temporary)
        # A completed output from a concurrent run is never overwritten.
        if output.exists() or output.is_symlink():
            raise FileExistsError(f"Output appeared while encoding: {output}")
        os.rename(temporary, output)
        published = True
        fsync_directory(output.parent)
        return manifest
    finally:
        if not published:
            shutil.rmtree(temporary, ignore_errors=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="mode", required=True)
    for mode in ("reference", "reads"):
        sub = subparsers.add_parser(mode)
        sub.add_argument("--input", type=Path, required=True, help="FASTA reference or FASTA/FASTQ reads, optionally gzip")
        sub.add_argument("--encoder", type=Path, required=True, help="Trusted trained TorchScript .pt artifact")
        sub.add_argument("--output", type=Path, required=True, help="New output directory; existing paths are rejected")
        sub.add_argument("--device", default="cpu", help="cpu, cuda, or cuda:N; requested CUDA never silently falls back")
        sub.add_argument("--window", type=int, default=512)
        sub.add_argument("--batch-size", type=int, default=256)
        sub.add_argument("--threads", type=int, default=1)
        sub.add_argument("--seed", type=int, default=17)
        sub.add_argument("--ambiguous", choices=("zero", "error"), default="zero")
        if mode == "reference":
            sub.add_argument("--stride", type=int, default=256)
        else:
            sub.add_argument("--seeds", type=int, default=8)
            sub.add_argument("--max-read-length", type=int, default=2_000_000)
            sub.add_argument("--short-reads", choices=("skip", "error"), default="skip")
    args = parser.parse_args(argv)
    try:
        manifest = encode_file(args)
    except (EngineError, ValueError, OSError, RuntimeError, ImportError, sqlite3.Error) as exc:
        parser.exit(2, f"encoder error: {exc}\n")
    print(canonical_bytes(manifest).decode("ascii"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
