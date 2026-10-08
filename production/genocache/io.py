"""Bounded-memory input validation and durable local-file primitives.

No sequence normalization happens here. The exact bytes that are hashed are
snapshotted and passed to minimap2, including headers, qualities, and gzip data.
"""

from __future__ import annotations

import contextlib
import gzip
import hashlib
import json
import os
from pathlib import Path
import stat
import time
from typing import Any, BinaryIO, Iterator
import zlib


CHUNK_BYTES = 1024 * 1024
DNA_ALPHABET = b"ACGTNRYSWKMBDHVacgtnryswkmbdhv"
QUALITY_ALPHABET = bytes(range(33, 127))


class EngineError(RuntimeError):
    """An expected, safely reportable production-engine error."""


class ValidationError(EngineError):
    """Invalid input or unsupported configuration."""


class IntegrityError(EngineError):
    """A snapshot, reference pack, or cached output failed verification."""


class ExecutionError(EngineError):
    """An external execution failed or exceeded its deadline."""


def check_deadline(deadline: float | None) -> None:
    if deadline is not None and time.perf_counter() >= deadline:
        raise ExecutionError("Execution deadline exceeded")


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False
    ).encode("utf-8")


def content_id(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def file_signature(path: Path) -> tuple[int, int, int, int, int]:
    """Detect ordinary changes/replacements while operating on an immutable file."""
    if path.is_symlink():
        raise IntegrityError(f"Symlinks are not allowed in immutable artifacts: {path.name}")
    value = path.stat()
    if not stat.S_ISREG(value.st_mode):
        raise IntegrityError(f"Expected a regular file: {path.name}")
    return _stat_signature(value)


def _stat_signature(value: os.stat_result) -> tuple[int, int, int, int, int]:
    return value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns


def hash_file(path: Path, deadline: float | None = None) -> dict[str, Any]:
    before = file_signature(path)
    digest = hashlib.sha256()
    byte_count = 0
    with path.open("rb") as stream:
        if _stat_signature(os.fstat(stream.fileno())) != before:
            raise IntegrityError(f"File changed before hashing: {path.name}")
        while chunk := stream.read(CHUNK_BYTES):
            check_deadline(deadline)
            digest.update(chunk)
            byte_count += len(chunk)
        after_fd = _stat_signature(os.fstat(stream.fileno()))
    if before != after_fd or before != file_signature(path):
        raise IntegrityError(f"File changed while hashing: {path.name}")
    return {"sha256": digest.hexdigest(), "byte_count": byte_count}


def snapshot_file(source: Path, target: Path, deadline: float | None = None) -> dict[str, Any]:
    """Take one verified snapshot; never execute against the caller's mutable path."""
    try:
        source = source.expanduser().resolve(strict=True)
        before = file_signature(source)
    except (OSError, RuntimeError) as exc:
        if isinstance(exc, EngineError):
            raise
        raise ValidationError("Input must be an existing readable regular file") from exc
    digest = hashlib.sha256()
    byte_count = 0
    try:
        with source.open("rb") as stream, target.open("xb") as output:
            if _stat_signature(os.fstat(stream.fileno())) != before:
                raise IntegrityError("Input changed before snapshotting; retry with an immutable file")
            while chunk := stream.read(CHUNK_BYTES):
                check_deadline(deadline)
                output.write(chunk)
                digest.update(chunk)
                byte_count += len(chunk)
            after_fd = _stat_signature(os.fstat(stream.fileno()))
            output.flush()
            os.fsync(output.fileno())
        if before != after_fd or before != file_signature(source) or byte_count != before[2]:
            raise IntegrityError("Input changed during snapshotting; retry with an immutable file")
        target.chmod(0o444)
        return {"sha256": digest.hexdigest(), "byte_count": byte_count}
    except OSError as exc:
        raise ValidationError("Unable to snapshot the input file") from exc


@contextlib.contextmanager
def _sequence_stream(path: Path) -> Iterator[tuple[BinaryIO, bool]]:
    try:
        with path.open("rb") as raw:
            compressed = raw.read(2) == b"\x1f\x8b"
            raw.seek(0)
            if compressed:
                with gzip.GzipFile(fileobj=raw, mode="rb") as stream:
                    yield stream, True
            else:
                yield raw, False
    except (gzip.BadGzipFile, EOFError, zlib.error, OSError) as exc:
        raise ValidationError("Unreadable input or damaged gzip stream") from exc


def _line_pieces(stream: BinaryIO, deadline: float | None) -> Iterator[tuple[bytes, bool, bool, int]]:
    """A line iterator that does not allocate an entire unwrapped chromosome."""
    at_start, line_number = True, 1
    while raw := stream.readline(CHUNK_BYTES):
        check_deadline(deadline)
        # A CRLF pair may straddle the bounded read. Consume at most one extra
        # byte so CR is only accepted as an actual line terminator.
        reached_eof_after_cr = False
        if len(raw) == CHUNK_BYTES and raw.endswith(b"\r"):
            following = stream.read(1)
            reached_eof_after_cr = not following
            raw += following
        at_end = raw.endswith(b"\n") or len(raw) < CHUNK_BYTES or reached_eof_after_cr
        value = raw[:-1] if raw.endswith(b"\n") else raw
        if at_end and value.endswith(b"\r"):
            value = value[:-1]
        yield value, at_start, at_end, line_number
        if at_end:
            line_number += 1
        at_start = at_end


def _header_name(value: bytes, at_end: bool, line_number: int, reference: bool) -> str:
    if not at_end or len(value) > CHUNK_BYTES:
        raise ValidationError(f"Header is too long at line {line_number}")
    body = value[1:]
    if not body or body[:1].isspace() or any((c < 32 and c != 9) or c > 126 for c in body):
        raise ValidationError(f"Expected a nonempty ASCII sequence identifier at line {line_number}")
    name = body.split(None, 1)[0]
    if (not reference and (len(name) > 254 or b"@" in name)) or (
        reference and name[:1] in (b"*", b"=")
    ):
        raise ValidationError(f"Sequence identifier is not SAM compatible at line {line_number}")
    return name.decode("ascii")


def validate_sequences(
    path: Path, *, reference: bool = False, deadline: float | None = None
) -> dict[str, Any]:
    """Validate strict DNA FASTA/FASTQ, streaming both bases and quality strings.

    References must be FASTA with unique first-token names. Reads can be FASTA
    or multiline FASTQ. Empty sequences, non-IUPAC DNA bases, malformed quality
    lengths, and non-SAM-compatible names are rejected before native execution.
    Memory is bounded by a chunk plus reference-name metadata, not sequence size.
    """
    sequence_count = base_count = current_length = quality_length = 0
    current_name: str | None = None
    sequence_format: str | None = None
    state = "header"
    contigs: list[dict[str, Any]] = []
    seen_names: set[str] = set()

    def finish_sequence() -> None:
        if current_length == 0:
            raise ValidationError("Empty sequences are not supported")
        if current_length >= 2**31:
            raise ValidationError("Individual sequences must be shorter than 2^31 bases")
        if reference:
            contigs.append({"name": current_name, "length": current_length})

    with _sequence_stream(path) as (stream, compressed):
        for value, at_start, at_end, line_number in _line_pieces(stream, deadline):
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
                        finish_sequence()
                    current_name = _header_name(value, at_end, line_number, reference)
                    if reference and current_name in seen_names:
                        raise ValidationError(f"Duplicate reference sequence identifier at line {line_number}")
                    if reference:
                        seen_names.add(current_name)
                    sequence_count += 1
                    current_length = 0
                else:
                    if value.translate(None, DNA_ALPHABET):
                        raise ValidationError(f"Non-IUPAC DNA character at line {line_number}")
                    current_length += len(value)
                    base_count += len(value)
                continue

            if state == "header":
                if not value:
                    continue
                if not at_start or not value.startswith(b"@"):
                    raise ValidationError(f"Expected a FASTQ header at line {line_number}")
                current_name = _header_name(value, at_end, line_number, False)
                current_length = quality_length = 0
                sequence_count += 1
                state = "sequence"
            elif state == "sequence":
                if at_start and value.startswith(b"+"):
                    if not at_end:
                        raise ValidationError(f"FASTQ separator is too long at line {line_number}")
                    if current_length == 0:
                        raise ValidationError("Empty sequences are not supported")
                    repeated_name = value[1:].split(None, 1)
                    if repeated_name and repeated_name[0] != current_name.encode("ascii"):
                        raise ValidationError(f"FASTQ separator identifier mismatch at line {line_number}")
                    state = "quality"
                else:
                    if value.translate(None, DNA_ALPHABET):
                        raise ValidationError(f"Non-IUPAC DNA character at line {line_number}")
                    current_length += len(value)
                    base_count += len(value)
            else:
                if value.translate(None, QUALITY_ALPHABET):
                    raise ValidationError(f"Invalid FASTQ quality character at line {line_number}")
                quality_length += len(value)
                if quality_length > current_length:
                    raise ValidationError(f"FASTQ sequence/quality length mismatch at line {line_number}")
                if quality_length == current_length and at_end:
                    finish_sequence()
                    state = "header"

    if not sequence_count:
        raise ValidationError("Input contains no sequences")
    if sequence_format == "fasta":
        finish_sequence()
    elif state == "quality" and quality_length == current_length:
        finish_sequence()
    elif state != "header":
        raise ValidationError("Truncated FASTQ record or sequence/quality length mismatch")
    result: dict[str, Any] = {
        "format": sequence_format,
        "gzip": compressed,
        "sequence_count": sequence_count,
        "base_count": base_count,
    }
    if reference:
        result["contigs"] = contigs
    return result


def write_json(path: Path, value: Any) -> None:
    with path.open("xb") as output:
        output.write(canonical_bytes(value) + b"\n")
        output.flush()
        os.fsync(output.fileno())


def read_json(path: Path) -> dict[str, Any]:
    try:
        file_signature(path)
        # Metadata should never approach input-data scale.
        if path.stat().st_size > 64 * 1024 * 1024:
            raise IntegrityError("Artifact manifest exceeds the metadata size limit")
        with path.open("rb") as stream:
            value = json.load(stream)
        if not isinstance(value, dict):
            raise IntegrityError("Artifact manifest is not an object")
        return value
    except (OSError, ValueError) as exc:
        raise IntegrityError("Missing or invalid artifact manifest") from exc


def fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _convert_main() -> None:
    """Internal subprocess boundary so sort/index operations have a deadline."""
    import sys
    import pysam

    if len(sys.argv) != 6 or sys.argv[1] != "convert":
        raise SystemExit("Internal usage: python -m genocache.io convert SAM BAM THREADS SORT_MEMORY")
    source, destination, threads, memory = sys.argv[2:]
    pysam.sort("-@", threads, "-m", memory, "-O", "BAM", "-o", destination, source)
    pysam.index("-c", "-@", threads, destination)
    pysam.quickcheck(destination)


if __name__ == "__main__":
    _convert_main()
