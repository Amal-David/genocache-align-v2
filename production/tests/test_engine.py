from __future__ import annotations

from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import gzip
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys

import pysam
import pytest

from genocache.engine import Engine, ExecutionError, IntegrityError, SUPPORTED_PRESETS, ValidationError
from genocache.io import CHUNK_BYTES, snapshot_file, validate_sequences


@pytest.fixture
def binary() -> Path:
    requested = os.environ.get("MINIMAP2", "minimap2")
    found = shutil.which(requested)
    if found is None:
        local = Path(__file__).resolve().parents[3] / "minimap2-src" / "minimap2"
        found = str(local) if local.is_file() else None
    if found is None:
        pytest.skip("Real minimap2 is required for native integration tests")
    return Path(found).resolve()


@pytest.fixture
def samples(tmp_path: Path) -> tuple[Path, Path]:
    rng = random.Random(814031)
    def dna(count):
        return "".join(rng.choices("ACGT", k=count))
    reference = dna(36000)
    other_chromosome = dna(9000)
    repeated = dna(4200)
    def reverse(value):
        return value.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    records = [
        ("forward", reference[1200:4500]),
        ("reverse", reverse(reference[7600:11000])),
        ("chimeric", reference[15000:17200] + other_chromosome[1500:3800]),
        ("unmapped", dna(3100)),
        ("repetitive", repeated[200:4000]),
        ("indel", reference[19000:20900] + dna(77) + reference[20900:22300]),
    ]
    reference_path = tmp_path / "reference.fa"
    reference_path.write_text(f">chr1\n{reference}\n>chr2\n{other_chromosome}\n>repeat_a\n{repeated}\n>repeat_b\n{repeated}\n")
    reads_path = tmp_path / "reads.fq"
    reads_path.write_text("".join(
        f"@{name} original description\n{sequence}\n+\n"
        + "".join(chr(40 + offset % 30) for offset in range(len(sequence))) + "\n"
        for name, sequence in records
    ))
    return reference_path, reads_path


@pytest.fixture
def engine(tmp_path: Path, binary: Path) -> Engine:
    return Engine(tmp_path / "cache", minimap2=binary)


def _records(path: Path) -> list[tuple]:
    with pysam.AlignmentFile(str(path), "rb" if path.suffix == ".bam" else "r") as stream:
        return [
            (record.query_name, record.flag, record.reference_name, record.reference_start,
             record.mapping_quality, record.cigarstring, record.next_reference_name,
             record.next_reference_start, record.template_length, record.query_sequence,
             tuple(record.query_qualities or ()), tuple(record.get_tags(with_value_type=True)))
            for record in stream.fetch(until_eof=True)
        ]


def _native(binary: Path, pack: dict, reads: Path, destination: Path) -> None:
    with destination.open("wb") as output:
        subprocess.run([str(binary), "-a", "-x", pack["preset"], "-t", "1", pack["index_path"], str(reads)],
                       stdout=output, stderr=subprocess.PIPE, check=True)


def test_streamed_multiline_fastq_and_gzip(tmp_path: Path) -> None:
    # A read larger than the parser chunk forces bounded continuation handling.
    sequence = b"ACGT" * (CHUNK_BYTES // 4 + 17)
    content = b"@long\n" + sequence + b"\n+long\n" + b"I" * len(sequence) + b"\n@second\nAC\nGT\n+\n@@\nII\n"
    source = tmp_path / "reads.any-extension"
    source.write_bytes(gzip.compress(content))
    result = validate_sequences(source)
    assert result == {"format": "fastq", "gzip": True, "sequence_count": 2, "base_count": len(sequence) + 4}


@pytest.mark.parametrize("content", [
    b"@read\nACGT\n+\nIII\n", b"@read\nACGT\n+\nIIIII\n",
    b"@read\nACGT\n+other\nIIII\n", b">read\nACGX\n", b">empty\n",
    b"@read\nACGT\n+\nII I\n", b"@bad@name\nACGT\n+\nIIII\n",
])
def test_invalid_reads_fail_before_execution(tmp_path: Path, content: bytes) -> None:
    source = tmp_path / "bad.input"
    source.write_bytes(content)
    with pytest.raises(ValidationError):
        validate_sequences(source)


def test_duplicate_reference_identifiers_are_rejected(tmp_path: Path) -> None:
    source = tmp_path / "duplicate.fa"
    source.write_text(">same first\nACGT\n>same second\nTGCA\n")
    with pytest.raises(ValidationError, match="Duplicate reference"):
        validate_sequences(source, reference=True)


def test_truncated_gzip_is_rejected(tmp_path: Path) -> None:
    source = tmp_path / "damaged.gz"
    source.write_bytes(gzip.compress(b">read\nACGT\n")[:-5])
    with pytest.raises(ValidationError, match="gzip"):
        validate_sequences(source)


def test_snapshot_is_independent_of_subsequent_source_edits(tmp_path: Path) -> None:
    source = tmp_path / "original.fa"
    snapshot = tmp_path / "snapshot.fa"
    source.write_text(">r\nACGT\n")
    identity = snapshot_file(source, snapshot)
    source.write_text(">r\nTGCA\n")
    assert snapshot.read_text() == ">r\nACGT\n"
    assert identity["byte_count"] == 8
    assert snapshot.stat().st_mode & 0o222 == 0


def test_missing_binary_and_unsupported_preset_fail_closed(tmp_path: Path, engine: Engine, samples) -> None:
    with pytest.raises(ValidationError, match="missing"):
        Engine(tmp_path / "missing-cache", minimap2="/does/not/exist/minimap2")
    with pytest.raises(ValidationError, match="Unsupported preset"):
        engine.build_index(samples[0], preset="splice")
    assert not list((engine.root / "packs").iterdir())


def test_reference_cache_invalidation_and_relocation(engine: Engine, samples, tmp_path: Path, binary: Path) -> None:
    reference, _ = samples
    first = engine.build_index(reference)
    second = engine.build_index(reference)
    assert not first["cache_hit"] and second["cache_hit"]
    assert first["reference_id"] == second["reference_id"]
    assert second["metrics"]["native_index_seconds"] == 0
    manifest = json.loads(Path(first["manifest_path"]).read_text())
    assert manifest["identity"]["index_options"]["index_batch_bases"] == first["input_stats"]["base_count"] + 1
    assert all(not Path(item["path"]).is_absolute() for item in manifest["files"].values())
    assert all(path.stat().st_mode & 0o222 == 0 for path in Path(first["path"]).iterdir())

    moved = Engine(tmp_path / "moved", minimap2=binary)
    shutil.copytree(first["path"], moved.root / "packs" / first["reference_id"])
    assert moved.inspect_reference(first["reference_id"])["index_path"].startswith(str(moved.root))

    reference.write_text(reference.read_text().replace("AC", "AG", 1))
    changed = engine.build_index(reference)
    assert changed["reference_id"] != first["reference_id"]
    other_preset = engine.build_index(reference, preset="map-hifi")
    assert other_preset["reference_id"] != changed["reference_id"]


@pytest.mark.parametrize("preset", sorted(SUPPORTED_PRESETS))
def test_native_sam_is_preserved_exactly_at_record_level(engine: Engine, samples, binary: Path, tmp_path: Path, preset: str) -> None:
    reference, reads = samples
    pack = engine.build_index(reference, preset=preset)
    aligned = engine.align(pack["reference_id"], reads, output_format="sam")
    baseline = tmp_path / "native.sam"
    _native(binary, pack, reads, baseline)
    records = _records(Path(aligned["output_path"]))
    assert records == _records(baseline)
    assert any(row[1] & 16 for row in records)  # Reverse strand.
    assert any(row[1] & 256 for row in records)  # Secondary mapping.
    assert any(row[1] & 2048 for row in records)  # Chimeric supplementary mapping.
    assert any(row[1] & 4 for row in records)  # Unmapped read retained.
    assert aligned["output_stats"]["primary_records"] == 6
    assert aligned["metrics"]["native_alignment_seconds"] > 0
    assert aligned["metrics"]["request_wall_seconds"] >= aligned["metrics"]["native_alignment_seconds"]
    assert not any(path.name.startswith("reads.") for path in Path(aligned["job_path"]).iterdir())


def test_bam_is_sorted_indexed_and_semantically_matches_native(engine: Engine, samples, binary: Path, tmp_path: Path) -> None:
    reference, reads = samples
    pack = engine.build_index(reference)
    aligned = engine.align(pack["reference_id"], reads, output_format="bam")
    baseline = tmp_path / "native.sam"
    _native(binary, pack, reads, baseline)
    assert Counter(_records(Path(aligned["output_path"]))) == Counter(_records(baseline))
    assert Path(aligned["index_path"]).name.endswith(".csi")
    with pysam.AlignmentFile(aligned["output_path"], "rb") as stream:
        assert stream.header.to_dict()["HD"]["SO"] == "coordinate"
        assert list(stream.fetch("chr1", 1000, 5000))
    warmed = engine.align(pack["reference_id"], reads, output_format="bam")
    assert warmed["cache_hit"] and warmed["job_id"] == aligned["job_id"]
    assert warmed["metrics"]["native_alignment_seconds"] == 0
    assert warmed["metrics"]["conversion_seconds"] == 0
    assert warmed["metrics"]["request_wall_seconds"] > 0


def test_exact_job_key_includes_names_qualities_format_and_execution(engine: Engine, samples, tmp_path: Path) -> None:
    reference, reads = samples
    pack = engine.build_index(reference)
    original = engine.align(pack["reference_id"], reads, output_format="sam")
    same_bytes = tmp_path / "other-name.fq"
    shutil.copyfile(reads, same_bytes)
    assert engine.align(pack["reference_id"], same_bytes, output_format="sam")["job_id"] == original["job_id"]

    renamed = tmp_path / "renamed.fq"
    renamed.write_text(reads.read_text().replace("@forward ", "@renamed ", 1))
    name_job = engine.align(pack["reference_id"], renamed, output_format="sam")
    assert not name_job["cache_hit"] and name_job["job_id"] != original["job_id"]
    assert _records(Path(name_job["output_path"]))[0][0] == "renamed"

    quality_changed = tmp_path / "qualities.fq"
    lines = reads.read_text().splitlines()
    lines[3] = "J" * len(lines[3])
    quality_changed.write_text("\n".join(lines) + "\n")
    quality_job = engine.align(pack["reference_id"], quality_changed, output_format="sam")
    assert quality_job["job_id"] != original["job_id"]
    assert _records(Path(quality_job["output_path"]))[0][10] != _records(Path(original["output_path"]))[0][10]

    compressed = tmp_path / "reads.gz"
    compressed.write_bytes(gzip.compress(reads.read_bytes()))
    gzip_job = engine.align(pack["reference_id"], compressed, output_format="sam")
    assert gzip_job["job_id"] != original["job_id"]
    assert _records(Path(gzip_job["output_path"])) == _records(Path(original["output_path"]))
    assert engine.align(pack["reference_id"], reads, threads=2, output_format="sam")["job_id"] != original["job_id"]
    assert engine.align(pack["reference_id"], reads, output_format="bam")["job_id"] != original["job_id"]


def test_corrupt_reference_pack_is_rejected(engine: Engine, samples) -> None:
    pack = engine.build_index(samples[0])
    index = Path(pack["index_path"])
    index.chmod(0o644)
    with index.open("ab") as stream:
        stream.write(b"corruption")
    with pytest.raises(IntegrityError, match="checksum"):
        engine.align(pack["reference_id"], samples[1], output_format="sam")
    with pytest.raises(IntegrityError, match="checksum"):
        engine.build_index(samples[0])
    assert not list((engine.root / "jobs").iterdir())


def test_corrupt_cached_alignment_is_not_returned(engine: Engine, samples) -> None:
    pack = engine.build_index(samples[0])
    first = engine.align(pack["reference_id"], samples[1], output_format="sam")
    output = Path(first["output_path"])
    output.chmod(0o644)
    with output.open("ab") as stream:
        stream.write(b"corruption")
    with pytest.raises(IntegrityError, match="checksum"):
        engine.align(pack["reference_id"], samples[1], output_format="sam")
    assert not list((engine.root / ".tmp").iterdir())


def test_different_binary_hash_cannot_reuse_pack(engine: Engine, samples, binary: Path, tmp_path: Path) -> None:
    pack = engine.build_index(samples[0])
    changed_binary = tmp_path / "different-minimap2"
    shutil.copyfile(binary, changed_binary)
    # ELF executables permit trailing bytes; version output is unchanged while
    # exact executable provenance changes and must invalidate the old pack.
    with changed_binary.open("ab") as stream:
        stream.write(b"\nprovenance-test\n")
    changed_binary.chmod(0o755)
    other_engine = Engine(engine.root, minimap2=changed_binary)
    with pytest.raises(IntegrityError, match="different minimap2"):
        other_engine.align(pack["reference_id"], samples[1], output_format="sam")


def test_simultaneous_identical_jobs_publish_once(engine: Engine, samples) -> None:
    pack = engine.build_index(samples[0])
    with ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(lambda _: engine.align(pack["reference_id"], samples[1], output_format="sam"), range(2)))
    assert results[0]["job_id"] == results[1]["job_id"]
    assert sorted(result["cache_hit"] for result in results) == [False, True]
    assert len(list((engine.root / "jobs").iterdir())) == 1
    assert not list((engine.root / ".tmp").iterdir())


def test_failed_execution_never_publishes_partial_job(engine: Engine, samples, monkeypatch) -> None:
    pack = engine.build_index(samples[0])

    def fail(command, *, cwd, stderr_path, stdout_path=None, deadline=None):
        stdout_path.write_bytes(b"@HD\tVN:1.6\n")
        stderr_path.write_text("simulated native failure")
        raise ExecutionError("simulated failure")

    monkeypatch.setattr(engine, "_execute", fail)
    with pytest.raises(ExecutionError):
        engine.align(pack["reference_id"], samples[1], output_format="sam")
    assert not list((engine.root / "jobs").iterdir())
    assert not list((engine.root / ".tmp").iterdir())


def test_deadline_cleans_staging(engine: Engine, samples) -> None:
    pack = engine.build_index(samples[0])
    with pytest.raises(ExecutionError, match="deadline"):
        engine.align(pack["reference_id"], samples[1], timeout=0.000001)
    assert not list((engine.root / "jobs").iterdir())
    assert not list((engine.root / ".tmp").iterdir())


def test_running_subprocess_timeout_discards_partial_sam(binary: Path, samples, tmp_path: Path) -> None:
    wrapper = tmp_path / "slow-minimap2"
    wrapper.write_text(
        f"#!{sys.executable}\n"
        "import os, sys, time\n"
        "if '-a' in sys.argv:\n"
        "    print('@HD\\tVN:1.6', flush=True)\n"
        "    time.sleep(60)\n"
        f"os.execv({str(binary)!r}, [{str(binary)!r}, *sys.argv[1:]])\n"
    )
    wrapper.chmod(0o755)
    slow_engine = Engine(tmp_path / "slow-cache", minimap2=wrapper)
    pack = slow_engine.build_index(samples[0])
    with pytest.raises(ExecutionError, match="partial outputs were discarded"):
        slow_engine.align(pack["reference_id"], samples[1], timeout=0.2, output_format="sam")
    assert not list((slow_engine.root / "jobs").iterdir())
    assert not list((slow_engine.root / ".tmp").iterdir())


def test_cli_json_and_version(binary: Path, tmp_path: Path) -> None:
    version = subprocess.run([sys.executable, "-m", "genocache.cli", "--version"], capture_output=True, text=True, check=True)
    assert version.stdout.startswith("genocache ")
    doctor = subprocess.run([sys.executable, "-m", "genocache.cli", "--root", str(tmp_path / "cli"),
                             "--minimap2", str(binary), "doctor"], capture_output=True, text=True, check=True)
    result = json.loads(doctor.stdout)
    assert result["ok"] and result["minimap2"]["sha256"]
