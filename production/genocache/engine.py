"""Conservative, exact-input batch alignment with immutable disk artifacts.

The reference index is cached on disk. Every cold job starts a fresh minimap2
process, which loads that index once for the batch. This is not an in-memory
minimap2 daemon and it never borrows an alignment from a similar read.
"""

from __future__ import annotations

import contextlib
from datetime import datetime, timezone
import fcntl
from functools import wraps
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from typing import Any, Iterator

from .io import (
    EngineError as EngineError,
    ExecutionError,
    IntegrityError,
    ValidationError,
    check_deadline,
    content_id,
    file_signature,
    fsync_directory,
    hash_file,
    read_json,
    snapshot_file,
    validate_sequences,
    write_json,
)


SCHEMA_VERSION = 1
SUPPORTED_PRESETS = frozenset({"map-ont", "lr:hq", "map-hifi", "map-pb"})
SORT_MEMORY = "256M"
ID_PATTERN = re.compile(r"^[0-9a-f]{64}$")
_NATIVE_SUMMARY = re.compile(
    r"Real time:\s*([\d.]+) sec; CPU:\s*([\d.]+) sec; Peak RSS:\s*([\d.]+) GB"
)


def _utcnow() -> str:
    return datetime.now(timezone.utc).isoformat()


def _measure_request(function):
    """Include lock release and staging cleanup in the returned wall time."""
    @wraps(function)
    def measured(*args, **kwargs):
        started = time.perf_counter()
        result = function(*args, **kwargs)
        result["metrics"]["request_wall_seconds"] = time.perf_counter() - started
        return result
    return measured


def _threads(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= 256:
        raise ValidationError("threads must be an integer between 1 and 256")
    return value


def _deadline(start: float, timeout: float | None) -> float | None:
    if timeout is None:
        return None
    if isinstance(timeout, bool) or not isinstance(timeout, (int, float)) or not math.isfinite(timeout) or timeout <= 0:
        raise ValidationError("timeout must be a finite positive number of seconds")
    return start + timeout


def _require_id(value: str) -> str:
    if not isinstance(value, str) or not ID_PATTERN.fullmatch(value):
        raise ValidationError("Artifact identifier must be a lowercase 64-character SHA-256 digest")
    return value


def _tail(path: Path, limit: int = 8192) -> str:
    with path.open("rb") as stream:
        stream.seek(max(0, path.stat().st_size - limit))
        return stream.read(limit).decode("utf-8", errors="replace")


def _native_report(path: Path) -> dict[str, float] | None:
    matches = list(_NATIVE_SUMMARY.finditer(_tail(path)))
    if not matches:
        return None
    real, cpu, rss = matches[-1].groups()
    return {"real_seconds": float(real), "cpu_seconds": float(cpu), "peak_rss_gb": float(rss)}


class Engine:
    """A single-machine, POSIX-local-filesystem engine.

    ``root`` must be on the worker's local filesystem. Do not point this class at
    a Modal Volume, NFS, an object-store mount, or a shared multi-host cache.
    Cloud workers may securely stage verified packs into their own local roots.
    """

    def __init__(self, root: Path, minimap2: str | Path = "minimap2") -> None:
        self.root = Path(root).expanduser().resolve()
        self.minimap2 = str(minimap2)
        self.root.mkdir(parents=True, exist_ok=True)
        for name in ("packs", "jobs", ".locks", ".tmp"):
            path = self.root / name
            if path.is_symlink():
                raise ValidationError("Engine directories must not be symlinks")
            path.mkdir(exist_ok=True)
        self._check_local_filesystem()
        self._binary_identity()

    def _check_local_filesystem(self) -> None:
        # Linux can reject common unsupported mounts. This is not a distributed
        # lock probe; the caller must still satisfy the documented local-root
        # requirement on platforms without mountinfo or unfamiliar filesystems.
        mountinfo = Path("/proc/self/mountinfo")
        if not mountinfo.exists():
            return
        selected: tuple[int, str] | None = None
        for line in mountinfo.read_text().splitlines():
            left, separator, right = line.partition(" - ")
            if not separator:
                continue
            fields = left.split()
            if len(fields) < 5 or not right.split():
                continue
            mount = Path(fields[4].replace("\\040", " ").replace("\\134", "\\"))
            if self.root == mount or mount in self.root.parents:
                candidate = (len(str(mount)), right.split()[0])
                if selected is None or candidate[0] > selected[0]:
                    selected = candidate
        if selected and (
            selected[1] in {"nfs", "nfs4", "cifs", "smb3", "smbfs", "ceph", "glusterfs", "lustre", "9p"}
            or selected[1].startswith("fuse")
        ):
            raise ValidationError("Engine root must be a local filesystem supporting POSIX locks and atomic rename")

    def _binary_identity(self) -> tuple[Path, dict[str, str]]:
        found = shutil.which(self.minimap2)
        if found is None:
            raise ValidationError("minimap2 executable is missing or not executable")
        binary = Path(found).resolve()
        before = file_signature(binary)
        fingerprint = hash_file(binary)
        try:
            result = subprocess.run(
                [str(binary), "--version"], check=True, capture_output=True, timeout=10,
                env={**os.environ, "LC_ALL": "C"},
            )
        except (OSError, subprocess.SubprocessError) as exc:
            raise ExecutionError("Unable to execute minimap2 --version") from exc
        version = result.stdout.decode("utf-8", errors="replace").strip()
        if not version or len(version) > 256 or "\n" in version:
            raise ValidationError("minimap2 returned an invalid version string")
        if before != file_signature(binary):
            raise IntegrityError("minimap2 binary changed during version verification")
        return binary, {"sha256": fingerprint["sha256"], "version": version}

    @contextlib.contextmanager
    def _lock(self, namespace: str, identifier: str, deadline: float | None) -> Iterator[float]:
        start = time.perf_counter()
        path = self.root / ".locks" / f"{namespace}-{identifier}.lock"
        descriptor = os.open(path, os.O_CREAT | os.O_RDWR | getattr(os, "O_NOFOLLOW", 0), 0o600)
        try:
            while True:
                check_deadline(deadline)
                try:
                    fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    break
                except BlockingIOError:
                    time.sleep(0.025)
            yield time.perf_counter() - start
        finally:
            fcntl.flock(descriptor, fcntl.LOCK_UN)
            os.close(descriptor)

    @contextlib.contextmanager
    def _staging(self, prefix: str) -> Iterator[Path]:
        stage = Path(tempfile.mkdtemp(prefix=prefix, dir=self.root / ".tmp"))
        try:
            yield stage
        finally:
            if stage.exists():
                shutil.rmtree(stage)

    def _execute(
        self, command: list[str], *, cwd: Path, stderr_path: Path,
        stdout_path: Path | None = None, deadline: float | None = None,
    ) -> float:
        check_deadline(deadline)
        start = time.perf_counter()
        with stderr_path.open("xb") as stderr, contextlib.ExitStack() as stack:
            stdout: Any = subprocess.DEVNULL
            if stdout_path is not None:
                stdout = stack.enter_context(stdout_path.open("xb"))
            try:
                process = subprocess.Popen(
                    command, cwd=cwd, stdout=stdout, stderr=stderr,
                    env={**os.environ, "LC_ALL": "C"}, start_new_session=True,
                )
            except OSError as exc:
                raise ExecutionError("Unable to start the alignment subprocess") from exc
            try:
                returncode = process.wait(timeout=None if deadline is None else max(0.001, deadline - time.perf_counter()))
            except BaseException as exc:
                # Cancel the entire process group, including any helper process.
                with contextlib.suppress(ProcessLookupError):
                    os.killpg(process.pid, signal.SIGKILL)
                process.wait()
                if isinstance(exc, subprocess.TimeoutExpired):
                    raise ExecutionError("Alignment execution deadline exceeded; partial outputs were discarded") from exc
                raise
            if returncode != 0:
                stderr.flush()
                detail = _tail(stderr_path, 2048).strip()
                raise ExecutionError(f"Alignment subprocess failed with exit code {returncode}: {detail}")
        check_deadline(deadline)
        return time.perf_counter() - start

    def _publish(self, stage: Path, target: Path) -> None:
        if target.exists() or target.is_symlink():
            raise IntegrityError("Refusing to overwrite an immutable artifact")
        for path in stage.iterdir():
            if path.is_symlink() or not path.is_file():
                raise IntegrityError("Unexpected file in artifact staging directory")
            with path.open("rb") as stream:
                os.fsync(stream.fileno())
            path.chmod(0o444)
        fsync_directory(stage)
        os.rename(stage, target)
        fsync_directory(target.parent)

    @staticmethod
    def _file_record(path: Path, deadline: float | None = None) -> dict[str, Any]:
        return {"path": path.name, **hash_file(path, deadline)}

    @staticmethod
    def _verify_files(
        directory: Path, files: dict[str, Any], expected_paths: dict[str, str], deadline: float | None,
    ) -> dict[str, tuple[int, int, int, int, int]]:
        if directory.is_symlink() or not directory.is_dir() or set(files) != set(expected_paths):
            raise IntegrityError("Incomplete or invalid immutable artifact")
        signatures = {}
        for label, name in expected_paths.items():
            entry = files[label]
            if not isinstance(entry, dict) or entry.get("path") != name:
                raise IntegrityError("Unsafe or unexpected artifact file path")
            path = directory / name
            try:
                actual = hash_file(path, deadline)
                if actual["sha256"] != entry["sha256"] or actual["byte_count"] != entry["byte_count"]:
                    raise IntegrityError(f"Artifact checksum mismatch: {name}")
                signatures[name] = file_signature(path)
            except (OSError, KeyError) as exc:
                raise IntegrityError(f"Missing or invalid artifact file: {name}") from exc
        return signatures

    @staticmethod
    def _assert_unchanged(directory: Path, signatures: dict[str, tuple[int, int, int, int, int]]) -> None:
        try:
            if any(file_signature(directory / name) != signature for name, signature in signatures.items()):
                raise IntegrityError("Immutable artifact changed during execution")
        except OSError as exc:
            raise IntegrityError("Immutable artifact disappeared during execution") from exc

    def _read_pack(
        self, reference_id: str, binary_identity: dict[str, str], deadline: float | None,
    ) -> tuple[dict[str, Any], dict[str, tuple[int, int, int, int, int]]]:
        directory = self.root / "packs" / _require_id(reference_id)
        manifest = read_json(directory / "manifest.json")
        try:
            identity = manifest["identity"]
            stats = manifest["input_stats"]
            if (
                manifest["schema_version"] != SCHEMA_VERSION
                or manifest["kind"] != "reference_pack"
                or manifest["reference_id"] != reference_id
                or content_id(identity) != reference_id
                or identity["schema_version"] != SCHEMA_VERSION
                or identity["preset"] not in SUPPORTED_PRESETS
                or identity["reference"]["format"] != "fasta"
                or identity["reference"]["sequence_count"] != stats["sequence_count"]
                or identity["reference"]["base_count"] != stats["base_count"]
                or identity["reference"]["contigs_sha256"] != content_id(stats["contigs"])
            ):
                raise IntegrityError("Reference pack identity or schema is invalid")
            if identity["minimap2"] != binary_identity:
                raise IntegrityError("Reference pack was built with a different minimap2 binary/version; rebuild it with this worker")
            name = "reference.fa.gz" if identity["reference"]["gzip"] else "reference.fa"
            files = manifest["files"]
            if (
                files["reference"]["sha256"] != identity["reference"]["sha256"]
                or files["reference"]["byte_count"] != identity["reference"]["byte_count"]
                or len(stats["contigs"]) != stats["sequence_count"]
                or sum(item["length"] for item in stats["contigs"]) != stats["base_count"]
                or identity["index_options"]["index_batch_bases"] != stats["base_count"] + 1
            ):
                raise IntegrityError("Reference pack metadata is inconsistent")
            signatures = self._verify_files(
                directory, files,
                {"reference": name, "index": "index.mmi", "native_stderr": "index.stderr.log"}, deadline,
            )
            with (directory / "index.mmi").open("rb") as stream:
                if stream.read(4) != b"MMI\x02":
                    raise IntegrityError("Reference index has an invalid minimap2 index signature")
            return manifest, signatures
        except (KeyError, TypeError, ValueError) as exc:
            raise IntegrityError("Reference pack manifest is malformed") from exc

    def _pack_result(self, manifest: dict[str, Any], *, cache_hit: bool, metrics: dict[str, Any]) -> dict[str, Any]:
        directory = self.root / "packs" / manifest["reference_id"]
        return {
            "reference_id": manifest["reference_id"], "cache_hit": cache_hit,
            "path": str(directory), "manifest_path": str(directory / "manifest.json"),
            "reference_path": str(directory / manifest["files"]["reference"]["path"]),
            "index_path": str(directory / "index.mmi"), "preset": manifest["identity"]["preset"],
            "input_stats": manifest["input_stats"], "metrics": metrics,
        }

    @_measure_request
    def inspect_reference(self, reference_id: str) -> dict[str, Any]:
        """Verify all pack checksums and compatibility before returning its paths."""
        start = time.perf_counter()
        _, binary_identity = self._binary_identity()
        manifest, _ = self._read_pack(reference_id, binary_identity, None)
        return self._pack_result(manifest, cache_hit=True, metrics={"request_wall_seconds": time.perf_counter() - start})

    @_measure_request
    def build_index(
        self, reference: Path, preset: str = "map-ont", threads: int = 1,
        timeout: float | None = None,
    ) -> dict[str, Any]:
        start = time.perf_counter()
        deadline = _deadline(start, timeout)
        threads = _threads(threads)
        if preset not in SUPPORTED_PRESETS:
            raise ValidationError(f"Unsupported preset; choose one of {', '.join(sorted(SUPPORTED_PRESETS))}")
        binary, binary_identity = self._binary_identity()
        binary_signature = file_signature(binary)
        metrics: dict[str, Any] = {}
        with self._staging("index-") as stage:
            snapshot = stage / "reference.input"
            phase = time.perf_counter()
            source = snapshot_file(Path(reference), snapshot, deadline)
            metrics["snapshot_seconds"] = time.perf_counter() - phase
            phase = time.perf_counter()
            stats = validate_sequences(snapshot, reference=True, deadline=deadline)
            metrics["validation_seconds"] = time.perf_counter() - phase
            name = "reference.fa.gz" if stats["gzip"] else "reference.fa"
            snapshot.rename(stage / name)
            identity = {
                "schema_version": SCHEMA_VERSION, "kind": "reference_pack",
                "reference": {**source, **{key: stats[key] for key in ("format", "gzip", "sequence_count", "base_count")},
                              "contigs_sha256": content_id(stats["contigs"])},
                "preset": preset, "minimap2": binary_identity,
                # One index part is required for valid single-pass SAM/MAPQ.
                "index_options": {"index_batch_bases": stats["base_count"] + 1, "threads": threads},
            }
            reference_id = content_id(identity)
            target = self.root / "packs" / reference_id
            with self._lock("pack", reference_id, deadline) as waited:
                metrics["lock_wait_seconds"] = waited
                if target.exists() or target.is_symlink():
                    phase = time.perf_counter()
                    manifest, _ = self._read_pack(reference_id, binary_identity, deadline)
                    metrics["integrity_verification_seconds"] = time.perf_counter() - phase
                    metrics["native_index_seconds"] = 0.0
                    metrics["original_build_metrics"] = manifest["build_metrics"]
                    metrics["request_wall_seconds"] = time.perf_counter() - start
                    return self._pack_result(manifest, cache_hit=True, metrics=metrics)
                command = [
                    str(binary), "-x", preset, "-t", str(threads), "-I", str(stats["base_count"] + 1),
                    "-d", "index.mmi", name,
                ]
                metrics["native_index_seconds"] = self._execute(
                    command, cwd=stage, stderr_path=stage / "index.stderr.log", deadline=deadline,
                )
                if binary_signature != file_signature(binary):
                    raise IntegrityError("minimap2 binary changed while building the index")
                index = stage / "index.mmi"
                if not index.is_file() or index.stat().st_size < 24:
                    raise ExecutionError("minimap2 did not create a complete index")
                with index.open("rb") as stream:
                    if stream.read(4) != b"MMI\x02":
                        raise IntegrityError("minimap2 produced an invalid index signature")
                metrics["native_reported"] = _native_report(stage / "index.stderr.log")
                phase = time.perf_counter()
                files = {
                    "reference": self._file_record(stage / name, deadline),
                    "index": self._file_record(index, deadline),
                    "native_stderr": self._file_record(stage / "index.stderr.log", deadline),
                }
                metrics["integrity_verification_seconds"] = time.perf_counter() - phase
                if files["reference"]["sha256"] != source["sha256"]:
                    raise IntegrityError("Reference snapshot changed during index construction")
                metrics["request_wall_seconds"] = time.perf_counter() - start
                manifest = {
                    "schema_version": SCHEMA_VERSION, "kind": "reference_pack", "reference_id": reference_id,
                    "identity": identity, "input_stats": stats, "files": files,
                    "created_at": _utcnow(),
                    "build_metrics": {**{key: value for key, value in metrics.items() if key != "request_wall_seconds"},
                                      "prepublish_wall_seconds": metrics["request_wall_seconds"]},
                    "command": ["minimap2", *command[1:]],
                }
                write_json(stage / "manifest.json", manifest)
                check_deadline(deadline)
                self._publish(stage, target)
                metrics["request_wall_seconds"] = time.perf_counter() - start
                return self._pack_result(manifest, cache_hit=False, metrics=metrics)

    @staticmethod
    def _dependency_identity() -> dict[str, str]:
        try:
            import pysam
        except ImportError as exc:
            raise ValidationError("pysam is required for alignment validation and BAM conversion") from exc
        return {"pysam": pysam.__version__, "samtools": pysam.__samtools_version__}

    @staticmethod
    def _validate_output(
        output: Path, output_format: str, reference_stats: dict[str, Any], reads_stats: dict[str, Any],
        deadline: float | None,
    ) -> dict[str, int]:
        import pysam

        totals = {"records": 0, "primary_records": 0, "mapped_primary_records": 0,
                  "unmapped_primary_records": 0, "secondary_records": 0, "supplementary_records": 0,
                  "primary_bases": 0}
        try:
            if output_format == "bam":
                pysam.quickcheck(str(output))
            with pysam.AlignmentFile(str(output), "rb" if output_format == "bam" else "r") as stream:
                expected = reference_stats["contigs"]
                if list(stream.references) != [item["name"] for item in expected] or list(stream.lengths) != [item["length"] for item in expected]:
                    raise IntegrityError("Alignment header does not match the reference pack")
                if output_format == "bam":
                    if stream.header.to_dict().get("HD", {}).get("SO") != "coordinate" or not stream.check_index():
                        raise IntegrityError("BAM must be coordinate sorted with a readable CSI index")
                for record in stream.fetch(until_eof=True):
                    check_deadline(deadline)
                    totals["records"] += 1
                    if record.is_secondary:
                        totals["secondary_records"] += 1
                    if record.is_supplementary:
                        totals["supplementary_records"] += 1
                    if not record.is_secondary and not record.is_supplementary:
                        totals["primary_records"] += 1
                        key = "unmapped_primary_records" if record.is_unmapped else "mapped_primary_records"
                        totals[key] += 1
                        totals["primary_bases"] += record.infer_read_length() or record.query_length or 0
            if totals["primary_records"] != reads_stats["sequence_count"] or totals["primary_bases"] != reads_stats["base_count"]:
                raise IntegrityError("Alignment is incomplete: primary-record count or total read length does not match the input")
        except (OSError, ValueError, pysam.SamtoolsError) as exc:
            raise IntegrityError("Native alignment output failed SAM/BAM validation") from exc
        return totals

    def _read_job(
        self, job_id: str, identity: dict[str, Any], deadline: float | None,
    ) -> dict[str, Any]:
        directory = self.root / "jobs" / job_id
        manifest = read_json(directory / "manifest.json")
        try:
            if (
                manifest["schema_version"] != SCHEMA_VERSION or manifest["kind"] != "alignment_job"
                or manifest["job_id"] != job_id or manifest["identity"] != identity
                or content_id(manifest["identity"]) != job_id
                or manifest["output_stats"]["primary_records"] != identity["reads"]["sequence_count"]
                or manifest["output_stats"]["primary_bases"] != identity["reads"]["base_count"]
            ):
                raise IntegrityError("Cached job identity or output metadata is invalid")
            output_format = identity["execution"]["output_format"]
            expected = {"alignment": f"alignment.{output_format}", "native_stderr": "alignment.stderr.log"}
            if output_format == "bam":
                expected.update({"index": "alignment.bam.csi", "conversion_stderr": "conversion.stderr.log"})
            self._verify_files(directory, manifest["files"], expected, deadline)
            return manifest
        except (KeyError, TypeError, ValueError) as exc:
            raise IntegrityError("Cached job manifest is malformed") from exc

    def _job_result(self, manifest: dict[str, Any], *, cache_hit: bool, metrics: dict[str, Any]) -> dict[str, Any]:
        directory = self.root / "jobs" / manifest["job_id"]
        return {
            "job_id": manifest["job_id"], "reference_id": manifest["identity"]["reference_id"],
            "cache_hit": cache_hit, "job_path": str(directory), "manifest_path": str(directory / "manifest.json"),
            "output_path": str(directory / manifest["files"]["alignment"]["path"]),
            "index_path": str(directory / manifest["files"]["index"]["path"]) if "index" in manifest["files"] else None,
            "output_format": manifest["identity"]["execution"]["output_format"],
            "input_stats": manifest["input_stats"], "output_stats": manifest["output_stats"], "metrics": metrics,
        }

    @_measure_request
    def align(
        self, reference_id: str, reads: Path, threads: int = 1, timeout: float | None = None,
        output_format: str = "bam",
    ) -> dict[str, Any]:
        start = time.perf_counter()
        deadline = _deadline(start, timeout)
        threads = _threads(threads)
        _require_id(reference_id)
        if output_format not in {"sam", "bam"}:
            raise ValidationError("output_format must be 'sam' or 'bam'")
        binary, binary_identity = self._binary_identity()
        binary_signature = file_signature(binary)
        dependency_identity = self._dependency_identity()
        metrics: dict[str, Any] = {}
        phase = time.perf_counter()
        pack, pack_signatures = self._read_pack(reference_id, binary_identity, deadline)
        metrics["reference_verification_seconds"] = time.perf_counter() - phase
        pack_directory = self.root / "packs" / reference_id
        with self._staging("align-") as stage:
            snapshot = stage / "reads.input"
            phase = time.perf_counter()
            source = snapshot_file(Path(reads), snapshot, deadline)
            metrics["snapshot_seconds"] = time.perf_counter() - phase
            phase = time.perf_counter()
            stats = validate_sequences(snapshot, deadline=deadline)
            metrics["validation_seconds"] = time.perf_counter() - phase
            name = ("reads.fa" if stats["format"] == "fasta" else "reads.fq") + (".gz" if stats["gzip"] else "")
            snapshot = snapshot.rename(stage / name)
            snapshot_signature = file_signature(snapshot)
            preset = pack["identity"]["preset"]
            identity = {
                "schema_version": SCHEMA_VERSION, "kind": "alignment_job", "reference_id": reference_id,
                "reference_index_sha256": pack["files"]["index"]["sha256"],
                "reads": {**source, **stats}, "minimap2": binary_identity,
                "execution": {"preset": preset, "threads": threads, "output_format": output_format,
                              "native_options": ["-a"], "timeout_seconds": timeout,
                              "sort_memory_per_thread": SORT_MEMORY if output_format == "bam" else None,
                              "bam_index_type": "CSI" if output_format == "bam" else None,
                              **dependency_identity},
            }
            job_id = content_id(identity)
            target = self.root / "jobs" / job_id
            with self._lock("job", job_id, deadline) as waited:
                metrics["lock_wait_seconds"] = waited
                if target.exists() or target.is_symlink():
                    phase = time.perf_counter()
                    manifest = self._read_job(job_id, identity, deadline)
                    metrics["output_verification_seconds"] = time.perf_counter() - phase
                    metrics.update({"native_alignment_seconds": 0.0, "conversion_seconds": 0.0,
                                    "output_validation_seconds": 0.0, "original_execution_metrics": manifest["execution_metrics"]})
                    self._assert_unchanged(pack_directory, pack_signatures)
                    metrics["request_wall_seconds"] = time.perf_counter() - start
                    return self._job_result(manifest, cache_hit=True, metrics=metrics)
                sam = stage / "alignment.sam"
                command = [str(binary), "-a", "-x", preset, "-t", str(threads), str(pack_directory / "index.mmi"), name]
                metrics["native_alignment_seconds"] = self._execute(
                    command, cwd=stage, stderr_path=stage / "alignment.stderr.log", stdout_path=sam, deadline=deadline,
                )
                self._assert_unchanged(pack_directory, pack_signatures)
                if binary_signature != file_signature(binary) or snapshot_signature != file_signature(snapshot):
                    raise IntegrityError("Executable or read snapshot changed during alignment")
                metrics["native_reported"] = _native_report(stage / "alignment.stderr.log")
                output = sam
                metrics["conversion_seconds"] = 0.0
                conversion_command = None
                if output_format == "bam":
                    output = stage / "alignment.bam"
                    conversion_command = [
                        sys.executable, "-m", "genocache.io", "convert", str(sam), str(output), str(threads), SORT_MEMORY,
                    ]
                    metrics["conversion_seconds"] = self._execute(
                        conversion_command, cwd=stage, stderr_path=stage / "conversion.stderr.log", deadline=deadline,
                    )
                phase = time.perf_counter()
                output_stats = self._validate_output(output, output_format, pack["input_stats"], stats, deadline)
                metrics["output_validation_seconds"] = time.perf_counter() - phase
                phase = time.perf_counter()
                files = {"alignment": self._file_record(output, deadline),
                         "native_stderr": self._file_record(stage / "alignment.stderr.log", deadline)}
                if output_format == "bam":
                    files.update({"index": self._file_record(stage / "alignment.bam.csi", deadline),
                                  "conversion_stderr": self._file_record(stage / "conversion.stderr.log", deadline)})
                    sam.unlink()
                metrics["output_verification_seconds"] = time.perf_counter() - phase
                # Exact input digests stay in the manifest; read data is not
                # retained as a second permanent copy in the job cache.
                snapshot.unlink()
                metrics["request_wall_seconds"] = time.perf_counter() - start
                manifest = {
                    "schema_version": SCHEMA_VERSION, "kind": "alignment_job", "job_id": job_id,
                    "identity": identity, "input_stats": stats, "output_stats": output_stats,
                    "files": files, "created_at": _utcnow(),
                    "execution_metrics": {**{key: value for key, value in metrics.items() if key != "request_wall_seconds"},
                                          "prepublish_wall_seconds": metrics["request_wall_seconds"]},
                    "command": ["minimap2", "-a", "-x", preset, "-t", str(threads),
                                f"packs/{reference_id}/index.mmi", name],
                    "conversion_command": ["python", "-m", "genocache.io", "convert", "alignment.sam", "alignment.bam",
                                           str(threads), SORT_MEMORY] if conversion_command else None,
                }
                write_json(stage / "manifest.json", manifest)
                check_deadline(deadline)
                self._publish(stage, target)
                metrics["request_wall_seconds"] = time.perf_counter() - start
                return self._job_result(manifest, cache_hit=False, metrics=metrics)
