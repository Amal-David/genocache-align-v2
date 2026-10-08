"""Small S3/Modal/AWS Batch bridge. Credentials come from the worker's identity.

Reference packs are immutable tar files with a required whole-file SHA-256.
All locks and native I/O stay on local disk. Large reads/results never traverse
an agent prompt or Modal function return value.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import json
import os
from pathlib import Path
import re
import shutil
import tarfile
import tempfile
import time
from urllib.parse import urlsplit
import uuid

from .engine import Engine
from .io import EngineError, IntegrityError, file_signature, hash_file


SHA = re.compile(r"^[0-9a-f]{64}$")
PACK_SCHEMA = "genocache.cloud-pack.v1"
JOB_SCHEMA = "genocache.cloud-job.v1"
MAX_OBJECT_MANIFEST_BYTES = 65536
MAX_PACK_MANIFEST_BYTES = 64 * 1024 * 1024
MAX_PACK_MEMBERS = 4
COPY_CHUNK_BYTES = 8 * 1024 * 1024
BUCKET = re.compile(r"^[a-z0-9][a-z0-9.-]{1,61}[a-z0-9]$")


def parse_s3(uri: str) -> tuple[str, str]:
    if not isinstance(uri, str) or any(ord(c) < 33 or ord(c) == 127 for c in uri):
        raise ValueError("S3 URI must be a string without whitespace or control characters")
    try:
        value = urlsplit(uri)
    except ValueError as exc:
        raise ValueError("invalid S3 URI") from exc
    key = value.path[1:] if value.path.startswith("/") else ""
    if (
        value.scheme != "s3" or not BUCKET.fullmatch(value.netloc) or not key
        or value.query or value.fragment or value.username or value.password
        or ".." in value.netloc or ".-" in value.netloc or "-." in value.netloc
        or "\\" in key or any(part in {"", ".", ".."} for part in key.rstrip("/").split("/"))
        or key.endswith("//")
    ):
        raise ValueError("expected s3://bucket/key with a valid bucket and unambiguous path, without credentials, query or fragment")
    return value.netloc, key


def require_sha(digest: str) -> str:
    if not isinstance(digest, str) or not SHA.fullmatch(digest):
        raise ValueError("a lowercase SHA-256 is required for each immutable input")
    return digest


def download_verified(client, uri: str, destination: Path, digest: str) -> None:
    require_sha(digest)
    bucket, key = parse_s3(uri)
    destination = Path(destination)
    if destination.exists() or destination.is_symlink():
        raise FileExistsError("Download destination already exists")
    descriptor, temporary_name = tempfile.mkstemp(prefix=".download-", dir=destination.parent)
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        client.download_file(bucket, key, str(temporary))
        if hash_file(temporary)["sha256"] != digest:
            raise ValueError("Downloaded object checksum does not match the job manifest")
        # Publish without overwriting an existing destination, including a
        # destination created concurrently after the initial existence check.
        os.link(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("Duplicate keys are not allowed in a manifest")
        result[key] = value
    return result


def _read_manifest(body) -> dict:
    data = bytearray()
    while len(data) <= MAX_OBJECT_MANIFEST_BYTES:
        chunk = body.read(min(8192, MAX_OBJECT_MANIFEST_BYTES + 1 - len(data)))
        if not chunk:
            break
        data.extend(chunk)
    if len(data) > MAX_OBJECT_MANIFEST_BYTES:
        raise ValueError("job/pack manifest exceeds 64 KiB")
    value = json.loads(data, object_pairs_hook=_unique_object)
    if not isinstance(value, dict):
        raise ValueError("manifest must be a JSON object")
    return value


def load_s3_json(client, uri: str) -> dict:
    bucket, key = parse_s3(uri)
    response = client.get_object(Bucket=bucket, Key=key)
    with response["Body"] as body:
        if response.get("ContentLength", 0) > MAX_OBJECT_MANIFEST_BYTES:
            raise ValueError("job/pack manifest exceeds 64 KiB")
        return _read_manifest(body)


def _upload_json(client, uri: str, value: dict) -> None:
    bucket, key = parse_s3(uri)
    body = (json.dumps(value, indent=2, allow_nan=False) + "\n").encode()
    if len(body) > MAX_OBJECT_MANIFEST_BYTES:
        raise ValueError("job/pack/result manifest exceeds 64 KiB")
    client.put_object(
        Bucket=bucket, Key=key,
        Body=body,
        ContentType="application/json",
    )


def export_pack(engine: Engine, reference_id: str, archive: Path) -> dict:
    """Archive only verified fixed-name pack files; no links or arbitrary trees."""
    pack = engine.inspect_reference(reference_id)
    path = Path(pack["path"])
    manifest = json.loads((path / "manifest.json").read_text())
    names = ["manifest.json"] + sorted(v["path"] for v in manifest["files"].values())
    archive = Path(archive)
    if archive.exists() or archive.is_symlink():
        raise FileExistsError("Archive destination already exists")
    if path == archive.parent.resolve() or path in archive.resolve().parents:
        raise ValueError("archive must be outside the immutable pack")
    before = {name: file_signature(path / name) for name in names}
    created = False
    try:
        with archive.open("xb") as output:
            created = True
            # GNU numeric fields support >8 GiB references/indexes without PAX
            # metadata. Fixed short member names need no extended headers.
            with tarfile.open(fileobj=output, mode="w", format=tarfile.GNU_FORMAT, dereference=True) as tar:
                for name in names:
                    source = path / name
                    if source.is_symlink() or not source.is_file() or Path(name).name != name:
                        raise ValueError("pack contains an unsafe member")
                    tar.add(source, arcname=f"{reference_id}/{name}", recursive=False)
            output.flush()
            os.fsync(output.fileno())
        if any(file_signature(path / name) != signature for name, signature in before.items()):
            raise IntegrityError("Reference pack changed during archive export")
        engine.inspect_reference(reference_id)
        archive_identity = hash_file(archive)
    except BaseException:
        if created:
            archive.unlink(missing_ok=True)
        raise
    return {
        "schema": PACK_SCHEMA, "reference_id": reference_id,
        **archive_identity,
        "preset": pack["preset"],
    }


def _extract_regular_pack(archive: Path, stage: Path, reference_id: str) -> None:
    """Read at most four fixed regular members, with bounded header/data reads.

    Generic TarFile iteration processes PAX/GNU metadata before returning a
    member. Parsing one 512-byte header first lets us reject those extensions,
    sparse files, links, and arbitrary member floods without allocating them.
    """
    allowed = {"manifest.json", "reference.fa", "reference.fa.gz", "index.mmi", "index.stderr.log"}
    seen = set()
    archive_size = archive.stat().st_size
    with archive.open("rb") as source:
        while True:
            header = source.read(tarfile.BLOCKSIZE)
            if len(header) != tarfile.BLOCKSIZE:
                raise ValueError("reference archive is missing its complete end marker")
            if header == b"\0" * tarfile.BLOCKSIZE:
                if source.read(tarfile.BLOCKSIZE) != b"\0" * tarfile.BLOCKSIZE:
                    raise ValueError("reference archive has an invalid end marker")
                while padding := source.read(COPY_CHUNK_BYTES):
                    if any(padding):
                        raise ValueError("unexpected data after the reference archive end marker")
                break
            if len(seen) >= MAX_PACK_MEMBERS:
                raise ValueError("unexpected pack member count")
            try:
                member = tarfile.TarInfo.frombuf(header, "utf-8", "strict")
            except (tarfile.HeaderError, UnicodeError) as exc:
                raise ValueError("invalid reference archive header") from exc
            filename = member.name.rsplit("/", 1)[-1]
            if (
                member.type not in {tarfile.REGTYPE, tarfile.AREGTYPE}
                or filename not in allowed or member.name != f"{reference_id}/{filename}"
                or filename in seen or member.size < 0
                or member.size > archive_size - source.tell()
            ):
                raise ValueError("unsafe, duplicate, or truncated reference archive member")
            if filename == "manifest.json" and member.size > MAX_PACK_MANIFEST_BYTES:
                raise ValueError("reference pack manifest exceeds 64 MiB")
            seen.add(filename)
            remaining = member.size
            with (stage / filename).open("xb") as output:
                while remaining:
                    chunk = source.read(min(COPY_CHUNK_BYTES, remaining))
                    if not chunk:
                        raise ValueError("truncated reference archive member")
                    output.write(chunk)
                    remaining -= len(chunk)
            padding_size = (-member.size) % tarfile.BLOCKSIZE
            if len(source.read(padding_size)) != padding_size:
                raise ValueError("truncated reference archive padding")
    if len(seen) != MAX_PACK_MEMBERS or seen - {"reference.fa", "reference.fa.gz"} != {"manifest.json", "index.mmi", "index.stderr.log"}:
        raise ValueError("reference archive inventory is invalid")


def import_pack(engine: Engine, archive: Path, reference_id: str, sha256: str) -> dict:
    """Reject traversal, links, duplicate members and unexpected filenames."""
    require_sha(reference_id)
    require_sha(sha256)
    archive = Path(archive)
    archive_signature = file_signature(archive)
    if hash_file(archive)["sha256"] != sha256:
        raise ValueError("reference archive SHA-256 mismatch")
    with engine._lock("pack", reference_id, None):
        target = engine.root / "packs" / reference_id
        if target.exists():
            return engine.inspect_reference(reference_id)
        stage = Path(tempfile.mkdtemp(prefix="cloud-pack-", dir=engine.root / ".tmp"))
        try:
            _extract_regular_pack(archive, stage, reference_id)
            if archive_signature != file_signature(archive):
                raise IntegrityError("Reference archive changed during import")
            # Validate before publication by using an isolated engine root.
            validation_root = stage.parent / f"validate-{uuid.uuid4().hex}"
            try:
                isolated = Engine(validation_root, engine.minimap2)
                validation_target = validation_root / "packs" / reference_id
                os.rename(stage, validation_target)
                isolated.inspect_reference(reference_id)
                stage = validation_target
                engine._publish(stage, target)
            finally:
                shutil.rmtree(validation_root, ignore_errors=True)
            return engine.inspect_reference(reference_id)
        finally:
            shutil.rmtree(stage, ignore_errors=True)


def validate_job(job: dict) -> dict:
    fields = {"schema", "reference_pack", "reads", "output_prefix", "threads", "timeout_seconds", "output_format"}
    if not isinstance(job, dict) or set(job) != fields or job["schema"] != JOB_SCHEMA:
        raise ValueError("cloud job must match the documented v1 schema exactly")
    pack, reads = job["reference_pack"], job["reads"]
    if not isinstance(pack, dict) or set(pack) != {"uri", "sha256", "reference_id"}:
        raise ValueError("reference_pack requires uri, sha256 and reference_id")
    if not isinstance(reads, dict) or set(reads) != {"uri", "sha256"}:
        raise ValueError("reads requires uri and sha256")
    for item in (pack, reads):
        parse_s3(item["uri"])
        require_sha(item["sha256"])
    require_sha(pack["reference_id"])
    parse_s3(job["output_prefix"])
    if type(job["threads"]) is not int or not 1 <= job["threads"] <= 256:
        raise ValueError("threads must be 1..256")
    if type(job["timeout_seconds"]) is not int or not 1 <= job["timeout_seconds"] <= 86400:
        raise ValueError("timeout_seconds must be 1..86400")
    if not isinstance(job["output_format"], str) or job["output_format"] not in {"bam", "sam"}:
        raise ValueError("output_format must be bam or sam")
    # A caller mutating its input dictionary must not change a validated job
    # while the worker is downloading or uploading its artifacts.
    return {**job, "reference_pack": dict(pack), "reads": dict(reads)}


def run_job(job: dict, *, root: Path, minimap2: str = "minimap2", client=None) -> dict:
    started = time.perf_counter()
    job = validate_job(job)
    if client is None:
        import boto3
        client = boto3.client("s3")
    engine = Engine(root, minimap2)
    worker_setup_seconds = time.perf_counter() - started
    pack = job["reference_pack"]
    phase = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix="genocache-transfer-") as temporary:
        work = Path(temporary)
        if (engine.root / "packs" / pack["reference_id"]).exists():
            engine.inspect_reference(pack["reference_id"])
            reference_source = "verified_local_pack"
        else:
            archive = work / "reference.tar"
            download_verified(client, pack["uri"], archive, pack["sha256"])
            import_pack(engine, archive, pack["reference_id"], pack["sha256"])
            archive.unlink()
            reference_source = "verified_s3_archive"
        reference_stage_seconds = time.perf_counter() - phase
        phase = time.perf_counter()
        reads = work / "reads.input"
        download_verified(client, job["reads"]["uri"], reads, job["reads"]["sha256"])
        reads_stage_seconds = time.perf_counter() - phase
        result = engine.align(
            pack["reference_id"], reads, threads=job["threads"],
            timeout=job["timeout_seconds"], output_format=job["output_format"],
        )
        phase = time.perf_counter()
        # Attempts have separate prefixes so simultaneous retries cannot mix artifacts.
        prefix = job["output_prefix"].rstrip("/") + f'/{result["job_id"]}/{uuid.uuid4().hex}'
        metadata = json.loads(Path(result["manifest_path"]).read_text())
        uploaded = {}
        for label, entry in metadata["files"].items():
            path = Path(result["job_path"]) / entry["path"]
            before = file_signature(path)
            identity = hash_file(path)
            if identity["sha256"] != entry["sha256"] or identity["byte_count"] != entry["byte_count"]:
                raise IntegrityError("Cached result changed before upload")
            uri = prefix + "/" + path.name
            bucket, key = parse_s3(uri)
            client.upload_file(str(path), bucket, key)
            if before != file_signature(path):
                raise IntegrityError("Cached result changed during upload")
            uploaded[label] = {"uri": uri, "sha256": entry["sha256"], "byte_count": entry["byte_count"]}
        # Strip local absolute paths from the cloud-facing return value.
        final = {
            "schema": "genocache.cloud-result.v1",
            "job_id": result["job_id"], "reference_id": result["reference_id"],
            "cache_hit": result["cache_hit"], "output_format": result["output_format"],
            "reference_source": reference_source,
            "input_stats": result["input_stats"], "output_stats": result["output_stats"],
            "files": uploaded, "engine_metrics": result["metrics"],
            "timings": {
                "worker_setup_seconds": worker_setup_seconds,
                "reference_stage_seconds": reference_stage_seconds,
                "reads_stage_seconds": reads_stage_seconds,
                "upload_seconds": time.perf_counter() - phase,
                "wall_before_manifest_upload_seconds": time.perf_counter() - started,
            },
            "provenance": metadata,
        }
        manifest_uri = prefix + "/result.json"
        final["manifest_uri"] = manifest_uri
        phase = time.perf_counter()
        _upload_json(client, manifest_uri, final)  # Commit marker comes last.
        # These completion measurements are returned to the caller. The stored
        # manifest accurately labels its earlier pre-upload wall boundary.
        final["manifest_upload_seconds"] = time.perf_counter() - phase
    final["request_wall_seconds"] = time.perf_counter() - started
    return final


@contextmanager
def s3_client():
    import boto3
    client = boto3.client("s3")
    try:
        yield client
    finally:
        client.close()


def _safe_cli_error(exc: Exception) -> str:
    """SDK error messages may contain signed URLs or sensitive request data."""
    try:
        from botocore.exceptions import BotoCoreError, ClientError
    except ImportError:
        BotoCoreError = ClientError = ()
    if isinstance(exc, ClientError):
        code = exc.response.get("Error", {}).get("Code", "Unknown")
        safe_code = code if isinstance(code, str) and re.fullmatch(r"[A-Za-z0-9_.-]{1,80}", code) else "Unknown"
        return f"S3 request failed ({safe_code}); check the worker identity and configured object access"
    if isinstance(exc, BotoCoreError):
        return f"S3 client failed ({type(exc).__name__}); check worker credentials, region and connectivity"
    if isinstance(exc, ImportError):
        return "Cloud dependencies are unavailable; install the package's cloud extra"
    if isinstance(exc, (ValueError, EngineError, tarfile.TarError)):
        return str(exc)[:2048]
    if isinstance(exc, OSError):
        return f"Local filesystem operation failed: {exc.strerror or type(exc).__name__}"
    return f"Cloud worker failed ({type(exc).__name__}); error details were omitted"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(os.environ.get("GENOCACHE_ROOT", "/tmp/genocache-worker")))
    parser.add_argument("--minimap2", default=os.environ.get("MINIMAP2", "minimap2"))
    subs = parser.add_subparsers(dest="command", required=True)
    export = subs.add_parser("export-pack")
    export.add_argument("--reference-id", required=True)
    export.add_argument("--destination", required=True, help="S3 prefix")
    run = subs.add_parser("run")
    run.add_argument("--manifest", required=True, help="Local JSON file or s3://bucket/job.json")
    args = parser.parse_args(argv)
    try:
        with s3_client() as client:
            if args.command == "export-pack":
                engine = Engine(args.root, args.minimap2)
                parse_s3(args.destination)
                with tempfile.TemporaryDirectory(prefix="genocache-export-") as directory:
                    archive = Path(directory) / "pack.tar"
                    pointer = export_pack(engine, args.reference_id, archive)
                    uri = args.destination.rstrip("/") + f'/{pointer["sha256"]}.tar'
                    bucket, key = parse_s3(uri)
                    client.upload_file(str(archive), bucket, key)
                    pointer["uri"] = uri
                    _upload_json(client, uri + ".json", pointer)
                print(json.dumps(pointer, indent=2))
            else:
                if args.manifest.startswith("s3://"):
                    job = load_s3_json(client, args.manifest)
                else:
                    with Path(args.manifest).open("rb") as body:
                        job = _read_manifest(body)
                result = run_job(job, root=args.root, minimap2=args.minimap2, client=client)
                print(json.dumps(result, indent=2, allow_nan=False))
        return 0
    except Exception as exc:
        parser.exit(2, f"genocache-cloud: {_safe_cli_error(exc)}\n")


if __name__ == "__main__":
    raise SystemExit(main())
