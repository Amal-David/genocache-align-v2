from __future__ import annotations

from contextlib import contextmanager
from copy import deepcopy
import hashlib
import io
import json
import os
from pathlib import Path
import random
import shutil
import tarfile
import time

import pysam
import pytest

from genocache import cloud
from genocache.engine import Engine, IntegrityError


class FakeS3:
    """Only local bytes; these tests never load or contact AWS credentials."""

    def __init__(self, objects=None, *, fail_upload_number=None):
        self.objects = dict(objects or {})
        self.events = []
        self.fail_upload_number = fail_upload_number
        self.upload_count = 0

    def download_file(self, bucket, key, destination):
        self.events.append(("download", key))
        Path(destination).write_bytes(self.objects[(bucket, key)])

    def upload_file(self, source, bucket, key):
        self.upload_count += 1
        self.events.append(("upload", key))
        if self.upload_count == self.fail_upload_number:
            raise OSError("simulated transfer failure")
        self.objects[(bucket, key)] = Path(source).read_bytes()

    def put_object(self, *, Bucket, Key, Body, ContentType):
        assert ContentType == "application/json"
        self.events.append(("put", Key))
        self.objects[(Bucket, Key)] = Body

    def get_object(self, *, Bucket, Key):
        content = self.objects[(Bucket, Key)]
        return {"Body": io.BytesIO(content), "ContentLength": len(content)}


@pytest.fixture(scope="module")
def binary() -> Path:
    found = shutil.which(os.environ.get("MINIMAP2", "minimap2"))
    if found is None:
        local = Path(__file__).resolve().parents[3] / "minimap2-src" / "minimap2"
        found = str(local) if local.is_file() else None
    if found is None:
        pytest.skip("Real minimap2 is required for reference pack integration tests")
    return Path(found).resolve()


@pytest.fixture(scope="module")
def native_pack(tmp_path_factory, binary):
    work = tmp_path_factory.mktemp("cloud-native-source")
    rng = random.Random(9871)
    sequence = "".join(rng.choices("ACGT", k=24000))
    reference = work / "reference.fa"
    reference.write_text(f">chr1\n{sequence}\n")
    reads = work / "reads.fq"
    first = sequence[2100:5100]
    second = sequence[10000:13500].translate(str.maketrans("ACGT", "TGCA"))[::-1]
    reads.write_text(f"@first\n{first}\n+\n{'I' * len(first)}\n@second\n{second}\n+\n{'J' * len(second)}\n")
    engine = Engine(work / "source-cache", binary)
    pack = engine.build_index(reference)
    archive = work / "reference.tar"
    pointer = cloud.export_pack(engine, pack["reference_id"], archive)
    return {"engine": engine, "pack": pack, "archive": archive, "pointer": pointer, "reads": reads}


def _job(native_pack) -> dict:
    return {
        "schema": cloud.JOB_SCHEMA,
        "reference_pack": {"uri": "s3://test-bucket/input/reference.tar", "sha256": native_pack["pointer"]["sha256"],
                           "reference_id": native_pack["pack"]["reference_id"]},
        "reads": {"uri": "s3://test-bucket/input/reads.fq", "sha256": hashlib.sha256(native_pack["reads"].read_bytes()).hexdigest()},
        "output_prefix": "s3://test-bucket/output", "threads": 1, "timeout_seconds": 30, "output_format": "bam",
    }


def _client(native_pack, **kwargs) -> FakeS3:
    return FakeS3({
        ("test-bucket", "input/reference.tar"): native_pack["archive"].read_bytes(),
        ("test-bucket", "input/reads.fq"): native_pack["reads"].read_bytes(),
    }, **kwargs)


def _members(native_pack):
    with tarfile.open(native_pack["archive"], "r:") as archive:
        return [(entry.name, archive.extractfile(entry).read(), tarfile.REGTYPE) for entry in archive]


def _archive(path: Path, members):
    with tarfile.open(path, "w", format=tarfile.GNU_FORMAT) as archive:
        for name, body, kind in members:
            entry = tarfile.TarInfo(name)
            entry.type = kind
            entry.size = len(body) if kind == tarfile.REGTYPE else 0
            if kind in {tarfile.SYMTYPE, tarfile.LNKTYPE}:
                entry.linkname = "../../outside"
            archive.addfile(entry, io.BytesIO(body) if entry.size else None)
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_real_pack_export_import_relocation_and_alignment(native_pack, binary, tmp_path):
    target = Engine(tmp_path / "worker", binary)
    pointer = native_pack["pointer"]
    imported = cloud.import_pack(target, native_pack["archive"], pointer["reference_id"], pointer["sha256"])
    assert Path(imported["path"]).parent == target.root / "packs"
    assert all(path.stat().st_mode & 0o222 == 0 for path in Path(imported["path"]).iterdir())
    result = target.align(imported["reference_id"], native_pack["reads"])
    with pysam.AlignmentFile(result["output_path"], "rb") as stream:
        records = list(stream.fetch("chr1"))
        assert [record.query_name for record in records] == ["first", "second"]
        assert records[1].is_reverse
    assert not list((target.root / ".tmp").iterdir())
    before = Path(imported["index_path"]).stat().st_ino
    cloud.import_pack(target, native_pack["archive"], pointer["reference_id"], pointer["sha256"])
    assert Path(imported["index_path"]).stat().st_ino == before


def test_wrong_archive_sha_publishes_nothing(native_pack, binary, tmp_path):
    target = Engine(tmp_path / "worker", binary)
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        cloud.import_pack(target, native_pack["archive"], native_pack["pack"]["reference_id"], "0" * 64)
    assert not list((target.root / "packs").iterdir())
    assert not list((target.root / ".tmp").iterdir())


@pytest.mark.parametrize("attack", ["traversal", "dot_component", "symlink", "hardlink", "unexpected", "duplicate", "too_many", "sparse", "pax"])
def test_malicious_archives_are_rejected_without_publication(attack, native_pack, binary, tmp_path):
    target = Engine(tmp_path / "worker", binary)
    entries = _members(native_pack)
    reference_id = native_pack["pack"]["reference_id"]
    if attack == "traversal":
        entries[0] = (f"{reference_id}/../outside", b"bad", tarfile.REGTYPE)
    elif attack == "dot_component":
        name, data, kind = entries[0]
        entries[0] = (name.replace("/", "/./", 1), data, kind)
    elif attack in {"symlink", "hardlink", "sparse", "pax"}:
        kinds = {"symlink": tarfile.SYMTYPE, "hardlink": tarfile.LNKTYPE,
                 "sparse": tarfile.GNUTYPE_SPARSE, "pax": tarfile.XHDTYPE}
        entries[0] = (entries[0][0], b"", kinds[attack])
    elif attack == "unexpected":
        entries[0] = (f"{reference_id}/unexpected.file", b"bad", tarfile.REGTYPE)
    elif attack == "duplicate":
        entries[1] = entries[0]
    else:
        entries.append(entries[0])
    archive = tmp_path / "attack.tar"
    digest = _archive(archive, entries)
    with pytest.raises((ValueError, IntegrityError)):
        cloud.import_pack(target, archive, reference_id, digest)
    assert not list((target.root / "packs").iterdir())
    assert not list((target.root / ".tmp").iterdir())
    assert not (tmp_path / "outside").exists()


def test_forged_huge_extension_header_is_rejected_before_read(native_pack, binary, tmp_path):
    target = Engine(tmp_path / "worker", binary)
    header = tarfile.TarInfo("pax")
    header.type = tarfile.XHDTYPE
    header.size = 1 << 40
    archive = tmp_path / "huge-header.tar"
    archive.write_bytes(header.tobuf(format=tarfile.GNU_FORMAT) + b"\0" * 1024)
    with pytest.raises(ValueError, match="unsafe"):
        cloud.import_pack(target, archive, native_pack["pack"]["reference_id"], hashlib.sha256(archive.read_bytes()).hexdigest())
    assert not list((target.root / ".tmp").iterdir())


def test_truncated_or_trailing_payload_is_rejected(native_pack, binary, tmp_path):
    for case, contents in (
        ("truncated", native_pack["archive"].read_bytes()[:600]),
        ("trailing", native_pack["archive"].read_bytes() + b"hidden payload"),
    ):
        target = Engine(tmp_path / case, binary)
        archive = tmp_path / f"{case}.tar"
        archive.write_bytes(contents)
        with pytest.raises(ValueError):
            cloud.import_pack(target, archive, native_pack["pack"]["reference_id"], hashlib.sha256(contents).hexdigest())
        assert not list((target.root / "packs").iterdir())
        assert not list((target.root / ".tmp").iterdir())


def test_binary_mismatch_import_never_publishes(native_pack, binary, tmp_path):
    changed_binary = tmp_path / "different-minimap2"
    shutil.copyfile(binary, changed_binary)
    with changed_binary.open("ab") as output:
        output.write(b"\nchanged-binary-identity\n")
    changed_binary.chmod(0o755)
    target = Engine(tmp_path / "worker", changed_binary)
    with pytest.raises(IntegrityError, match="different minimap2"):
        cloud.import_pack(target, native_pack["archive"], native_pack["pack"]["reference_id"], native_pack["pointer"]["sha256"])
    assert not list((target.root / "packs").iterdir())
    assert not list((target.root / ".tmp").iterdir())


def test_failed_export_removes_partial_archive(native_pack, tmp_path, monkeypatch):
    original = tarfile.TarFile.add
    calls = 0

    def failing_add(self, *args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("simulated archive write failure")
        return original(self, *args, **kwargs)

    monkeypatch.setattr(tarfile.TarFile, "add", failing_add)
    archive = tmp_path / "partial.tar"
    with pytest.raises(OSError):
        cloud.export_pack(native_pack["engine"], native_pack["pack"]["reference_id"], archive)
    assert not archive.exists()


def test_cloud_output_manifest_is_last_and_local_reuse_is_explicit(native_pack, binary, tmp_path):
    job = _job(native_pack)
    client = _client(native_pack)
    worker = tmp_path / "worker"
    result = cloud.run_job(job, root=worker, minimap2=str(binary), client=client)
    assert result["reference_source"] == "verified_s3_archive"
    assert client.events[-1] == ("put", cloud.parse_s3(result["manifest_uri"])[1])
    stored = json.loads(client.objects[cloud.parse_s3(result["manifest_uri"])])
    assert stored["manifest_uri"] == result["manifest_uri"]
    assert "request_wall_seconds" not in stored
    assert result["request_wall_seconds"] >= stored["timings"]["wall_before_manifest_upload_seconds"]
    assert str(worker) not in json.dumps(stored)
    for entry in stored["files"].values():
        data = client.objects[cloud.parse_s3(entry["uri"])]
        assert hashlib.sha256(data).hexdigest() == entry["sha256"]
        assert len(data) == entry["byte_count"]

    # A verified local content identity is authoritative; no unused archive
    # URI/SHA claim is made when the worker does not download that archive.
    job["reference_pack"]["uri"] = "s3://test-bucket/missing.tar"
    job["reference_pack"]["sha256"] = "0" * 64
    client.events.clear()
    warmed = cloud.run_job(job, root=worker, minimap2=str(binary), client=client)
    assert warmed["reference_source"] == "verified_local_pack" and warmed["cache_hit"]
    assert [key for action, key in client.events if action == "download"] == ["input/reads.fq"]
    assert warmed["manifest_uri"] != result["manifest_uri"]


def test_failed_artifact_upload_has_no_success_manifest(native_pack, binary, tmp_path):
    client = _client(native_pack, fail_upload_number=2)
    with pytest.raises(OSError, match="transfer failure"):
        cloud.run_job(_job(native_pack), root=tmp_path / "worker", minimap2=str(binary), client=client)
    assert not any(action == "put" for action, _ in client.events)
    assert not any(key.endswith("result.json") for _, key in client.objects)


def test_wrong_reads_sha_rejected_and_download_partials_removed(native_pack, binary, tmp_path):
    client = _client(native_pack)
    destination = tmp_path / "reads.fq"
    with pytest.raises(ValueError, match="checksum"):
        cloud.download_verified(client, "s3://test-bucket/input/reads.fq", destination, "0" * 64)
    assert not destination.exists()
    assert not list(tmp_path.glob(".download-*"))
    job = _job(native_pack)
    job["reads"]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="checksum"):
        cloud.run_job(job, root=tmp_path / "worker", minimap2=str(binary), client=client)
    assert not any(action in {"upload", "put"} for action, _ in client.events)
    assert not list((tmp_path / "worker" / "jobs").iterdir())


@pytest.mark.parametrize("uri", [
    "https://test-bucket/file", "s3://test-bucket", "s3://user:secret@test-bucket/key",
    "s3://test-bucket/key?token=secret", "s3://test-bucket/key#secret", "s3://test-bucket:9000/key",
    "s3://test-bucket/../key", "s3://test-bucket/a/./b", "s3://test-bucket//key",
    "s3://test-bucket/a\\b", "s3://test-bucket/a\nb",
])
def test_ambiguous_or_credential_bearing_s3_uris_are_rejected(uri):
    with pytest.raises(ValueError):
        cloud.parse_s3(uri)


def test_job_schema_is_strict_and_freezes_input(native_pack):
    original = _job(native_pack)
    cases = []
    extra = deepcopy(original)
    extra["aws_access_key_id"] = "secret"
    cases.append(extra)
    nested = deepcopy(original)
    nested["reads"]["credentials"] = "secret"
    cases.append(nested)
    for field, value in (("threads", True), ("threads", 0), ("timeout_seconds", 1.5),
                         ("output_format", []), ("schema", "other")):
        case = deepcopy(original)
        case[field] = value
        cases.append(case)
    for case in cases:
        with pytest.raises(ValueError):
            cloud.validate_job(case)
    validated = cloud.validate_job(original)
    original["reads"]["uri"] = "s3://test-bucket/changed"
    assert validated["reads"]["uri"] == "s3://test-bucket/input/reads.fq"


def test_object_manifest_limit_and_duplicate_json_keys():
    for content in (b" " * (cloud.MAX_OBJECT_MANIFEST_BYTES + 1), b'{"schema":"a","schema":"b"}', b"[]"):
        client = FakeS3({("test-bucket", "job.json"): content})
        with pytest.raises(ValueError):
            cloud.load_s3_json(client, "s3://test-bucket/job.json")

    class ShortBody(io.BytesIO):
        def __init__(self, value):
            super().__init__(value)
            self.bytes_read = 0

        def read(self, requested=-1):
            assert 0 < requested <= 8192
            result = super().read(min(requested, 1024))
            self.bytes_read += len(result)
            return result

    body = ShortBody(b" " * (cloud.MAX_OBJECT_MANIFEST_BYTES + 100))
    with pytest.raises(ValueError, match="64 KiB"):
        cloud._read_manifest(body)
    assert body.bytes_read == cloud.MAX_OBJECT_MANIFEST_BYTES + 1


def test_cli_boto_error_is_sanitized_without_traceback(monkeypatch, capsys):
    botocore = pytest.importorskip("botocore.exceptions")
    secret = "AWS_SECRET_ACCESS_KEY=DO_NOT_PRINT_THIS"

    @contextmanager
    def unavailable():
        raise botocore.ClientError({"Error": {"Code": "AccessDenied", "Message": secret}}, "GetObject")
        yield  # pragma: no cover

    monkeypatch.setattr(cloud, "s3_client", unavailable)
    with pytest.raises(SystemExit) as raised:
        cloud.main(["run", "--manifest", "s3://test-bucket/job.json"])
    output = capsys.readouterr()
    assert raised.value.code == 2 and not output.out
    assert "AccessDenied" in output.err and secret not in output.err and "Traceback" not in output.err


def test_cloud_wall_metric_includes_transfer_directory_cleanup(native_pack, binary, tmp_path, monkeypatch):
    actual_temporary_directory = cloud.tempfile.TemporaryDirectory
    cleanup_done = []

    class MeasuredCleanup(actual_temporary_directory):
        def __exit__(self, *args):
            result = super().__exit__(*args)
            time.sleep(0.03)
            cleanup_done.append(time.perf_counter())
            return result

    monkeypatch.setattr(cloud.tempfile, "TemporaryDirectory", MeasuredCleanup)
    before = time.perf_counter()
    result = cloud.run_job(_job(native_pack), root=tmp_path / "worker", minimap2=str(binary), client=_client(native_pack))
    assert cleanup_done
    assert result["request_wall_seconds"] >= cleanup_done[-1] - before - 0.005
