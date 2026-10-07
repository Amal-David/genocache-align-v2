"""Authenticated Modal functions; no public web endpoint is created.

Run from production/: modal deploy deploy/modal_app.py
Requires a Modal secret named genocache-s3 with an appropriately scoped AWS
identity. Native packs must be built with this same container's minimap2 binary.
"""
from pathlib import Path

import modal


ROOT = Path(__file__).resolve().parents[1]
app = modal.App("genocache-production")
image = modal.Image.from_dockerfile(ROOT / "Dockerfile", context_dir=ROOT)


@app.function(
    image=image,
    secrets=[modal.Secret.from_name("genocache-s3")],
    cpu=8,
    memory=65536,
    ephemeral_disk=262144,
    timeout=86400,
    max_containers=2,
    scaledown_window=300,
)
def align(manifest_uri: str) -> dict:
    import boto3
    from genocache.cloud import load_s3_json, run_job

    client = boto3.client("s3")
    job = load_s3_json(client, manifest_uri)
    if job.get("threads", 0) > 8:
        raise ValueError("this deployment reserves eight CPUs; threads must be <= 8")
    return run_job(job, root=Path("/tmp/genocache-worker"), client=client)


@app.function(
    image=image,
    secrets=[modal.Secret.from_name("genocache-s3")],
    cpu=8,
    memory=65536,
    ephemeral_disk=262144,
    timeout=86400,
    max_containers=1,
)
def prepare_reference(reference_uri: str, reference_sha256: str, destination: str, preset: str = "map-ont") -> dict:
    import json
    import tempfile
    import boto3
    from genocache.cloud import download_verified, export_pack, parse_s3
    from genocache import Engine

    client = boto3.client("s3")
    with tempfile.TemporaryDirectory(prefix="genocache-reference-") as folder:
        folder = Path(folder)
        reference = folder / "reference.input"
        download_verified(client, reference_uri, reference, reference_sha256)
        engine = Engine(Path("/tmp/genocache-worker"))
        pack = engine.build_index(reference, preset=preset, threads=8, timeout=82800)
        archive = folder / "pack.tar"
        pointer = export_pack(engine, pack["reference_id"], archive)
        uri = destination.rstrip("/") + f'/{pointer["sha256"]}.tar'
        bucket, key = parse_s3(uri)
        client.upload_file(str(archive), bucket, key)
        pointer["uri"] = uri
        client.put_object(
            Bucket=bucket, Key=key + ".json",
            Body=(json.dumps(pointer, indent=2) + "\n").encode(),
            ContentType="application/json",
        )
        return pointer


@app.local_entrypoint()
def main(manifest_uri: str):
    import json
    print(json.dumps(align.remote(manifest_uri), indent=2))
