# Deploy the GenoCache batch worker

Updated 2026-10-07. Commands run from `production/`. The worker executes native
minimap2 and returns BAM/SAM artifacts plus provenance. GPU embedding experiments
are a separate optional workload; this CPU worker does not require Torch or
model weights.

Local tests exercised real native mapping, reference-pack relocation, checksum
failures, retry isolation and S3 publication order. Modal 1.6.1 imports this app
definition, and botocore validates the Batch job-definition schema. Docker image
construction and live cloud execution were not available in the development
environment. Build the image and complete the smoke job below before submitting
a large sample.

## 1. Build and check the image

```bash
docker build --platform linux/amd64 -t genocache:0.1.0 .
docker run --rm genocache:0.1.0 genocache doctor
```

The Docker build context must be `production/`. The image contains only the
installed `genocache` package and pinned direct runtime dependencies. It does
not contain repository examples, benchmark scripts, the optional Torch adapter,
or genomic inputs. Native minimap2 is compiled from commit
`3c28777e7e2dcc90f825de1b9f17a89cca7d4452` (v2.31).

After pushing the image to your registry, record its immutable digest in the
deployment and experiment receipt. Reference preparation and mapping must use
the **same compiled binary**. A locally compiled v2.31 executable is not assumed
byte-identical to the container's v2.31 executable. Packs are rejected when the
binary SHA differs. Rebuild packs inside the selected release image.

The default Docker user is uid 10001. For a bind-mounted local smoke run, use
your host uid/gid or give that uid write access to the intended directory:

```bash
mkdir -p local-data/container-work
docker run --rm --user "$(id -u):$(id -g)" \
  -v "$PWD/local-data/container-work:/work" \
  -v /data/smoke:/data:ro genocache:0.1.0 \
  genocache --root /work/cache index \
  --reference /data/reference.fa --preset map-ont --threads 2 \
  > local-data/container-reference.json
```

Read `reference_id` from that JSON and invoke `genocache align` with the same
mounts, image and cache directory. A second identical invocation should report
`cache_hit: true`; verify that output counts and artifact hashes agree. Use a
small independently characterized sample before scaling the deployment.

## 2. S3 inputs and completed results

The worker consumes a strict JSON manifest; see
[`examples/cloud-job.json`](../examples/cloud-job.json). Replace every placeholder:

```json
{
  "schema": "genocache.cloud-job.v1",
  "reference_pack": {
    "uri": "s3://YOUR_BUCKET/reference-packs/ARCHIVE_DIGEST.tar",
    "sha256": "64_lowercase_hex_characters",
    "reference_id": "reference_id_returned_by_preparation"
  },
  "reads": {
    "uri": "s3://YOUR_BUCKET/inputs/reads.fastq.gz",
    "sha256": "64_lowercase_hex_characters"
  },
  "output_prefix": "s3://YOUR_BUCKET/outputs",
  "threads": 8,
  "timeout_seconds": 7200,
  "output_format": "bam"
}
```

SHA-256 is over the exact uploaded bytes. Gzip and uncompressed versions of the
same biological sequence have different digests. A reference pack contains the
reference snapshot, native index, manifest and index-build log. Export through
`genocache-cloud export-pack` or `genocache.cloud.export_pack`; arbitrary tar
archives and legacy index folders do not satisfy the import contract.

Each completed attempt gets its own prefix:

```text
s3://YOUR_BUCKET/outputs/JOB_ID/ATTEMPT_ID/result.json
```

Read `files` inside the final manifest for exact BAM/SAM, CSI and log URIs,
checksums and sizes. Do not infer completion from an existing BAM object: the
worker uploads `result.json` last. A failed upload leaves no success manifest.
The same job can have multiple completed attempts; consumers should select one
complete attempt by `job_id` and their orchestration policy.

The returned cloud result also includes `reference_source`, distinguishing a
fully verified local pack from a newly downloaded verified archive. The archive
digest is checked on download; a matching local pack is verified by its own
identity and contents. Cloud wall time includes staging, mapping, upload and
temporary cleanup. The stored manifest labels its earlier pre-manifest-upload
timing boundary explicitly.

## 3. Modal

Use an authenticated Modal workspace and a Secret named `genocache-s3` containing
an AWS identity and region accepted by boto3. Give it read access to the selected
job/reference/read prefixes and write access to output/reference-pack prefixes.
Configure credentials in the [Modal Secrets panel](https://modal.com/secrets)
or your existing secret-management workflow; no credentials belong in the job
JSON or source tree. See [Modal Secrets](https://modal.com/docs/guide/secrets).

Deploy the two authenticated functions:

```bash
modal deploy deploy/modal_app.py
```

No public HTTP endpoint is created. The app defines `prepare_reference` and
`align` using the same Dockerfile image. The initial allocation is eight CPUs,
64 GiB RAM and 256 GiB ephemeral disk, with at most two alignment containers.
It is a pilot configuration, not a guarantee that a whole-genome sample fits.
Modal ignores Docker's `USER` instruction and runs these functions as root;
see [Modal image compatibility](https://modal.com/docs/guide/existing-images).

Prepare the reference in that exact deployed image:

```python
import json
from pathlib import Path
import modal

prepare = modal.Function.from_name("genocache-production", "prepare_reference")
pointer = prepare.remote(
    "s3://YOUR_BUCKET/inputs/reference.fa",
    "REPLACE_WITH_REFERENCE_SHA256",
    "s3://YOUR_BUCKET/reference-packs",
    preset="map-ont",
)
Path("local-data/reference-pointer.json").write_text(json.dumps(pointer, indent=2))
```

Copy the returned `uri`, `sha256` and `reference_id` into the cloud job manifest.
Upload that manifest to its chosen S3 jobs prefix, then invoke the deployed
function from Codex, a scheduler, or Python:

```python
import json
from pathlib import Path
import modal

align = modal.Function.from_name("genocache-production", "align")
result = align.remote("s3://YOUR_BUCKET/jobs/smoke.json")
Path("local-data/cloud-result.json").write_text(json.dumps(result, indent=2))
print(result["manifest_uri"])
```

For a development invocation that builds/uses the local app definition:

```bash
modal run deploy/modal_app.py --manifest-uri s3://YOUR_BUCKET/jobs/smoke.json
```

The native cache lives in `/tmp/genocache-worker` on each container. A warm
container can reuse it; a replacement container will stage the immutable pack
again. There is no promise of a cache hit across containers or redeployments.
The worker deliberately uses local POSIX locks. Do not put its root on a Modal
Volume or network filesystem; [Modal Volumes](https://modal.com/docs/guide/volumes)
do not provide the shared locking semantics this cache requires.

## 4. AWS Batch on EC2

Use an existing x86 EC2 Batch compute environment/queue, a private ECR repository,
and an IAM job role scoped to the selected S3 objects. The ECS/node configuration
also needs its normal ECR pull and CloudWatch log permissions. These are account
resources supplied by the deployment; the package does not invent account IDs,
create a network, or embed access keys.

Push the image from step 1 to your ECR repository and use its digest. Edit
[`deploy/aws-batch-job-definition.json`](../deploy/aws-batch-job-definition.json):

- Replace account, region, image digest, job-role ARN and default manifest URI.
- The example reserves 8 vCPUs and 61,440 MiB, leaving memory headroom on a nominal
  64 GiB machine. Select a larger node if the reference/sample demands it.
- Provision **container-visible** temporary storage. This definition writes to
  `/tmp` and has no host bind mount. An EC2 NVMe disk is not automatically the
  backing store of Docker `/tmp`; configure sufficient ECS/Docker storage.
  If adding a dedicated local-disk bind mount, set both the command's `--root`
  and `TMPDIR` onto that mounted disk with suitable write permissions.
- `threads` in the submitted manifest must fit the reserved CPUs. Native I/O,
  sorting and SDK threads add some overhead; benchmark the chosen allocation.

Register and submit:

```bash
aws batch register-job-definition \
  --cli-input-json file://deploy/aws-batch-job-definition.json
aws batch submit-job \
  --job-name genocache-smoke \
  --job-queue YOUR_EXISTING_QUEUE \
  --job-definition genocache-production \
  --parameters manifest=s3://YOUR_BUCKET/jobs/smoke.json
```

Use a pinned job-definition revision for repeatable subsequent runs. AWS applies
the `Ref::manifest` command parameter at submission; see its
[job-definition documentation](https://docs.aws.amazon.com/batch/latest/userguide/job_definition_parameters.html).

Prepare reference packs on an authorized machine/container using the exact ECR
image and the local `genocache index` command from step 1. With that same engine
root and an AWS identity available through the existing SDK credential chain,
export them through:

```bash
genocache-cloud --root /worker-local/cache export-pack \
  --reference-id RETURNED_REFERENCE_ID \
  --destination s3://YOUR_BUCKET/reference-packs
```

Run that command **inside the release image** or an environment with its exact
native binary. If preparing in Modal and executing in independently rebuilt AWS
images, verify the binary hashes first; equal source versions are insufficient.
For cross-platform reuse, pull the same immutable image rather than recompiling
two copies. Keep native indexing as a one-time preparation job per reference,
preset and binary identity.

## 5. Resource and reliability limits

Measure peak RAM, temporary disk, input/output transfer, first-run latency and
warm throughput on the target sample. Reference-pack extraction and verification
temporarily duplicate some files; mapping retains a SAM while sorting to BAM.
Compressed read size alone is not an adequate disk estimate. A dense vector
index additionally retains full vectors for exact auditing, even when an IVF-PQ
index is used; account for both storage and memory during construction.

The default worker is a single-job process/function. Local atomic publication
and locks handle concurrent requests on one local filesystem; workers do not
coordinate through distributed locks. AWS retries or cold Modal containers can
repeat work. Output prefixes isolate attempts and the final manifest identifies
a completed result. There is no automatic local-cache eviction; monitor storage
and retire idle workers or remove completed caches under an orchestration policy.

Choose a batch size that fits the allocation and preserves the intended whole
read records. Future cross-worker result lookup, resident native indexes and
hot-region admission need their own measured implementation. They are not
present behavior of this release.

Before a large run, the cloud smoke job should establish: image starts; native
pack imports under the matching binary; S3 permissions work; counts and sampled
records agree with the native baseline; BAM/CSI download and open successfully;
and a completed manifest is present. Record the image digest, job-definition
revision or Modal deployment version, input digests, wall time and billed
resources alongside that result. Then use [CODEX_RUNBOOK.md](CODEX_RUNBOOK.md) to
evaluate embedding candidates against the same workload.
