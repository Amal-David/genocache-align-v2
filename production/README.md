# GenoCache worker and embedding retrieval

GenoCache now provides a runnable batch alignment tool with native minimap2
output, content-addressed reference packs, verified exact-input result caching,
S3 job manifests, and Modal/AWS deployment definitions. A separate embedding
path turns an explicitly supplied encoder into batched reference/query vectors
and position-aware candidate intervals, with truth-based evaluation.

**Every novel alignment batch currently runs through minimap2.** The embedding
path is an experimental candidate generator. It does not emit SAM, assign MAPQ,
or bypass native mapping. No novel-read speed advantage over minimap2 or
Parabricks is claimed.

## What is available

| Component | Implemented behavior |
| --- | --- |
| `genocache` | Build/inspect immutable reference packs; align FASTA/FASTQ, including gzip; preserve native SAM or create sorted BAM + CSI; structured JSON results |
| Python `Engine` | Same functions for a worker, agent, scheduler, or notebook |
| Exact job cache | Full input bytes, reference identity, native binary and configuration determine a job; verified hits reuse output |
| `genocache-cloud` | Verify S3 inputs, safely import reference packs, execute jobs locally, upload artifacts before the result manifest |
| Modal / AWS Batch | Authenticated Modal functions and an EC2 Batch job definition for the same worker |
| Encoder adapter | Trusted TorchScript `[B,4,L] -> [B,D]`, streaming reference windows, both query orientations, bounded batches, hashes and timings |
| `genocache-retrieval` | Explicit normalized inner product; exact Flat or IVF-PQ; compact positions; batch search; optional exhaustive-search audit |
| Candidate regions | Interval uncertainty, one vote per query interval, distinct loci and strands, explicit truncation diagnostics |
| `genocache-evaluate` | Origin/strand evaluation, full-span candidate recall, missing/unmapped denominators, per-stratum uncertainty, MAPQ diagnostics |
| `genocache-model` | Index memory, window coverage, Amdahl/routing bounds, error budgets, cache break-even, confidence intervals |

The native engine uses disk-resident `.mmi` packs loaded once per uncached batch.
It does not keep a native index resident across separate jobs. The vector index
can be loaded once and reused through its Python object; separate CLI processes
pay their own load/verification cost. These costs appear in measurements.

## Local installation

Use Linux, Python **3.12 or later**, and local POSIX storage. Commands in this
file run from `production/` unless stated otherwise.

```bash
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-runtime.txt
python -m pip install -e '.[test]'
```

Install minimap2 v2.31 from the pinned source commit:

```bash
git clone https://github.com/lh3/minimap2.git /tmp/genocache-minimap2
git -C /tmp/genocache-minimap2 checkout 3c28777e7e2dcc90f825de1b9f17a89cca7d4452
make -C /tmp/genocache-minimap2 -j4
export MINIMAP2=/tmp/genocache-minimap2/minimap2
genocache --root ./local-data/cache doctor
```

The build requires a C compiler, make and zlib development headers. The supplied
Dockerfile performs that build in its first stage. Use the same binary artifact
to prepare and consume a reference pack: the engine checks its SHA-256, not just
the displayed version.

## Align a batch now

```bash
mkdir -p local-data
genocache --root ./local-data/cache index \
  --reference /data/reference.fa --preset map-ont --threads 8 \
  --timeout 7200 > local-data/reference.json
```

Choose the preset for the reads. Supported DNA presets are `map-ont`, `lr:hq`,
`map-hifi`, and `map-pb`. See the [minimap2 manual](https://lh3.github.io/minimap2/minimap2.html).
Spliced RNA mapping and paired short-read libraries are outside this release.

The JSON contains a `reference_id`. Run mapping with it:

```bash
REFERENCE_ID=$(python -c 'import json; print(json.load(open("local-data/reference.json"))["reference_id"])')
genocache --root ./local-data/cache align \
  --reference-id "$REFERENCE_ID" --reads /data/reads.fastq.gz \
  --threads 8 --timeout 7200 --output-format bam > local-data/job.json
```

`job.json` contains `job_id`, `output_path`, `index_path`, `manifest_path`, input
and output counts, measured stages, and `cache_hit`. Repeating the identical
request verifies and reuses that job. Changed read names, qualities, compression,
ordering, mapper binary or relevant configuration intentionally miss the cache.
The worker checks the full input and output artifacts; a hit is not zero-I/O.

Python integration:

```python
from pathlib import Path
from genocache import Engine

engine = Engine(Path("/worker-local/genocache"), minimap2="minimap2")
pack = engine.build_index(Path("/data/reference.fa"), preset="map-ont", threads=8)
job = engine.align(pack["reference_id"], Path("/data/reads.fastq.gz"),
                   threads=8, timeout=7200, output_format="bam")
print(job["output_path"])
```

Timeouts govern native execution and associated waits; they do not impose a
deadline on every filesystem byte operation. Cloud job limits separately bound
the overall container/function lifetime. Check available local disk before
large runs because snapshots, SAM, sorting temporaries and BAM coexist.

## Embeddings and the proposed change in approach

The next hypothesis is **coarse genomic localization followed by verified local
work**, with global fallback and explicit ambiguity. It is not enough to replace
every tiny exact seed lookup with a neural call. A dense embedding index can
increase memory and search work, while a wrong candidate caps final accuracy.

The implemented retrieval path makes that hypothesis measurable:

1. Supply a model with a verified architecture and preprocessing contract.
2. Encode reference windows once and query chunks in batches, preserving both
   query orientations and original coordinates.
3. Establish true-locus recall with exact vector search before compressing or
   approximating the index.
4. Convert retrieved windows into compatible diagonal **intervals**, retaining
   distinct loci and independent query support.
5. Evaluate complete candidate containment and difficult strata alongside the
   native mapping result.

Follow [ENCODER_CONTRACT.md](docs/ENCODER_CONTRACT.md) for executable vector and
search commands. No weights are shipped. Quality-dependent, tokenized, or other
architectures require an explicit validated wrapper. The adapter is optional;
the alignment worker has no Torch dependency.

An illustrative calculation is ready to run:

```bash
genocache-model examples/model-assumptions.json \
  --output local-data/model-report.json
```

The example assumes a 3.1 Gbp reference at stride 32, dimension 128 and FP32:
approximately 96.875 million vectors and 49.6 GB of raw vector payload for one
orientation. It also illustrates `alpha=0.10` and time-weighted fallback `0.75`,
which permit only `1/(0.10+0.75) = 1.18x` idealized speed. These are assumptions,
not measurements. [RESEARCH_MODEL.md](docs/RESEARCH_MODEL.md) gives the derivations,
limitations and current primary research sources.

## Evidence from this implementation

Verified on Linux on 7 October 2026: **214 public tests passed**, plus **173
subtests**, with Ruff clean. The final cloud acceptance run passed 217 tests,
including three additional private canonical-checkpoint diagnostic checks.
The public suite exercises native minimap2 v2.31, pysam,
CPU FAISS and CPU Torch. Twelve tests reproduce archived defects; a green legacy
audit means those old failures were reproduced. New-engine tests separately
assert corrected behavior. TorchScript emitted three deprecation warnings in
the installed Torch release.

The committed [integration receipt](reports/local-benchmark.json) is a single
small synthetic run with an 800 kb reference, 136 reads, substitutions, indels,
both orientations, three repeat copies and eight unrelated negative reads:

| Measurement | Result |
| --- | ---: |
| Correct primary start (within 100 bp) and strand among known-origin reads | 128 / 128 |
| Unrelated negative reads left unmapped | 8 / 8 |
| Secondary records preserved | 64 |
| BAM record multiset equal to independent native mapping + sort/index | Yes |
| Native equivalent-output wall time, already indexed | 0.228 s |
| GenoCache fresh batch, already indexed | 0.265 s |
| Identical verified cached batch | 0.023 s |

The fresh wrapper adds work. The exact repeat avoids remapping and sorting; it
does not measure the hit rate of a real workload or accelerate a novel read.
Tiny single-run timings are noisy and are not capacity estimates. The 128/128
result has a 95% Wilson lower bound of about 97.1%; it establishes no whole-genome,
variant-calling, clinical, or calibrated neural-mapping accuracy claim.

Reproduce the suite and receipt from `production/`:

```bash
MINIMAP2=/path/to/minimap2 pytest tests -q
ruff check .
python scripts/benchmark_local.py \
  --workdir /tmp/genocache-new-empty-benchmark \
  --minimap2 /path/to/minimap2 --report local-data/benchmark.json
```

Tests that need an absent optional backend skip explicitly. To exercise the
Torch adapter test, install the appropriate trusted PyTorch CPU/CUDA package in
the environment. The benchmark refuses a nonempty work directory so an old cache
cannot masquerade as a cold run.

## Deployment and operating scope

Use [DEPLOYMENT.md](docs/DEPLOYMENT.md) for Docker, Modal and AWS Batch commands.
The Docker image and authenticated Modal functions were exercised against real
S3 on 7 October 2026. A separate AWS EC2/SSM worker passed the same synthetic
mapping and artifact-download checks. Both providers preserved the independent
native record multiset, returned verified exact-repeat cache hits, placed
128/128 known-origin reads within 100 bp on the correct strand, and left all
eight unrelated controls unmapped. The AWS Batch definition has not been
deployed. These bounded synthetic calls do not qualify real WGS workloads.

The pilot used temporary AWS login credentials in a Modal secret. The exported
session expires independently of local CLI refresh. Refresh the secret and
redeploy before another bounded call; durable scoped authentication or a refresh
mechanism is required before unattended long jobs. See the deployment guide.

Local locks are appropriate for one machine's POSIX filesystem. Each worker
uses a local cache; S3 holds immutable packs and completed output manifests.
There is no distributed result lookup, cross-worker deduplication, eviction
scheduler, public HTTP service, or learned acceptance route in this release.
Cloud retries may repeat computation, but attempt-specific output prefixes
prevent mixed output. Consumers use the final `result.json` as the commit marker.

Use the [Codex runbook](docs/CODEX_RUNBOOK.md) to operate the tool and select the
next experiment from evidence. It includes concrete stop/go conditions so a
failed representation, ANN loss, region bug, and expensive fallback lead to
different actions.
