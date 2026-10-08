# Codex runbook: operate, diagnose, and decide

Updated 2026-10-07. All commands below run **from `production/`**, with Python
3.12+, the package installed, and a compatible minimap2 binary on `PATH`.
Use [DEPLOYMENT.md](DEPLOYMENT.md) for installation and cloud operation.

## 1. Keep the current production contract

`genocache.engine` runs native minimap2 on every novel input batch. An exact
cached batch reuses a verified artifact. There is currently **no neural accepted
alignment path**, selective neural fallback, per-read result cache, or automatic
hot-region cache. The native reference index resides on disk and loads for each
cold batch. Embedding commands run separately and produce shadow evidence.

Preserve native SAM records, including CIGAR, MAPQ, sequence/quality, flags,
secondary/supplementary records, and tags. BAM output adds sorting and CSI
indexing. A locally perfect sequence match or a high cosine score does not prove
the correct genomic origin; never manufacture SAM or MAPQ from proposals.
See [engine.py](../genocache/engine.py) and [LEGACY_AUDIT.md](LEGACY_AUDIT.md).

## 2. Register one experiment, then establish the native result

Choose one causal hypothesis, one changed mechanism, a target speedup `S`, an
error budget, minimum recall by stratum, and a compute/time budget **before**
examining results. Save these choices and the exact executed commands in the
new run directory. Keep input/truth/model bytes fixed during a comparison.

```bash
set -euo pipefail
export GC_ROOT=/data/genocache-worker
export GC_REFERENCE=/data/reference.fa
export GC_READS=/data/reads.fastq.gz
export GC_TRUTH=/data/truth.jsonl
export GC_RUN="$(mktemp -d /data/genocache-experiment-XXXXXXXX)"
python -m pip freeze > "$GC_RUN/environment.txt"
git rev-parse HEAD > "$GC_RUN/revision.txt"
tar --exclude='__pycache__' -czf "$GC_RUN/source.tgz" \
  genocache scripts pyproject.toml requirements-runtime.txt
cp "$GC_TRUTH" "$GC_RUN/truth.jsonl"
python -m genocache.cli --root "$GC_ROOT" doctor > "$GC_RUN/doctor.json"
python -m genocache.cli --root "$GC_ROOT" index \
  --reference "$GC_REFERENCE" --preset map-ont --threads 4 > "$GC_RUN/index.json"
GC_REF_ID=$(python -c 'import json,sys; print(json.load(open(sys.argv[1]))["reference_id"])' "$GC_RUN/index.json")
python -m genocache.cli --root "$GC_ROOT" align \
  --reference-id "$GC_REF_ID" --reads "$GC_READS" --threads 4 \
  --output-format bam > "$GC_RUN/alignment.json"
GC_BAM=$(python -c 'import json,sys; print(json.load(open(sys.argv[1]))["output_path"])' "$GC_RUN/alignment.json")
python -m genocache.evaluate alignments --alignment "$GC_BAM" \
  --truth "$GC_RUN/truth.jsonl" --tolerance 100 --output "$GC_RUN/mapping.json"
```

Choose the supported preset appropriate to the reads. Declare position tolerance;
this evaluator measures primary start/strand agreement, not CIGAR or variant
accuracy. Inspect `cache_hit` and `metrics`: separate fresh jobs, repeated exact
batches, index construction, conversion, validation, and transfer costs.

## 3. Truth and strata must be independent

Truth JSONL has unique read IDs and zero-based forward-reference coordinates:

```json
{"read_id":"r1","read_length":2048,"reference_span":2075,"stratum":"tandem_repeat","origins":[{"chrom":"chr1","start":100000,"strand":"-"}]}
{"read_id":"negative1","read_length":1800,"stratum":"known_unmappable","origins":[]}
```

Retrieval requires positive `read_length` for every read and positive
`reference_span` for every read with origins. Omit the span on deliberately
unmappable reads; zero is invalid. Multiple valid origins are allowed. Anchor
read lengths must match truth. Missing and unmapped reads stay in denominators.
Do not label a read unmappable merely because the baseline failed to map it.

Use independently generated origins or an independently validated benchmark.
Minimap2 agreement is a useful regression check, not independent truth. Hold out
whole samples or genomic intervals so overlapping training windows cannot leak
into evaluation. Define strata from independent truth: unique regions, repeats,
homopolymers, segmental duplications, read length/error classes, and difficult
intervals. `stratum` is one label per read; use declared combined labels or
separate immutable truth views for multiple stratification schemes.

## 4. Run embeddings as a bounded shadow experiment

Supply a trusted, trained TorchScript artifact under the explicit
[ENCODER_CONTRACT.md](ENCODER_CONTRACT.md). The adapter accepts `[B,4,L]` A/C/G/T
one-hot input and `[B,D]` output. Preprocessing and length must match training.
Do not substitute random weights or a contract-test model. Start with a small
reference panel: FlatIP duplicates dense vectors in its resident index, while
the pack also retains dense vectors on disk for exact audits.

```bash
export GC_ENCODER=/models/trained_encoder.torchscript.pt
python scripts/encode_torchscript.py reference --input "$GC_REFERENCE" \
  --encoder "$GC_ENCODER" --output "$GC_RUN/reference-vectors" \
  --window 512 --stride 256 --batch-size 256 --device cpu
python scripts/encode_torchscript.py reads --input "$GC_READS" \
  --encoder "$GC_ENCODER" --output "$GC_RUN/query-vectors" \
  --window 512 --seeds 8 --batch-size 256 --device cpu
GC_REF_SHA=$(python -c 'import json,sys; print(json.load(open(sys.argv[1]))["reference_sha256"])' "$GC_RUN/reference-vectors/manifest.json")
GC_ENC_SHA=$(python -c 'import json,sys; print(json.load(open(sys.argv[1]))["encoder_sha256"])' "$GC_RUN/reference-vectors/manifest.json")
GC_QUERY_ENC_SHA=$(python -c 'import json,sys; print(json.load(open(sys.argv[1]))["encoder_sha256"])' "$GC_RUN/query-vectors/manifest.json")
python -m genocache.retrieval build --vectors "$GC_RUN/reference-vectors/vectors.npy" \
  --windows "$GC_RUN/reference-vectors/windows.jsonl" --output "$GC_RUN/vector-pack" \
  --reference-sha256 "$GC_REF_SHA" --encoder-sha256 "$GC_ENC_SHA" \
  --kind flat > "$GC_RUN/vector-index.json"
python -m genocache.retrieval search --pack "$GC_RUN/vector-pack" \
  --queries "$GC_RUN/query-vectors/vectors.npy" --metadata "$GC_RUN/query-vectors/queries.jsonl" \
  --encoder-sha256 "$GC_QUERY_ENC_SHA" --k 32 --batch-size 256 --audit-exact 32 \
  --output "$GC_RUN/anchors.jsonl" > "$GC_RUN/retrieval-timing.json"
python - <<'PY'
import json, os
from pathlib import Path
p = Path(os.environ["GC_RUN"])
pack = json.loads((p / "index.json").read_text())
contigs = {c["name"]: c["length"] for c in pack["input_stats"]["contigs"]}
(p / "contigs.json").write_text(json.dumps(contigs))
PY
python -m genocache.evaluate retrieval --anchors "$GC_RUN/anchors.jsonl" \
  --truth "$GC_RUN/truth.jsonl" --contigs "$GC_RUN/contigs.json" \
  --positional-slop 256 --max-candidates 8 --min-seeds 2 \
  --output "$GC_RUN/candidate-recall.json"
```

These are declared example sampling settings, not recall guarantees. The adapter
encodes both query orientations, retains original forward-read offsets, and
counts short reads that it skips. Keep those reads in truth and native mapping.
Shifted overlapping windows can fail retrieval even when a true locus exists.

`exact_vector_audit.top_k_overlap` compares ANN neighbor IDs with dense cosine
neighbors; it is **not biological locus recall**. Flat search is the small-panel
oracle. `candidate-recall.json` instead requires a same-contig, same-strand
proposal containing the **entire** independent reference span. Keep both metrics.
Move to IVF-PQ only after exact-vector biological recall is adequate. Its builder
requires `D % pq_m == 0` and at least `max(39*nlist, 39*256)` training rows within
the 262,144-row training cap. Tune one of K, nprobe, or compression at a time.

## 5. Diagnose the first stage that loses the true locus

For every miss, join original query intervals, reference window rows, exact top-K
IDs, ANN IDs, and proposal diagnostics by read ID. Reuse `VectorIndex.exact_search`
on a bounded sample; it is not a genome-scale serving strategy.

| Earliest failure | Evidence to retain | Next experiment |
|---|---|---|
| No usable seed/window representation | Short-read counts, window grid and true interval | Change seed/window coverage or the trained representation; preserve the native route. |
| Exact cosine top-K lacks the true locus | True-window rank versus hard negative loci | Test shift/RC/error invariance and hard-negative training; FAISS tuning cannot repair exact-ranking failure. |
| Exact retrieves the locus; ANN drops it | Exact/ANN ID overlap and biological recall at fixed K | Change nprobe or compression within the latency/memory budget. |
| Hits survive; region proposals lose it | Strand, interval uncertainty, distinct intervals, coverage and truncation reasons | Fix coordinates or proposal budget/slop; do not treat every window start as a point anchor. |
| Correct region exists; local verification fails | Full alignment coordinates, clipping, edit operations and competing loci | Investigate a separate verifier; current shadow tools do not implement an accepted alignment route. |
| Accuracy holds; speed fails | End-to-end phase times and fallback subset cost | Optimize the measured dominant cost or close the branch. |

Locate hot intervals using independently known origins and aggregate missed
reads, baseline time, candidate count, and truncation by stratum/interval. Keep
request frequency separate from difficult-region cost. Admission to a future
hot-region cache needs measured reuse and miss cost; frequency alone is not a
correctness signal. These reports/admission decisions are not automated today.

## 6. Apply the mathematical stop/go gates

For the same novel-read workload and timing boundary, let `T0` be native time,
`alpha = T_front/T0`, and `f_time = T_native(fallback subset)/T0`. Include encoding,
search, proposals, verification, transfers and merge work in `T_front`.

\[
S_{\mathrm{possible}}=1/(\alpha+f_{\mathrm{time}}),\qquad
\alpha+f_{\mathrm{time}}\le 1/S_{\mathrm{target}}.
\]

Use measured fallback **time**, not fallback read fraction. If alpha exceeds
`1/S_target`, that design fails even with zero fallback. Today all novel reads
still take native mapping, so `f_time=1`; adding shadow retrieval cannot establish
an end-to-end speedup. Overlapping pipelines require critical-path measurements.

For a future accepted fraction `a`, test `a*e_accept + (1-a)*e_fallback <= epsilon`
with independently measured conditional errors. Use per-stratum confidence
intervals, not zero observed failures as proof of zero risk. Candidate recall
only bounds the possibility of success; it does not calibrate acceptance/MAPQ.
Failing any predeclared accuracy gate blocks promotion, even if average speed wins.

Save measured inputs in `model-input.json` using [RESEARCH_MODEL.md](RESEARCH_MODEL.md)
and run `python -m genocache.model_cli "$GC_RUN/model-input.json" --output "$GC_RUN/model-report.json"`.
The example assumptions file contains illustrations, not measurements. Keep
unknown quantities unknown. If a gate fails, record the failure and stop that
configuration. Open another experiment only with a different causal hypothesis;
do not rebrand the same result or widen correctness tolerance after seeing it.

## 7. Seal the evidence and state its limits

Save `decision.json`: hypothesis, changed setting, acceptance thresholds, timing
boundary, hardware, seeds, data splits, per-stratum results, gate outcomes, and
the next permitted action. Retain input/model versions and native job/reference
manifests at their recorded hashes; worker paths alone are not durable records.
Seal each completed directory and start a new directory for every revision:

```bash
python - <<'PY'
import json, os
from pathlib import Path
from genocache.retrieval import file_sha256
p = Path(os.environ["GC_RUN"])
assert (p / "decision.json").is_file(), "Write the experiment decision first"
assert not (p / "checksums.json").exists(), "A sealed run must not be rewritten"
files = {str(f.relative_to(p)): file_sha256(f) for f in sorted(p.rglob("*")) if f.is_file()}
(p / "checksums.json").write_text(json.dumps(files, indent=2) + "\n")
PY
chmod -R a-w "$GC_RUN"
```

Run `python -m pytest tests -q` and `ruff check genocache scripts tests` after code
changes. The original build-time suite recorded 212 tests. The repaired Linux
cloud acceptance run passed 214 public tests plus three private checkpoint
diagnostic checks, with 173 subtests; this is regression evidence,
not exhaustive validation. Inspect skips and native dependency availability.
The [local receipt](../reports/local-benchmark.json) covers an 800 kb synthetic
reference: 128/128 known origins were correct and 8 known negatives stayed
unmapped; native alignment record multisets matched. Fresh GenoCache work was
slower than the equivalent native mapping/conversion pipeline. Its warm timing
replays the identical batch and says nothing about novel-read throughput or
real-workload cache hit rate. Human-genome accuracy, learned-model gains, GPU
throughput, and cloud end-to-end performance remain separate experiments.
