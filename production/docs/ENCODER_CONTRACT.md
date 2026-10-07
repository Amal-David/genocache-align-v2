# Connecting a trained GenoCache encoder to retrieval

Updated 2026-10-07. This adapter converts an explicitly supplied, trusted
TorchScript encoder into the vector and coordinate files consumed by
genocache-retrieval. No model weights are distributed with this public package.
An existing checkpoint must pass the architecture and preprocessing contract
below before being exported; the existence of a weights file is not model
qualification.

The sequence, metadata, batching, validation, and artifact code is implemented.
On 2026-10-07, all 21 adapter tests passed with PyTorch 2.14.1+cpu, including real
TorchScript CPU loading and batched inference on a small pooling-module contract
fixture. That test also verified rejection of zero output vectors. Artifact
publication tests use an explicitly labeled runtime stub to isolate file and
metadata behavior. The real TorchScript test skips explicitly in installations
without Torch. These checks establish the adapter's data and runtime contract.
A separately recovered checkpoint also underwent a bounded synthetic diagnostic
outside this public package; that was direct model execution, not a validated
export through this adapter. Held-out genomic candidate recall, biological
alignment accuracy, and GPU throughput remain unmeasured.

## Encoder contract

| Property | Required behavior |
|---|---|
| Artifact | A trusted TorchScript module saved with torch.jit.save or ScriptModule.save |
| Input | One float32 tensor with shape [B, 4, L] |
| Channel order | A, C, G, T |
| Length | Exactly the CLI window size L; default 512 |
| Batch dimension | Support B from 1 through the selected batch size, including a final partial batch |
| Ambiguous IUPAC DNA | All four channels zero by default, or reject the input with --ambiguous error |
| Output | One floating tensor with shape [B, D], with D positive and constant across batches |
| Numeric constraints | Output must remain finite after float32 conversion; float32 L2 norms must be finite and greater than 1e-12 |
| Evaluation behavior | Deterministic evaluation without stochastic augmentation, batch-dependent pooling across examples, or hidden persistent state |
| Coordinates | Every row refers to a zero-based, half-open interval on the original forward sequence |

The adapter calls model.eval() and executes under torch.inference_mode(). These
are separate controls: evaluation mode changes modules such as dropout, while
inference mode disables autograd bookkeeping. Its inference call, input transfer,
and output transfer operate on whole batches. It never calls the model once per
seed or transfers each seed's output separately.

PyTorch's documented TorchScript loader accepts serialized ScriptModules and
supports explicit device remapping. A TorchScript artifact is executable model
code; only use an artifact whose origin you trust. See the primary
[torch.jit.load documentation](https://docs.pytorch.org/docs/2.9/generated/torch.jit.load.html)
and [PyTorch security policy](https://github.com/pytorch/pytorch/security/policy).
The [inference-mode documentation](https://docs.pytorch.org/docs/2.9/generated/torch.autograd.grad_mode.inference_mode.html)
describes the runtime context used here.

The encoder's training preprocessing must match this table. A model trained with
another channel order, k-mer tokens, ambiguity probabilities, special tokens,
different pooling, or different sequence lengths requires its own explicit
wrapper and validation. Changing its input representation without checking the
training contract can yield plausible vectors with poor genomic retrieval.

## Sequence and window policies

Reference input must be FASTA. Read input may be FASTA or multiline FASTQ.
Both may be gzip compressed; detection uses file bytes rather than the extension.
Headers use their first token as the identifier, with SAM-compatible validation.
Empty sequences, non-IUPAC DNA, broken quality lengths, duplicate identifiers,
and malformed inputs fail before encoding. FASTQ quality scores are validated
but are not input to this encoder contract.

The reference mode emits only full windows, in contig order followed by ascending
start coordinate. It resets the stride grid at every contig boundary and emits
only the forward orientation. For a contig of length G, window length L, and
stride s, its row count is

\[
N_G=\max\left(0,\left\lfloor\frac{G-L}{s}\right\rfloor+1\right).
\]

The starts are 0, s, 2s, and so on while start + L is at most G. There are no
cross-contig windows, no padding, and no extra off-grid terminal window.
Contigs shorter than L contribute zero rows and are counted as skipped.
If stride exceeds L, internal gaps also exist. Use the actual census in the
manifest when sizing this index: the general mathematical model's optional
right-clipped tiling estimate is a different window policy.

The reads mode selects up to eight distinct full seed intervals per read by
default. For read length R at least L, let

\[
K=\min(\text{seeds},R-L+1).
\]

If K is one, it uses the midpoint start floor((R-L)/2). Otherwise the starts are

\[
q_i=\left\lfloor\frac{i(R-L)}{K-1}\right\rfloor,\quad 0\leq i<K.
\]

This includes the two ends and removes duplicate windows when a read is only
slightly longer than L. Each seed is encoded twice, first as its forward sequence
and then as its reverse complement. Thus the usual maximum is 16 query vectors
per read. The two orientations share seed_id so they cannot be counted as two
independent intervals.

Both orientations retain the same original forward-read query_start and
query_end. For example, an interval [400, 912) in a 2,000-base read keeps those
coordinates when its reverse complement is encoded. Its query_strand becomes
minus; the coordinates do not become [1088, 1600).

Reads shorter than L are skipped and counted by default. They have no query rows
and must remain eligible for the native alignment path. Set --short-reads error
to require every read to have an encoder seed. There is no silent padding or
length substitution. The default maximum accepted read length is 2,000,000 bases;
increase --max-read-length deliberately if the data requires it. The parser
checks the limit before accumulating an oversized read.

Uniform seeds are a baseline sampling policy, not a guarantee of covering every
informative region. They can miss a small useful interval in a long or repetitive
read. Likewise, an overlapping reference window does not by itself imply that an
encoder will retrieve it. Truth-based candidate recall, stratified by read length,
error profile, repeats, and difficult regions, must decide whether this policy is
adequate.

## Files, memory, and provenance

The output directory is created only after a successful run. It contains:

| File | Contents |
|---|---|
| vectors.npy | N by D little-endian float32 dense array, in metadata row order |
| windows.jsonl | Reference mode: chrom, start, end, strand |
| queries.jsonl | Reads mode: read_id, read_length, read_ordinal, query_start, query_end, seed_id, query_strand |
| manifest.json | Input and encoder hashes, configuration, runtime, row counts, dimensions, file hashes, identities, and timings |

Exactly one of the two metadata files is present. Vectors are stored without
normalization; the retrieval importer normalizes them for cosine-style inner
product search. Query metadata fields match VectorIndex.anchors directly.

Input and encoder bytes are copied into verified, temporary snapshots before
use. The manifest's reference_sha256 or reads_sha256 hashes the exact input bytes,
including compression and headers. encoder_sha256 hashes the exact artifact
loaded. Temporary snapshots are removed after the output checksums are computed.
An uncompressed FASTA and a gzip encoding of the same bases intentionally have
different input identities.

The first parsing pass counts rows and checks unique names using a disk-backed
SQLite table with a bounded page cache. The second pass produces windows and
metadata as a stream. A disk-backed NumPy array is allocated after the first
batch reveals D. No list of all embeddings, all read sequences, or all contig
strings per vector is accumulated. Reference buffering is bounded by one input
line chunk plus one window. Read buffering is bounded by the explicit maximum
read length. Sequence batches, input tensors, output batches, and model weights
are additional live memory. Operating-system file cache may also grow as the
memmap is written.

Provision local disk for the input snapshot, model snapshot, identifier database,
N x D x 4 vector payload, JSONL metadata, and later retrieval-pack files. JSONL
metadata can be substantial at genome scale. This adapter streams it; the
retrieval importer converts reference positions into numeric contig IDs and a
separate contig dictionary. Building the reference index is an amortized offline
cost; encoding reads remains part of per-request cost.

The manifest records two identities:

- recipe_id hashes the input/model digests, algorithm configuration, and runtime
  settings, excluding paths and elapsed times.
- output_id additionally hashes the vector and metadata files. Runs with the
  same recipe and identical output bytes have the same output_id even though
  their timestamps or elapsed times differ.

The adapter requests deterministic PyTorch algorithms, sets the random seed
(default 17), disables TF32, and records device, Torch, CUDA, NumPy, and thread
settings. This does not promise bitwise equality across hardware or releases, or
repair a stateful supplied model. PyTorch documents those reproducibility limits
in its [reproducibility guidance](https://docs.pytorch.org/docs/2.14/notes/randomness.html).
Use the output file hashes to detect actual differences.

The timing object reports:

- encoding_seconds: summed batch input packing, device transfer, model execution,
  output transfer, and numeric validation. CUDA output transfer is blocking once
  per batch, so completed model work is included.
- model_load_seconds: runtime import/setup and loading the trusted artifact.
- full_elapsed_seconds: work from function entry through source snapshotting,
  both input passes, model setup, encoding, vector/metadata writing and checksums,
  and temporary input cleanup. The final small manifest write and directory
  publication occur immediately afterward.

Measure caller wall time as well for end-to-end comparisons, then add index
loading/search, candidate processing, native fallback, and output writing.
Reporting only encoding_seconds would omit real request costs.

## Run the adapter and reuse its artifacts

Commands below assume the repository root and an installed production package.
Install PyTorch separately for the selected CPU or CUDA environment; importing
the production package or asking the script for help does not require Torch.
The repository's optional torch dependency group is also available. No command
below downloads weights, starts training, or chooses a random replacement model.

Encode a reference once:

~~~bash
python production/scripts/encode_torchscript.py reference \
  --input /data/reference.fa \
  --encoder /models/trained_encoder.torchscript.pt \
  --output /data/reference-vectors \
  --window 512 --stride 256 --batch-size 256 --device cuda:0
~~~

Build the exact vector index first so ANN approximation can later be measured
separately from encoder recall:

~~~bash
python - <<'PY'
import json
import subprocess
from pathlib import Path

source = Path("/data/reference-vectors")
manifest = json.loads((source / "manifest.json").read_text())
subprocess.run([
    "genocache-retrieval", "build",
    "--vectors", str(source / "vectors.npy"),
    "--windows", str(source / "windows.jsonl"),
    "--output", "/data/reference-pack",
    "--reference-sha256", manifest["reference_sha256"],
    "--encoder-sha256", manifest["encoder_sha256"],
    "--kind", "flat",
], check=True)
PY
~~~

Encode reads with the same exported model and preprocessing:

~~~bash
python production/scripts/encode_torchscript.py reads \
  --input /data/reads.fastq.gz \
  --encoder /models/trained_encoder.torchscript.pt \
  --output /data/query-vectors \
  --window 512 --seeds 8 --max-read-length 2000000 \
  --batch-size 256 --device cuda:0
~~~

Retrieve candidate windows while carrying every seed's orientation and interval:

~~~bash
python - <<'PY'
import json
import subprocess
from pathlib import Path

source = Path("/data/query-vectors")
manifest = json.loads((source / "manifest.json").read_text())
subprocess.run([
    "genocache-retrieval", "search",
    "--pack", "/data/reference-pack",
    "--queries", str(source / "vectors.npy"),
    "--metadata", str(source / "queries.jsonl"),
    "--encoder-sha256", manifest["encoder_sha256"],
    "--output", "/data/candidate-hits.jsonl",
    "--k", "32", "--batch-size", "256", "--device", "cpu",
], check=True)
PY
~~~

The search example uses CPU because the default retrieval extra is FAISS CPU.
Request CUDA search only in an environment with compatible GPU FAISS. Encoder
CUDA and FAISS CUDA are separate capabilities. Both paths fail explicitly if a
requested GPU capability is unavailable.

The reference and query encoder digests must match. Also compare their window,
channel-order, and ambiguity policies before a run; retrieval's current import
interface checks the artifact digest and vector dimension but does not consume
this adapter's whole preprocessing manifest.

Use a new output directory for changed settings. Preserve and reuse the encoded
reference and the retrieval pack for repeated requests. IVF-PQ can replace the
flat index after an exact-vector audit has measured the additional ANN recall
loss. Candidate hits remain evidence for evaluation and routing; they are not a
SAM alignment or calibrated mapping quality.

## Exporting an existing state_dict checkpoint

A state_dict is parameter data, not a portable forward function. Its original
architecture, constructor settings, channel order, feature pooling, normalization,
and checkpoint key layout are required. Restore these from the training code and
configuration. A mismatched architecture is a hard error; do not use partial
parameter loading to make an unrelated network appear compatible.

The following is an export outline for the owner of the actual trained model.
The imported module, constructor configuration, and checkpoint layout must come
from that model's training repository:

~~~python
import torch
from your_training_package import GenomeEncoder  # actual original architecture

model = GenomeEncoder(**original_training_configuration)
state = torch.load("trained_state_dict.pt", map_location="cpu", weights_only=True)
# If the original file is a checkpoint dictionary, select its documented
# state_dict key here. Do not guess the architecture or discard mismatched keys.
model.load_state_dict(state, strict=True)
model.eval()

# Add an explicit wrapper first if the trained model's public forward contract
# differs from [B, 4, L] -> [B, D]. Preserve the trained computation exactly.
exported = torch.jit.script(model)
exported.save("trained_encoder.torchscript.pt")
~~~

Some architectures need an export-specific wrapper or a validated tracing path.
Check several real windows, both orientations, ambiguity behavior, and batch
sizes one, a full batch, and a final partial batch against the original model.
Check output dimension and numerical agreement before using the export to build
a genome index. TorchScript is used here as an explicit artifact bridge; it does
not remove the need to recover the actual architecture and trained checkpoint.
The tested PyTorch 2.14.1 runtime emits deprecation warnings for TorchScript
tracing and loading while executing this contract successfully. A future
torch.export adapter should be treated as a separate, versioned artifact
contract; this adapter does not silently reinterpret an existing checkpoint.

Run the local contract tests with:

~~~bash
PYTHONPATH=production python -m unittest discover \
  -s production/tests -p test_encoder_adapter.py -v
~~~

The optional real TorchScript test uses a small pooling module only as an input/
output fixture. It is never a replacement genomic encoder. A successful real
trained-model run, exact-vector candidate recall measurement, full end-to-end timing,
and held-out genomic accuracy evaluation are still required before the learned
retrieval path can be allowed to influence production alignments.
