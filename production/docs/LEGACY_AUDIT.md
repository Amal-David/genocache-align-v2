# Audit of the archived GenoCache alignment pipelines

Audit date: 2026-10-07. Source baseline: repository commit `bcec02b` (November 2025 backup). Paths and line numbers below refer to that archived baseline, not to the new `production/` implementation.

## Finding

The archived code does not establish a production-quality neural aligner or an end-to-end speed advantage over minimap2. There are concrete score, strand, chaining, coordinate, output, and measurement errors. Training quality may also be a problem, but the available evidence does not support treating more training as the only remaining work.

The useful result from the neural experiments is a candidate-localization hypothesis. A candidate window within 10 kb of the origin is not a base-level alignment, and a plausible chain is not a calibrated mapping decision. The replacement production service must keep those contracts separate.

No original path was modified as part of this audit. The new service should not import these historical alignment wrappers.

## What the stored evidence actually establishes

| Evidence | Observation | What it measures | What it does not establish |
| --- | --- | --- | --- |
| `genocache-v4.1-production/validation/results/baseline_results.txt` | 100 rows; 26 labeled correct; 78 labeled mapped | Stored chromosome-only outcomes | Base coordinates, strand, CIGAR validity, downstream accuracy |
| `genocache-v4.1-production/validation/results/improved_results.txt` | 100 rows; 26 labeled correct; 74 labeled mapped | Same chromosome-only outcome after the documented change | An improvement in correctness |
| `genocache-v4.1-production/validation/results/giab_comparison.txt` | GenoCache 21/100 mapped, 108 s, 0.93 reads/s; minimap2 approximately 94/100 mapped | A historical report of mapped counts on real reads | Ground-truth accuracy; complete rerunnable timing comparison |
| `genocache-v4.1-production/development/training/nal_chr22_PRODUCTION_BACKUP/documentation/RESULTS.md:8-22` | Reports 76.5% versus 75.5% on 200 synthetic reads of about 2 kb | Coarse position recovery under a 10 kb tolerance | Whole-genome recall, exact alignment, statistically established superiority |
| The same report, lines 58-62 | Reports about 35 ms/read versus 260 ms/read | A documentation claim | A reproduced, identical-output end-to-end speedup |
| `GenoCache_Version6.1/COMPREHENSIVE_SUMMARY.md:287-290` | Claims about 10,000 reads/s and 5 GB inference memory | An unsupported capacity estimate in the summary | A measured production capacity |

The 76.5% versus 75.5% comparison amounts to two reads out of 200. A paired per-read comparison and uncertainty analysis are needed before interpreting it as superiority. The two methods are also evaluated at different levels of completeness: the archived Chr22 test performs candidate seeding and chaining, without final alignment.

The Chr22 `test_chr22_nal.py:79-106` accepts the best chain directly. Its success criterion at lines 124-134 allows a start-position error of 10,000 bp. It neither generates an exact CIGAR nor checks strand. The referenced per-read final log, model weights, indexes, and full reference are not present in this source-only backup, so this audit did not rerun that historical neural benchmark.

The Chr22 read generator (`generate_chr22_test_reads.py:56-65`) samples the chromosome without excluding N-rich regions, uses only forward-strand reads, and does not fix a random seed. Thus the low minimap2 figure cannot be interpreted as a general long-read accuracy ceiling; its read content, preset, reference and primary-alignment selection need inspection in a fresh benchmark.

## Blocking correctness defects

### 1. WFA costs are ranked in the wrong direction

`genocache-v4.1-production/genocache_core/fast_alignment.py:75-80` converts the WFA score into an absolute, positive cost. In that contract, zero is a perfect match and larger values are worse. `genocache_core/extend_phase.py:87-99` sorts that value in descending order and rejects values below 100.

An isolated execution of the actual archived definitions, with the documented WFA score contract, chooses a wrong candidate of cost 200 over a perfect candidate of cost zero. If the perfect candidate is the only candidate, it returns no alignment. This is a component-contract reproduction using a controlled backend double, not a claim about newly measured biological accuracy.

The adapter needs a typed distinction between `edit_cost` (smaller is better) and an alignment reward (larger is better). A backend switch must not silently change the ordering contract.

### 2. The requested alignment mode and actual genomic origin are discarded

`genocache_core/fast_alignment.py:21-33` stores `mode='semi-global'` without passing mode or ends-free settings to WFA. Lines 85-95 report `start_ref=0` and the full window's last coordinate, without extracting aligned endpoints. Lines 114-118 add 5 kb of flank and an additional read length to a candidate region. Lines 134-136 convert the hardcoded window origin into a genomic position.

These operations can report the search-window origin as the read origin and charge alignment against unnecessary flanks. Increasing retrieval recall or training time cannot fix those coordinates.

### 3. Neighbor multiplicity is mistaken for independent read support

`development/training/nal_chr22_PRODUCTION_BACKUP/alignment/chaining_nal.py:138-156,191-225` assigns one vote to every retrieved neighbor and sums these votes. Distinct query seeds are calculated only in the output metadata at line 225. There is no cap of one vote per query seed per candidate region.

The executable audit supplies 32 nearby reference windows retrieved by one query seed. The actual archived chainer emits score 32 with `num_seeds=1`. Such a score can dominate a region supported by several genuinely independent query positions. Dense overlapping indexes amplify this problem: the top K results may be many versions of one locus.

The older `genocache_core/adaptive_seeding.py:231-283` has a different failure: it groups all retained hits for a chromosome into one chain and rejects the entire group if any adjacent pair violates its tolerance. Multiple candidate loci on one chromosome can therefore erase a real local chain. It also discards all but ten hits per seed at line 234, regardless of the configured retrieval K.

### 4. Strand metadata is absent or ignored

The Chr22 index builder records only `strand='+'` (`indexing/build_chr22_index.py:103-110`). Its query path does not encode reverse complements (`alignment/seeding_chr22.py:47-62`). The chainer tries both diagonal formulas for every anchor (`alignment/chaining_nal.py:89-98,138-156`) without restricting anchors to their strand.

The executable audit consequently obtains both a positive and a negative chain from the same positive-only one-seed evidence. This is not a valid way to infer reverse-complement matches. The core V4 seeder also does not retain strand at all.

Any new implementation needs explicit orientation throughout indexing, retrieval, candidate clustering, local alignment, and SAM emission. A safe approach indexes the forward reference and queries both the read and its reverse complement, preserving oriented query coordinates.

### 5. Failed alignment is replaced by fabricated output

`development/training/nal_aligned/align_nal.py:318-365` returns a full-length `len(read)M` CIGAR with score zero on reference failure, missing WFA, empty CIGAR, and WFA exceptions. The caller at lines 232-243 still returns `mapped=True`, the unrefined chain coordinate, and a MAPQ derived from chaining. With a single chain of six votes, lines 265-269 produce MAPQ 60 even if WFA is absent.

The three-tier path has the same fallback (`test/align_adaptive_3tier.py:148-177`). A CIGAR `M` operation may contain mismatches; inventing it without alignment hides missing insertions, deletions and endpoints. Failure must instead trigger a valid fallback backend, an unmapped record, or a surfaced job failure, according to the actual failure type.

### 6. SAM fields and quality evidence are unreliable

`genocache_align.py:195-203` hardcodes MAPQ 60, writes flag zero, and writes an internal start position without adding one for SAM. It discards FASTQ qualities while loading reads at lines 147-149 and limits `@SQ` records to 50 reference sequences at lines 162-163.

The alternate writer in `genocache_core/sam_output.py` does add one to POS, but it manufactures `I` qualities at lines 261-262, ignores mismatches inside `M` operations when computing NM at lines 30-44, treats all `M` bases as matches for identity at lines 123-129, and never sets reverse strand in `format_sam_line` at lines 243-246. Its MAPQ is an uncalibrated heuristic. One correct helper does not repair the other code path.

## Retrieval, indexing and runtime causes

### Dense indexing has a large cost before search even starts

For reference length G, stride s, embedding dimension d and b bytes per component, a single-strand dense vector table occupies approximately `(G / s) * d * b` bytes.

For a 3 billion base reference, stride 32, dimension 128 and FP32, this is about 93.75 million vectors and 48 GB of vector data, excluding FAISS overhead, metadata, model, reference, and temporary copies. Indexing both strands doubles the vector table. FP16 or product quantization changes this arithmetic but does not eliminate preprocessing and recall costs.

The Chr22 builder accumulates all embeddings and then calls `vstack` (`indexing/build_chr22_index.py:100-142`). It stores a Python dictionary per locus and serializes the chromosome name as NumPy `U50` for every vector (lines 167-179): 200 bytes per locus just for that Unicode chromosome field, about 18.75 GB at the illustrative full-genome scale. The older production builder similarly collects embeddings, chromosome strings, and positions before concatenation (`development/indexing/build_production_index.py:89-145`) and creates more copies while training and adding the index (lines 181-188).

The replacement should stream shards, use integer contig IDs and coordinate arrays, persist a checksummed manifest, and measure peak memory rather than infer it from compressed index bytes alone.

### The successful Chr22 search path is not GPU-batched search

The Chr22 builder constructs a CPU `faiss.IndexFlatIP` (`indexing/build_chr22_index.py:149-164`). `alignment/seeding_chr22.py:28-30` loads it without transferring it to GPU. Lines 48-62 encode one seed at a time on GPU, transfer the result to CPU, and then perform one search. It repeats that for every read.

This path pays repeated Python, kernel-launch, synchronization and search overhead. The older core query loop similarly operates one seed at a time (`genocache_core/adaptive_seeding.py:100-109,166-172,306-309`) and recomputes a full rescue set at lines 335-343. The three-tier path caches separate full seeder objects by K (`test/align_adaptive_3tier.py:80-90`), even though K is a query parameter, creating avoidable model/index duplication when K changes.

Batched read encoding, batched query submission, persistent indexes, and incremental rescue are justified candidates for optimization. Their benefit must be measured separately from correctness and without assuming a GPU speedup in advance.

### Artifact identity is not validated

`genocache_align.py:60-71` loads index and metadata without checking reference identity, encoder weights, vector dimension, row count, normalization, index metric, or coordinate convention. `development/indexing/build_production_index.py:216-224` hardcodes embedding dimension 128 in metadata rather than deriving and validating it. The same builder constructs IVF-PQ without an explicit metric argument at lines 175-176 despite comments and callers assuming inner-product scores; metric semantics must be checked in the actual FAISS runtime before reuse.

The newer NAL index path explicitly sets inner-product metric (`nal_aligned/build_index_nal.py:201-204`), so this particular omission should not be attributed indiscriminately to every version. The backup lacks its production indexes; their actual dimensions and metrics were not inspected in this audit.

There is no verified alignment-result cache or bounded hot-region admission policy in the audited latest paths. `reference_cache` is an in-memory dictionary of whole chromosomes loaded by scanning FASTA. The project name should not be mistaken for evidence that a validated result-cache implementation already exists.

## Measurement errors that explain repeated false progress

1. **Mapping is often called accuracy.** `validation/scripts/04_run_improved.py:157-160` checks only chromosome. `nal_aligned/test/speed_accuracy_sweep.py:209-221` labels a read mapped whenever a chain exists, without checking a known origin.
2. **Stage time is apportioned rather than measured.** `speed_accuracy_sweep.py:191-194` assigns 70% of a combined operation to encoding and 30% to FAISS. Those percentages are assumptions, not profiler measurements.
3. **Refinement is absent from the three-tier timer.** `align_adaptive_3tier.py:207-220` stops timing before the actual WFA call at line 233.
4. **Rescue costs are undercounted.** The same script computes its weighted average from the tier at which the read exits (lines 273-283), omitting preceding tiers. The executable audit gives one read 100 ms in tier 1 and 200 ms in tier 2; the summary reports 200 ms even though the cumulative time is 300 ms.
5. **Backend and workload completeness differ.** A neural candidate localizer is compared with a mapper producing base-level alignments; the backup does not contain a comparable warm, preindexed minimap2 timing for the claimed 7.4x figure.
6. **The newest version is not an implemented full pipeline.** `GenoCache_Version6.1/README.md:24-38` depicts indexing and alignment programs, but only its training programs and documents exist in this backup. The training loader also drops entire reference sequences with at least 10% Ns or at most 1 Mb (`training/train_full_genome.py:42,55`), which is not equivalent to coverage of the entire reference.

The older root `FINAL_COMPREHENSIVE_SUMMARY.md` explicitly calls a four-case mock test definitive proof of 90-95% expected real-data accuracy and production readiness. A mock can establish that a supplied score selects a supplied candidate. It cannot establish candidate recall, real alignment score semantics, correct coordinates, safe output, throughput, or performance on biological repeats.

## New implementation contracts

The mathematical distinction to preserve is:

- `candidate_recall = P(true locus belongs to proposed candidate set)`;
- `conditional_refinement_accuracy = P(correct alignment | true locus was proposed)`;
- without an independent fallback, total correctness cannot exceed candidate recall.

A local exact alignment verifies sequence compatibility with a candidate. It does not prove uniqueness across the genome. Even a perfect local match can have a perfect competitor elsewhere. An embedding similarity or locally verified cache hit therefore cannot, by itself, justify a high MAPQ or skipping all other plausible loci.

For immediate service reliability, a proven global mapper should own final alignment output. Exact-result caches require reference, mapper version, full configuration, and input identity in the key. Where tie-breaking depends on read names, that identity must also be included or its equivalence explicitly tested. Similarity caches should supply candidate hints, not silently reuse another read's CIGAR. The experimental embedding lane should initially record candidate recall, downstream agreement and timing alongside authoritative output.

For a minimal new candidate chainer:

- Require explicit `(contig_id, orientation, query_seed_id, query_offset, reference_offset)` anchors.
- Transform query offsets into oriented coordinates before computing `diagonal = reference_offset - oriented_query_offset`.
- Group by contig and orientation, cluster compatible diagonal intervals, and cap support at one vote per query seed per candidate.
- Keep distinct genomic candidate windows; do not merge unrelated regions across a chromosome or let raw K inflate confidence.
- Reject candidate intervals outside the reference, validate finite similarities, and deduplicate seed coordinates before counting support.
- Treat diagonal clustering as candidate generation. A genuine base-level aligner supplies final starts, ends, CIGAR and edit cost.

## Reproducing the audit

Run from the repository root:

```bash
python -m unittest discover -s production/tests -p test_legacy_contracts.py -v
```

The same tests are pytest-compatible:

```bash
pytest production/tests/test_legacy_contracts.py -q
```

The suite has no third-party runtime requirement and executes actual archived class/function definitions in private namespaces. A controlled WFA score double isolates the backend contract; it does not replace a biological benchmark. No global `sys.modules` entries are changed. Tests assert observed historical failure behavior and intentionally do not use `xfail`. A green audit suite means those failures were successfully reproduced. Acceptance tests for the new service must assert the corrected behavior separately.
