# GenoCache: an auditable model for alignment acceleration

Research checked **7 October 2026 (UTC)**. Calculations below are analytical estimates, not measured GenoCache speedups. `genocache.model.analyze_model` does not run alignment, inspect a GPU, or certify a model. Missing inputs and unknown results are represented by JSON `null`; explicit zero means the caller supplied zero.

## Decision

Use learned retrieval to propose candidate intervals, and promote it into an accepting fast path only after measuring candidate recall, mapping errors, and the full cost of fallback. Keep a pinned conventional engine available. Accelerating vector lookup cannot repair a representation that retrieves the wrong genomic locus; verifying one candidate cannot establish that another repeat copy was not missed.

The useful change is to optimize **verified mapping per unit time/cost subject to an error budget**, with each stage separately observable. A faster embedding encoder, lower training loss, high vector-neighbor recall, and high mapping accuracy are different outcomes.

## What the primary sources establish

| Source and date | Reported evidence | Consequence for this project |
| --- | --- | --- |
| [Embed-Search-Align (ESA)](https://arxiv.org/html/2309.11087v6), published in Bioinformatics, March 2025 | Learned embeddings can retrieve human-genome fragments for 250-base reads; the paper reports roughly 99% alignment accuracy under its evaluation definitions. Its implementation processes about 10,000 reads/minute versus roughly 1 million for Bowtie. Its real PacBio reads are cropped to 250 bases and candidates undergo Smith–Waterman alignment. | Evidence that learned retrieval is possible, not evidence of a production long-read speed advantage. Read length, per-chromosome top-K, success definition, and alignment work must accompany any comparison. |
| [NeurALigner](https://openreview.net/forum?id=aZjdkgvPbh), ICLR 2026 submission, initially 2025 | The indexed primary abstract describes GPU embedding-based seed retrieval for long-read alignment. Direct access encountered OpenReview's verification wall during this review. | Relevant architecture, but quantitative claims, final acceptance status, code availability, and reproducibility were not verified here. Do not insert an assumed speedup from its abstract. |
| [LOCALE paper](https://pubmed.ncbi.nlm.nih.gov/42182481/) and [authors' code](https://github.com/ryansynk/locale), May 2026 | Trains local-alignment embeddings with overlapping crops and substitution/indel augmentation. The reported task retrieves SRA accessions, with 62.4% recall under the paper's 10% mutation setting; it is not a long-read-to-SAM benchmark. Code and a checkpoint are available. | Useful training and diagnostic reference. Local containment/ranking is a better objective than generic biological semantics or global edit-distance regression. |
| [minimap2 release notes](https://github.com/lh3/minimap2/blob/master/NEWS.md) and [manual](https://github.com/lh3/minimap2/blob/master/minimap2.1), v2.31, 19 May 2026 | v2.31 fixes alignment-flag and inversion-alignment bugs. `lr:hq` targets accurate long reads; `map-ont` targets noisier long reads. | Pin the version and select the preset by chemistry/error profile. A comparison against an old or unsuitable preset is not a valid present-day win. |
| [mm2-fast](https://github.com/bwa-mem2/mm2-fast), paper 2022; repository checked 2026 | Authors report up to 1.8× CPU acceleration. The repository identifies v2.24 as its comparison baseline and requires effectively unlimited chain skip for the optimized chaining equivalence check. | A useful CPU comparator, with its exact version and parameters recorded. |
| [mm2-gb](https://github.com/Minimap2onGPU/mm2-gb), paper 2024 | Authors report 2.57–5.33× gains on 10–100 kb ONT reads and 1.87× on 100–300 kb reads with MI210 versus their 32-core mm2-fast baseline. | Hardware and read-length-specific reported results, not a promised GenoCache or NVIDIA speedup. |
| [NVIDIA Parabricks minimap2 documentation](https://docs.nvidia.com/clara/parabricks/tool-reference/tools/minimap2), documentation v4.7.1, checked 7 October 2026 | Describes accelerated KSW2 and sorted BAM/CRAM output; its CPU equivalence recipe uses minimap2 v2.31 with a custom `map-pbmm2` preset and a KSW2 loop-bound fix. | Compare equivalent versions, presets, fixes, and output work. Mapping-only PAF time cannot be compared directly with sorted BAM production. |
| [WFA-GPU](https://github.com/quim0/WFA-GPU), [paper](https://academic.oup.com/bioinformatics/article/39/12/btad701/7425447), 2023 | CUDA gap-affine global pairwise alignment, with exact and heuristic modes. | An extension backend candidate, not a whole-genome candidate finder. Its scoring and boundary conditions must match the intended mapping semantics. |
| [Winnowmap](https://github.com/marbl/Winnowmap), Winnowmap2 paper 2022 | Weighted minimizers and confidently alignable substrings address repeat mapping. The authors explain that the highest pairwise score may select the wrong repeat copy under a non-reference allele. | Evaluate repeat strata and preserve ambiguity/split evidence; “best score among retrieved candidates” is not proof of correct origin. |
| [Faiss index documentation](https://github.com/facebookresearch/faiss/wiki/Faiss-indexes) and [cuVS CAGRA](https://docs.nvidia.com/cuvs/user-guide/api-guides/indexing-guide/cagra), checked 2026 | Faiss documents vector, ID, graph and PQ storage; cuVS supplies GPU ANN. | Existing ANN implementations are appropriate acceleration components. Their nearest-vector recall must be evaluated separately from true-locus recall. |

These sources support a hybrid investigation. None substitutes for a benchmark of this repository on the user's reads and hardware.

## 1. Candidate recall is the first accuracy ceiling

Let `C` mean the candidate set includes at least one acceptable true locus, `H` mean chaining preserves it, `V` mean verification succeeds, and `S` mean selection/reporting is correct. Then, for a pipeline whose successful path requires these events:

\[
P(\mathrm{correct})=P(C)P(H\mid C)P(V\mid C,H)P(S\mid C,H,V).
\]

This is conditional probability, not an independence assumption. The API's `*_given_*` fields must use conditioning on **all preceding successful stages**.

Measure at least these distinct quantities:

1. ANN recall against exact nearest-neighbor search on the same vectors.
2. Biological candidate recall against known origin or an explicitly defined acceptable-locus set.
3. Chain survival conditional on a correct candidate being present.
4. Final mapping errors, MAPQ calibration, aligned bases, primary/secondary/supplementary behavior, and downstream variant results.

For two candidate generators with recalls `r_a` and `r_b`, without assumptions about their dependence:

\[
\max(r_a,r_b)\le r_{a\cup b}\le\min(1,r_a+r_b).
\]

Two 99% systems can both miss exactly the same 1%. Obtain union recall from paired reads, not the independence shortcut `1 - (1-r_a)(1-r_b)`.

## 2. Index scale, window boundaries and positional uncertainty

The module counts one right-clipped reference window at starts `0, s, 2s, ... < contig_length`. For lengths `G_c`, stride `s`, and explicitly supplied orientation count `o`:

\[
N=o\sum_c\left\lceil G_c/s\right\rceil.
\]

If only `genome_bases` is supplied, the result treats it as one contig and labels that approximation. For uniform full-window length `F` and maximum reference span `q*` of a query chunk, every valid integer starting position is covered when `s <= F - q* + 1`. Include indel expansion in `q*`; do not confuse tiling coverage with retrieval recall.

**A retrieved window is not an exact anchor.** A query chunk starting at query coordinate `q` and contained in a reference window `[p,p+F)` produces a possible diagonal interval approximately `[p-q, p+F-q*-q]`. Using `p-q` as an exact anchor can reject valid chains. Cluster intervals, then refine local coordinates, or propagate uncertainty into chaining. Retain multiple chains when structural variation changes the diagonal.

For dimension `d` and `b` bytes/component, vector payload is `N*d*b`. Example calculations for a 3.1 Gbp reference, 256 dimensions and fp16:

| Stride | One-orientation rows | One-orientation vector payload | Two-orientation payload |
| --- | ---: | ---: | ---: |
| 1,000 bp | 3.1 million | 1.59 GB | 3.17 GB |
| 500 bp | 6.2 million | 3.17 GB | 6.35 GB |
| 100 bp | 31 million | 15.87 GB | 31.74 GB |

GB is decimal. Add IDs, metadata, graph links, encoder weights, scratch and allocator overhead. The calculator deliberately reports **subtotals**, not guaranteed peak VRAM. Its HNSW estimate includes only the base layer's `2*M` links per row. Uniform PQ with `m` subquantizers and `b_pq` bits requires `ceil(m*b_pq/8)` bytes per row and `d*2**b_pq` codebook scalars, plus IDs and other bookkeeping. Retaining original vectors for refinement can remove much of the apparent memory saving.

Dense exact vector search costs approximately `2*N*d` floating-point operations per query vector, excluding selection and traffic. A read with `u` query chunks multiplies that by `u`. IVF estimates use `N*nprobe/nlist` scanned rows under balanced-list assumptions, plus centroid comparisons. They are not worst-case or recall guarantees; PQ lookup arithmetic is different.

## 3. Stage acceleration and selective routing

With measured baseline stage times `t_i` and explicit stage speedups `s_i`:

\[
T_0=\sum_i t_i,\qquad T_1=\sum_i t_i/s_i+T_{\mathrm{new\ overhead}}.
\]

If a stage is unchanged, supply `speedup: 1`. Omitted speedups remain unknown. The Amdahl limit for eliminating an accelerated baseline fraction `p` is `1/(1-p)`, holding remaining work fixed. Concurrent pipelines require a critical-path or throughput analysis instead of adding overlapping wall times.

For a learned first pass with conventional fallback:

\[
\alpha=T_{\mathrm{new\ front}}/T_0,\quad
f_t=\frac{\sum_{r\in\mathrm{fallback}}t_0(r)}{\sum_r t_0(r)},\quad
\mathrm{speedup}=\frac1{\alpha+f_t}.
\]

`front_end` must include all new non-fallback work, including candidate verification and relevant transfer/merge costs. `fallback_time_fraction` is **time weighted**. A small number of repetitive or very long reads may consume most fallback time. `fallback_read_fraction` is retained for diagnosis but never substituted for `f_t`.

For a caller-supplied goal `S`, feasibility requires `f_t <= 1/S - alpha`. A negative right side rejects the idea even with zero fallback. No speed target is baked into the module. Illustratively, `alpha=0.10`, `f_t=0.75` predicts only 1.18× speed, even if just 5% of reads fall back.

Keep both baseline and candidate indexes warm when measuring algorithmic gains. Reference-index persistence can improve deployment latency, but conventional aligners can reuse indexes too.

## 4. Accepted-route and candidate-miss error budgets

With accepted read fraction `a`, error rate `e_a` conditional on acceptance, and fallback-subset error rate `e_f`:

\[
e_{\mathrm{total}}=a e_a+(1-a)e_f.
\]

For a chosen total-error budget `epsilon`, the accepted-route allowance is `[epsilon-(1-a)e_f]/a` when `a>0`. A negative numerator means fallback alone exhausts the budget. Use the actual fallback subset's error rate; the overall baseline error rate is not interchangeable.

If candidate recall is `r` and the probability of accepting despite a candidate miss is `q`, the miss-induced wrong-acceptance contribution is `(1-r)*q`. This is only one error mechanism, not the total error, and must not be added again to the mixture above. `maximum_miss_error_probability` is a separate budget for this mechanism.

Cosine margins require calibration against missed alternative loci. Exact alignment establishes candidate sequence consistency, not correctness of genomic origin or high MAPQ.

## 5. Cache value depends on the workload

Let `l` be lookup time, `h` the hit rate, `t_h` and `t_m` the post-lookup hit/miss costs, `B` additional cache build cost, and `J` the amortization count:

\[
T_{\mathrm{cache}}=l+h t_h+(1-h)t_m+B/J.
\]

Against uncached cost `t_0`, when `t_m>t_h`, strict speedup requires:

\[
h>\frac{l+t_m+B/J-t_0}{t_m-t_h}.
\]

The module handles equal or more expensive hits without dividing by zero or assuming a favorable inequality. Unknown build cost stays unknown; explicitly supply `build_seconds: 0` when measuring only an already amortized cache.

Three hit-rate bounds have different assumptions:

- A cold trace of `n` requests and `D` distinct exact keys has at most `(n-D)/n` hits, even with unlimited capacity, because first occurrences miss.
- A stationary uniform space of `K` keys with `C` entries has hit rate at most `min(1,C/K)`.
- Under a supplied key-frequency histogram, an oracle **static** cache holds the `C` most frequent keys; its hit fraction is the sum of those counts divided by all counts. This is not a universal bound on an adaptive cache in a temporally changing trace.

Cache reference indexes by reference digest, tool/version and index parameters. Cache alignment results only with every behavior-affecting input represented in the key. Frequent reference regions do not imply identical noisy reads; exact-result cache reuse must be measured.

## 6. Evidence and confidence intervals

For `k` observed counted events in `n` trials, the module returns a two-sided [Wilson score interval](https://www.itl.nist.gov/div898/handbook/prc/section2/prc241.htm). `kind: "error"` means `successes` counts errors; `kind: "recall"` means it counts successful retrievals. For zero errors in 100 trials, the 95% upper bound is approximately 3.70%, not zero. Zero trials produce an unknown rate. These binomial intervals do not address related reads, donor correlation, selection bias, or distribution drift.

Evaluate the same reference, read set, scoring, output format and timing boundary. Report cold startup separately from warm throughput, and measure cost per completed sample as well as latency. Use known-origin simulated data and independent real samples. Real-read agreement with minimap2 is concordance, not ground truth.

Use [GIAB/GA4GH stratifications](https://github.com/genome-in-a-bottle/genome-stratifications) to inspect homopolymers, tandem repeats, segmental duplications, low mappability and difficult regions. Their BED resources support downstream true/false positive and false-negative analysis. Hold out complete genomic regions or samples when training: randomly splitting overlapping windows can overstate generalization.

## Python / CLI integration

```python
from genocache.model import analyze_model

report = analyze_model({
    "index": {
        "genome_bases": 3_100_000_000,
        "stride_bases": 1000,
        "window_bases": 1250,
        "max_query_reference_span": 250,
        "orientations": 1,
        "embedding_dim": 256,
        "dtype": "float16",
        "id_bytes": 8,
        "metadata_bytes_per_row": 16
    },
    "routing": {
        "front_end_fraction": 0.10,
        "fallback_time_fraction": 0.75,
        "fallback_read_fraction": 0.05,
        "target_speedup": 2
    },
    "union_recall": {"recall_a": 0.99, "recall_b": 0.99},
    "observations": {
        "accepted_mapping_errors": {
            "kind": "error", "successes": 0, "trials": 100,
            "confidence": 0.95
        }
    }
})
# Values above are examples, not measured GenoCache performance.
```

`analyze_model` accepts `ModelInputs` or a nested mapping. The exported frozen dataclasses document the complete schema: `IndexInputs`, `StageTiming`, `RoutingInputs`, `ErrorBudgetInputs`, `UnionRecallInputs`, `CacheInputs`, and `Observation`. Supported top-level keys are `index`, `stages`, `new_overhead_seconds`, `routing`, `error_budget`, `union_recall`, `cache`, and `observations`.

Probabilities must be in `[0,1]`, counts nonnegative integers, strides/dimensions/speedup factors positive, and supplied floating values finite. Unknown keys, contradictory totals, and incompatible PQ dimensions fail with `ValueError`. JSON serialization with `allow_nan=False` is supported. Infinite idealized speed has no finite number and is represented with `null` plus an explicit indicator.

Formula validation:

```bash
PYTHONPATH=production python -m unittest discover -s production/tests -p test_model.py -v
```

These tests verify calculations and input contracts. They are not evidence of whole-genome accuracy, GPU compatibility or a production throughput improvement.
