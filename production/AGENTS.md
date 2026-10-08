# GenoCache implementation guidance

Work in this package; the sibling version directories are archived experiments.
Read `README.md` and `docs/CODEX_RUNBOOK.md` before changing the pipeline.

## Contracts to preserve

- Native minimap2 owns final SAM records. Preserve flags, read names, qualities,
  coordinates, CIGAR, MAPQ, secondary/supplementary records and tags.
- Embedding similarity and region proposals are candidate evidence. They do not
  establish genomic uniqueness or authorize a fabricated alignment record.
- Reference packs and cached jobs bind exact input bytes, native binary,
  configuration and output artifacts. Keep verification and atomic publication.
- Use worker-local POSIX storage for locks and cache roots. Cloud attempts have
  separate prefixes; the completed result manifest is uploaded last.
- Model/index identity includes the actual checkpoint and preprocessing contract.
  Missing or incompatible weights are errors, not a reason to use random weights.
- Never include private source, checkpoints, credentials or real genomic samples
  in this public repository. Small deterministic synthetic fixtures are suitable.

## Evidence-driven changes

Start from one causal hypothesis and a bounded experiment. Separate exact vector
ranking, ANN approximation, biological locus recall, region retention, final
mapping and end-to-end time. Keep misses and unmapped reads in denominators.
Record full preceding rescue costs and time-weighted fallback work. Use the
mathematical model and declared accuracy/speed gates to close failed branches.
Do not change a success tolerance after inspecting results to claim a win.

Run focused tests for the changed contracts. The acceptance suite is `pytest
tests -q`; lint is `ruff check .`. Optional backend skips are not verification of
that backend. The synthetic benchmark is `scripts/benchmark_local.py` and requires
a new empty work directory. Cloud deployment still requires image construction
and an authenticated smoke job in the target account. Report what actually ran.
