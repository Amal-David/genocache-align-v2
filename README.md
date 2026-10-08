# GenoCache

The current runnable implementation is in **[production/](production/README.md)**.
It provides verified minimap2 batch jobs, reusable reference packs, exact-input
result caching, Modal/AWS workers, and an independently evaluated embedding
retrieval path.

The other version directories are the November 2025 research backup. Their
historical accuracy and speed claims are not acceptance criteria for this
implementation. The [source audit](production/docs/LEGACY_AUDIT.md) reproduces
specific score, strand, coordinate, chaining, output, and timing defects.

- [Install and run](production/README.md)
- [Cloud deployment](production/docs/DEPLOYMENT.md)
- [Codex operating and experiment runbook](production/docs/CODEX_RUNBOOK.md)
- [Mathematical model and research](production/docs/RESEARCH_MODEL.md)
- [Encoder artifact contract](production/docs/ENCODER_CONTRACT.md)
- [Measured synthetic integration receipt](production/reports/local-benchmark.json)

New work should use the `production/` package. Its documentation distinguishes
implemented behavior, measured evidence, and experiments that still need to pass.
