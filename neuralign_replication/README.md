# Neuralign replication — GenoCache Align team workspace

## Immediate goal
1. Reproduce the core experiments and codebase of *Neuralign* (replicate results).
2. Produce a minimal runnable demo on Nebius H100.
3. After replication, iterate to integrate vector-index seeding & long-read mapping.

## Current status
- Created workspace on Nebius VM and S3 bucket: s3://hackathon-team-13/neuralign_replication
- Owner / contact: Siril Arockiam (arockiam)

## Collaboration
- Use this folder on the Nebius VM for code and small artifacts.
- Use `aws s3 sync` to share larger artifacts (models/indexes) with teammates.

## Next steps (to be executed after team sync)
- Pull Neuralign code/paper & identify core components to run (model, datasets, eval).
- Create reproducible environment (environment.yml provided).
