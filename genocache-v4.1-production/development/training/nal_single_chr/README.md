# Single Chromosome NAL Experiment

## Goal
Train NAL on ONE chromosome (Chr22) only to verify:
1. Can we achieve high accuracy on a controlled subset?
2. Does curriculum learning help?
3. Do adaptive seed sizes (256bp + 512bp) improve coverage?

## Parameters
- Chromosome: Chr22 (smallest autosome, ~50Mbp)
- Batches: 4000-8000 (vs 62 in original)
- Curriculum Learning: Start 5% error → 15% error
- Adaptive Seeds: 256bp + 512bp
- Model: 128D embeddings

## Expected Outcome
If training/index are correct:
- Should achieve >95% accuracy on Chr22 reads
- minimap2 comparison should show <10% FP rate

If still poor accuracy:
- Indicates fundamental model/approach issue
- Need different architecture or method
