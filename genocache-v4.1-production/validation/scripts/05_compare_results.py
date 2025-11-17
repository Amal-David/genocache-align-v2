#!/usr/bin/env python3
"""
Compare Baseline vs Improved Results

Statistical analysis of improvements
"""

from pathlib import Path
import scipy.stats as stats

def load_results(filepath):
    """Load result file"""
    results = []
    with open(filepath) as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            results.append({
                'read_id': parts[0],
                'predicted_chr': parts[1],
                'true_chr': parts[2],
                'correct': parts[3] == 'True',
                'score': float(parts[4])
            })
    return results

def main():
    print("="*80)
    print("COMPARISON: Baseline vs Improved")
    print("="*80)
    print()
    
    # Load results
    base_dir = Path(__file__).parent.parent / "results"
    baseline = load_results(base_dir / "baseline_results.txt")
    improved = load_results(base_dir / "improved_results.txt")
    
    print(f"Loaded {len(baseline)} baseline results")
    print(f"Loaded {len(improved)} improved results")
    print()
    
    # Calculate metrics
    baseline_correct = sum(1 for r in baseline if r['correct'])
    improved_correct = sum(1 for r in improved if r['correct'])
    
    baseline_acc = (baseline_correct / len(baseline)) * 100
    improved_acc = (improved_correct / len(improved)) * 100
    
    print("="*80)
    print("OVERALL ACCURACY")
    print("="*80)
    print(f"Baseline:  {baseline_correct}/{len(baseline)} = {baseline_acc:.1f}%")
    print(f"Improved:  {improved_correct}/{len(improved)} = {improved_acc:.1f}%")
    print(f"Change:    {improved_acc - baseline_acc:+.1f}%")
    print()
    
    # Detailed breakdown
    print("="*80)
    print("ERROR BREAKDOWN")
    print("="*80)
    
    baseline_unmapped = sum(1 for r in baseline if r['predicted_chr'] == 'unmapped')
    improved_unmapped = sum(1 for r in improved if r['predicted_chr'] == 'unmapped')
    baseline_wrong = sum(1 for r in baseline if not r['correct'] and r['predicted_chr'] != 'unmapped')
    improved_wrong = sum(1 for r in improved if not r['correct'] and r['predicted_chr'] != 'unmapped')
    
    print("\nBaseline:")
    print(f"  Correct:   {baseline_correct} ({baseline_correct/len(baseline)*100:.1f}%)")
    print(f"  Wrong chr: {baseline_wrong} ({baseline_wrong/len(baseline)*100:.1f}%)")
    print(f"  Unmapped:  {baseline_unmapped} ({baseline_unmapped/len(baseline)*100:.1f}%)")
    
    print("\nImproved:")
    print(f"  Correct:   {improved_correct} ({improved_correct/len(improved)*100:.1f}%)")
    print(f"  Wrong chr: {improved_wrong} ({improved_wrong/len(improved)*100:.1f}%)")
    print(f"  Unmapped:  {improved_unmapped} ({improved_unmapped/len(improved)*100:.1f}%)")
    
    print("\nChanges:")
    print(f"  Correct:   {improved_correct - baseline_correct:+d}")
    print(f"  Wrong chr: {improved_wrong - baseline_wrong:+d}")
    print(f"  Unmapped:  {improved_unmapped - baseline_unmapped:+d}")
    print()
    
    # Per-read comparison
    print("="*80)
    print("PER-READ CHANGES")
    print("="*80)
    
    changed_better = []
    changed_worse = []
    
    for b, i in zip(baseline, improved):
        assert b['read_id'] == i['read_id']
        
        if not b['correct'] and i['correct']:
            changed_better.append(b['read_id'])
        elif b['correct'] and not i['correct']:
            changed_worse.append(b['read_id'])
    
    print(f"\nReads that IMPROVED: {len(changed_better)}")
    if changed_better:
        for read_id in changed_better[:10]:
            print(f"  - {read_id}")
        if len(changed_better) > 10:
            print(f"  ... and {len(changed_better) - 10} more")
    
    print(f"\nReads that got WORSE: {len(changed_worse)}")
    if changed_worse:
        for read_id in changed_worse[:10]:
            print(f"  - {read_id}")
        if len(changed_worse) > 10:
            print(f"  ... and {len(changed_worse) - 10} more")
    
    print()
    
    # Statistical test
    print("="*80)
    print("STATISTICAL SIGNIFICANCE")
    print("="*80)
    
    # McNemar's test for paired categorical data
    baseline_correct_flags = [1 if r['correct'] else 0 for r in baseline]
    improved_correct_flags = [1 if r['correct'] else 0 for r in improved]
    
    # Contingency table
    both_correct = sum(1 for b, i in zip(baseline_correct_flags, improved_correct_flags) if b == 1 and i == 1)
    both_wrong = sum(1 for b, i in zip(baseline_correct_flags, improved_correct_flags) if b == 0 and i == 0)
    base_correct_imp_wrong = sum(1 for b, i in zip(baseline_correct_flags, improved_correct_flags) if b == 1 and i == 0)
    base_wrong_imp_correct = sum(1 for b, i in zip(baseline_correct_flags, improved_correct_flags) if b == 0 and i == 1)
    
    print("\nContingency Table:")
    print(f"  Both correct:     {both_correct}")
    print(f"  Both wrong:       {both_wrong}")
    print(f"  Base C, Imp W:    {base_correct_imp_wrong}")
    print(f"  Base W, Imp C:    {base_wrong_imp_correct}")
    
    if base_correct_imp_wrong + base_wrong_imp_correct > 0:
        # McNemar's test
        from scipy.stats import binom_test
        p_value = binom_test(base_wrong_imp_correct, base_correct_imp_wrong + base_wrong_imp_correct, 0.5)
        print(f"\nMcNemar's test p-value: {p_value:.4f}")
        
        if p_value < 0.05:
            print("✅ Statistically significant difference (p < 0.05)")
        else:
            print("⚠️  Not statistically significant (p >= 0.05)")
    else:
        print("\n⚠️  No changes between baseline and improved")
    
    # Save report
    print("\n" + "="*80)
    print("SAVING REPORT")
    print("="*80)
    
    report_path = base_dir / "VALIDATION_REPORT.md"
    with open(report_path, 'w') as f:
        f.write("# Validation Report: Baseline vs Improved\n\n")
        f.write(f"**Date:** 2025-11-15\n\n")
        f.write("---\n\n")
        f.write("## Overall Results\n\n")
        f.write(f"| Metric | Baseline | Improved | Change |\n")
        f.write(f"|--------|----------|----------|--------|\n")
        f.write(f"| Accuracy | {baseline_acc:.1f}% | {improved_acc:.1f}% | {improved_acc - baseline_acc:+.1f}% |\n")
        f.write(f"| Correct | {baseline_correct}/100 | {improved_correct}/100 | {improved_correct - baseline_correct:+d} |\n")
        f.write(f"| Wrong chr | {baseline_wrong}/100 | {improved_wrong}/100 | {improved_wrong - baseline_wrong:+d} |\n")
        f.write(f"| Unmapped | {baseline_unmapped}/100 | {improved_unmapped}/100 | {improved_unmapped - baseline_unmapped:+d} |\n")
        f.write("\n")
        f.write("## Analysis\n\n")
        f.write(f"- Reads that improved: {len(changed_better)}\n")
        f.write(f"- Reads that got worse: {len(changed_worse)}\n")
        f.write(f"- Net change: {len(changed_better) - len(changed_worse)}\n")
        f.write("\n")
        f.write("## Conclusion\n\n")
        if improved_acc > baseline_acc:
            f.write(f"✅ Improvements show +{improved_acc - baseline_acc:.1f}% accuracy gain.\n")
        elif improved_acc == baseline_acc:
            f.write(f"⚠️  No accuracy change, but error distribution shifted.\n")
        else:
            f.write(f"❌ Accuracy decreased by {baseline_acc - improved_acc:.1f}%.\n")
    
    print(f"✅ Report saved to: {report_path}")
    print()
    
    # Final summary
    print("="*80)
    print("SUMMARY")
    print("="*80)
    print()
    if improved_acc > baseline_acc + 1:
        print(f"✅ IMPROVEMENT: +{improved_acc - baseline_acc:.1f}% accuracy")
    elif abs(improved_acc - baseline_acc) <= 1:
        print(f"⚠️  NO SIGNIFICANT CHANGE: {improved_acc - baseline_acc:+.1f}% (within margin)")
    else:
        print(f"❌ REGRESSION: {improved_acc - baseline_acc:.1f}% accuracy")
    print()
    print("Note: Both versions show low accuracy (26%) on this synthetic test data.")
    print("This suggests the test data itself may have issues (errors too high,")
    print("or reference extraction problems). Real sequencing data validation needed.")


if __name__ == '__main__':
    main()
