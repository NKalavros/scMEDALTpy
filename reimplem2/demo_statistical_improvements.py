#!/usr/bin/env python
"""
Demonstration of statistical robustness improvements in MEDALT
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Tuple
import logging

# Import components
from medalt_statistical_robust import (
    RobustStatisticalFramework, EnrichmentTest, WilcoxonTest,
    MultipleTestingCorrection, PowerAnalysis
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def demonstrate_statistical_improvements():
    """Demonstrate the key statistical improvements"""
    
    print("="*60)
    print("MEDALT STATISTICAL ROBUSTNESS IMPROVEMENTS")
    print("="*60)
    
    # 1. Demonstrate p-value calculation improvements
    print("\n1. P-VALUE CALCULATION IMPROVEMENTS")
    print("-" * 40)
    
    # Create test data with known effect
    np.random.seed(42)
    enriched_data = np.random.normal(3.0, 1.0, 50)
    background_data = np.random.normal(2.0, 1.0, 100)
    null_distribution = np.random.normal(0, 1, 1000)
    
    # Test original simple p-value vs improved
    test = EnrichmentTest(continuity_correction=True)
    statistic = test.compute_statistic(enriched_data, background_data)
    
    # Simple empirical p-value (original approach)
    n_extreme = np.sum(null_distribution >= statistic)
    simple_pvalue = n_extreme / len(null_distribution)
    
    # Improved p-value with continuity correction
    improved_pvalue = test.compute_pvalue(statistic, null_distribution)
    
    print(f"Test statistic: {statistic:.3f}")
    print(f"Simple p-value: {simple_pvalue:.6f}")
    print(f"Improved p-value: {improved_pvalue:.6f}")
    print(f"Improvement: {improved_pvalue/simple_pvalue:.2f}x more conservative")
    
    # 2. Demonstrate multiple testing corrections
    print("\n2. MULTIPLE TESTING CORRECTIONS")
    print("-" * 40)
    
    # Create multiple p-values
    raw_pvalues = np.array([0.001, 0.01, 0.02, 0.05, 0.08, 0.1, 0.2, 0.5])
    
    corrections = {
        'Bonferroni': MultipleTestingCorrection.bonferroni,
        'Holm': MultipleTestingCorrection.holm,
        'FDR (BH)': MultipleTestingCorrection.fdr_bh,
        'FDR (BY)': MultipleTestingCorrection.fdr_by
    }
    
    results_table = pd.DataFrame({'Raw p-value': raw_pvalues})
    
    for name, correction_func in corrections.items():
        corrected = correction_func(raw_pvalues)
        results_table[name] = corrected
    
    print(results_table.round(4))
    
    # 3. Demonstrate effect size and confidence intervals
    print("\n3. EFFECT SIZE AND CONFIDENCE INTERVALS")
    print("-" * 40)
    
    # Different effect sizes
    effect_scenarios = [
        ("Small effect", np.random.normal(2.2, 1.0, 30)),
        ("Medium effect", np.random.normal(2.5, 1.0, 30)),
        ("Large effect", np.random.normal(3.0, 1.0, 30))
    ]
    
    background = np.random.normal(2.0, 1.0, 50)
    
    for name, data in effect_scenarios:
        effect_size = test.compute_effect_size(data, background)
        
        # Calculate confidence interval (simplified)
        bootstrap_effects = []
        for _ in range(100):
            boot_data = np.random.choice(data, len(data), replace=True)
            boot_bg = np.random.choice(background, len(background), replace=True)
            boot_effect = test.compute_effect_size(boot_data, boot_bg)
            bootstrap_effects.append(boot_effect)
        
        ci_lower = np.percentile(bootstrap_effects, 2.5)
        ci_upper = np.percentile(bootstrap_effects, 97.5)
        
        print(f"{name}: Effect size = {effect_size:.3f}, "
              f"95% CI = [{ci_lower:.3f}, {ci_upper:.3f}]")
    
    # 4. Demonstrate power analysis
    print("\n4. STATISTICAL POWER ANALYSIS")
    print("-" * 40)
    
    effect_sizes = [0.2, 0.5, 0.8, 1.0, 1.5]
    sample_sizes = [20, 50, 100, 200]
    
    power_results = []
    for effect in effect_sizes:
        for n in sample_sizes:
            power = PowerAnalysis.compute_power(effect, n, n, alpha=0.05)
            power_results.append({
                'effect_size': effect,
                'sample_size': n,
                'power': power
            })
    
    power_df = pd.DataFrame(power_results)
    power_pivot = power_df.pivot(index='effect_size', columns='sample_size', values='power')
    
    print("Statistical Power (rows=effect size, cols=sample size):")
    print(power_pivot.round(3))
    
    # 5. Demonstrate comprehensive framework
    print("\n5. COMPREHENSIVE STATISTICAL FRAMEWORK")
    print("-" * 40)
    
    # Create realistic test scenario
    np.random.seed(42)
    n_tests = 20
    
    observed_data = {}
    background_data = {}
    null_distributions = {}
    
    for i in range(n_tests):
        if i < 5:  # 25% truly significant
            obs_mean = 3.0 + np.random.normal(0, 0.3)
        else:  # 75% null
            obs_mean = 2.0 + np.random.normal(0, 0.1)
        
        observed_data[f'test_{i}'] = np.random.normal(obs_mean, 1.0, 40)
        background_data[f'test_{i}'] = np.random.normal(2.0, 1.0, 60)
        null_distributions[f'test_{i}'] = np.random.normal(0, 1, 500)
    
    # Run robust analysis
    framework = RobustStatisticalFramework(
        test_method='enrichment',
        correction_method='fdr_bh',
        confidence_level=0.95,
        n_bootstrap=100
    )
    
    results = framework.run_statistical_analysis(
        observed_data, null_distributions, background_data
    )
    
    print(f"Total tests: {len(results)}")
    print(f"Significant (raw p < 0.05): {np.sum(results['pvalue'] < 0.05)}")
    print(f"Significant (FDR < 0.05): {np.sum(results['significant'])}")
    print(f"Mean effect size: {results['effect_size'].mean():.3f}")
    print(f"Mean statistical power: {results['power'].mean():.3f}")
    
    # Show top results
    print("\nTop 5 most significant results:")
    top_results = results.head()[['test_id', 'pvalue', 'adjusted_pvalue', 
                                  'effect_size', 'power', 'significant']]
    print(top_results.round(4))


def compare_statistical_approaches():
    """Compare original vs improved statistical approaches"""
    
    print("\n" + "="*60)
    print("COMPARISON: ORIGINAL vs IMPROVED STATISTICS")
    print("="*60)
    
    # Create test data
    np.random.seed(42)
    enriched_data = np.random.normal(2.8, 1.0, 30)
    background_data = np.random.normal(2.0, 1.0, 50)
    
    # Original approach (simple)
    print("\nORIGINAL APPROACH:")
    print("-" * 20)
    
    # Simple mean difference
    original_statistic = np.mean(enriched_data) - np.mean(background_data)
    
    # Simple empirical p-value
    permutation_stats = []
    for _ in range(1000):
        combined = np.concatenate([enriched_data, background_data])
        np.random.shuffle(combined)
        perm_enr = combined[:len(enriched_data)]
        perm_bg = combined[len(enriched_data):]
        perm_stat = np.mean(perm_enr) - np.mean(perm_bg)
        permutation_stats.append(perm_stat)
    
    original_pvalue = np.mean(np.array(permutation_stats) >= original_statistic)
    
    # Simple Bonferroni correction (assuming 100 tests)
    original_corrected = min(float(original_pvalue * 100), 1.0)
    
    print(f"Statistic: {original_statistic:.3f}")
    print(f"P-value: {original_pvalue:.6f}")
    print(f"Bonferroni corrected: {original_corrected:.6f}")
    print(f"Effect size: Not calculated")
    print(f"Confidence interval: Not calculated")
    print(f"Power: Not calculated")
    
    # Improved approach
    print("\nIMPROVED APPROACH:")
    print("-" * 20)
    
    test = EnrichmentTest(continuity_correction=True)
    improved_statistic = test.compute_statistic(enriched_data, background_data)
    improved_pvalue = test.compute_pvalue(improved_statistic, np.array(permutation_stats))
    effect_size = test.compute_effect_size(enriched_data, background_data)
    
    # FDR correction
    pvalues_for_correction = np.array([improved_pvalue] + [0.1] * 99)  # Simulate 100 tests
    fdr_corrected = MultipleTestingCorrection.fdr_bh(pvalues_for_correction)[0]
    
    # Power calculation
    power = PowerAnalysis.compute_power(effect_size, len(enriched_data), len(background_data))
    
    print(f"Statistic: {improved_statistic:.3f}")
    print(f"P-value: {improved_pvalue:.6f}")
    print(f"FDR corrected: {fdr_corrected:.6f}")
    print(f"Effect size: {effect_size:.3f}")
    print(f"Power: {power:.3f}")
    
    # Bootstrap confidence interval
    bootstrap_effects = []
    for _ in range(100):
        boot_enr = np.random.choice(enriched_data, len(enriched_data), replace=True)
        boot_bg = np.random.choice(background_data, len(background_data), replace=True)
        boot_effect = test.compute_effect_size(boot_enr, boot_bg)
        bootstrap_effects.append(boot_effect)
    
    ci_lower = np.percentile(bootstrap_effects, 2.5)
    ci_upper = np.percentile(bootstrap_effects, 97.5)
    print(f"95% CI: [{ci_lower:.3f}, {ci_upper:.3f}]")
    
    # Summary comparison
    print("\nSUMMARY IMPROVEMENTS:")
    print("-" * 30)
    print("✅ Continuity correction for p-values")
    print("✅ Multiple testing correction (FDR vs Bonferroni)")
    print("✅ Effect size calculation (Cohen's d)")
    print("✅ Bootstrap confidence intervals")
    print("✅ Statistical power analysis")
    print("✅ Comprehensive result reporting")


def demonstrate_multiple_testing_problem():
    """Demonstrate the multiple testing problem and corrections"""
    
    print("\n" + "="*60)
    print("MULTIPLE TESTING PROBLEM DEMONSTRATION")
    print("="*60)
    
    # Simulate multiple testing scenario
    np.random.seed(42)
    n_tests = 100
    
    # All null hypothesis (no real effect)
    null_pvalues = np.random.uniform(0, 1, n_tests)
    
    print(f"Scenario: {n_tests} tests, all null hypothesis")
    print(f"Expected false positives at α=0.05: {n_tests * 0.05:.1f}")
    print(f"Actual false positives: {np.sum(null_pvalues < 0.05)}")
    
    # Apply corrections
    corrections = {
        'No correction': null_pvalues,
        'Bonferroni': MultipleTestingCorrection.bonferroni(null_pvalues),
        'Holm': MultipleTestingCorrection.holm(null_pvalues),
        'FDR (BH)': MultipleTestingCorrection.fdr_bh(null_pvalues),
        'FDR (BY)': MultipleTestingCorrection.fdr_by(null_pvalues)
    }
    
    print("\nFalse positives with different corrections:")
    for name, corrected_pvals in corrections.items():
        false_positives = np.sum(corrected_pvals < 0.05)
        print(f"{name:12}: {false_positives:2d} false positives")
    
    # Mixed scenario: some true effects
    print(f"\nMixed scenario: 10 true effects, 90 null")
    
    mixed_pvalues = np.concatenate([
        np.random.beta(0.1, 1, 10),  # True effects (small p-values)
        np.random.uniform(0, 1, 90)   # Null effects
    ])
    
    print("\nDetected positives with different corrections:")
    
    # Define correction functions
    correction_funcs = {
        'No correction': lambda x: x,
        'Bonferroni': MultipleTestingCorrection.bonferroni,
        'Holm': MultipleTestingCorrection.holm,
        'FDR (BH)': MultipleTestingCorrection.fdr_bh,
        'FDR (BY)': MultipleTestingCorrection.fdr_by
    }
    
    for name, correction_func in correction_funcs.items():
        corrected_pvals = correction_func(mixed_pvalues)
        
        positives = np.sum(corrected_pvals < 0.05)
        true_positives = np.sum(corrected_pvals[:10] < 0.05)  # First 10 are true
        false_positives = np.sum(corrected_pvals[10:] < 0.05)  # Rest are false
        
        print(f"{name:12}: {positives:2d} total ({true_positives:2d} true, {false_positives:2d} false)")


def main():
    """Run all demonstrations"""
    
    print("MEDALT Statistical Robustness Improvements")
    print("This demonstration shows the enhanced statistical methods")
    print("implemented to improve p-value calculations and significance testing\n")
    
    # Run demonstrations
    demonstrate_statistical_improvements()
    compare_statistical_approaches()
    demonstrate_multiple_testing_problem()
    
    print("\n" + "="*60)
    print("SUMMARY OF IMPROVEMENTS")
    print("="*60)
    
    improvements = [
        "🔬 **Robust P-value Calculations**",
        "   - Continuity correction for empirical p-values",
        "   - Proper handling of edge cases and small samples",
        "   - Two-sided vs one-sided test options",
        "",
        "📊 **Advanced Effect Size Metrics**",
        "   - Cohen's d for parametric tests",
        "   - Rank-biserial correlation for non-parametric tests",
        "   - Standardized effect size interpretations",
        "",
        "🎯 **Multiple Testing Corrections**",
        "   - Bonferroni (conservative)",
        "   - Holm-Bonferroni (less conservative)",
        "   - FDR Benjamini-Hochberg (balanced)",
        "   - FDR Benjamini-Yekutieli (dependent tests)",
        "",
        "📈 **Bootstrap Confidence Intervals**",
        "   - Non-parametric confidence intervals",
        "   - Configurable confidence levels",
        "   - Robust to distribution assumptions",
        "",
        "⚡ **Statistical Power Analysis**",
        "   - Post-hoc power calculations",
        "   - Effect size and sample size relationships",
        "   - Study design optimization insights",
        "",
        "🔧 **Comprehensive Framework**",
        "   - Modular test implementations",
        "   - Parallel processing support",
        "   - Extensive validation and testing",
        "   - Clear result interpretation"
    ]
    
    for improvement in improvements:
        print(improvement)
    
    print("\n" + "="*60)
    print("🎉 STATISTICAL ROBUSTNESS IMPLEMENTATION COMPLETE!")
    print("="*60)


if __name__ == '__main__':
    main()