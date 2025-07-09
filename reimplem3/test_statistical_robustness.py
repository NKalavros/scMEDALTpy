#!/usr/bin/env python
"""
Test suite for statistical robustness improvements
"""

import unittest
import numpy as np
import pandas as pd
import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from medalt_statistical_robust import (
    EnrichmentTest, WilcoxonTest, MultipleTestingCorrection,
    BootstrapConfidenceInterval, PowerAnalysis, RobustStatisticalFramework,
    StatisticalResult
)

import logging
logging.getLogger().setLevel(logging.WARNING)  # Suppress debug logs


class TestStatisticalTests(unittest.TestCase):
    """Test individual statistical test implementations"""
    
    def setUp(self):
        """Set up test data"""
        np.random.seed(42)
        
        # Create test data with known properties
        self.enriched_data = np.random.normal(3.0, 1.0, 50)  # Mean=3, high values
        self.normal_data = np.random.normal(2.0, 1.0, 50)    # Mean=2, baseline
        self.depleted_data = np.random.normal(1.0, 1.0, 50)  # Mean=1, low values
        
        # Background data
        self.background_data = np.random.normal(2.0, 1.0, 100)
    
    def test_enrichment_test_basic(self):
        """Test basic enrichment test functionality"""
        test = EnrichmentTest()
        
        # Test with enriched data
        statistic = test.compute_statistic(self.enriched_data, self.background_data)
        effect_size = test.compute_effect_size(self.enriched_data, self.background_data)
        
        # Should be positive for enriched data
        self.assertGreater(statistic, 0)
        self.assertGreater(effect_size, 0)
        
        # Test with depleted data
        statistic_dep = test.compute_statistic(self.depleted_data, self.background_data)
        effect_size_dep = test.compute_effect_size(self.depleted_data, self.background_data)
        
        # Should be negative for depleted data
        self.assertLess(statistic_dep, 0)
        self.assertLess(effect_size_dep, 0)
    
    def test_enrichment_test_edge_cases(self):
        """Test enrichment test edge cases"""
        test = EnrichmentTest()
        
        # Identical data
        identical_data = np.full(20, 2.0)
        statistic = test.compute_statistic(identical_data, identical_data)
        effect_size = test.compute_effect_size(identical_data, identical_data)
        
        self.assertEqual(statistic, 0.0)
        self.assertEqual(effect_size, 0.0)
        
        # Very small samples
        small_obs = np.array([2.5, 2.6])
        small_bg = np.array([2.0, 2.1])
        
        statistic = test.compute_statistic(small_obs, small_bg)
        self.assertIsInstance(statistic, float)
    
    def test_wilcoxon_test_basic(self):
        """Test Wilcoxon test functionality"""
        test = WilcoxonTest()
        
        # Test with different distributions
        statistic = test.compute_statistic(self.enriched_data, self.background_data)
        effect_size = test.compute_effect_size(self.enriched_data, self.background_data)
        
        # Should detect difference
        self.assertNotEqual(statistic, 0.0)
        self.assertIsInstance(effect_size, float)
    
    def test_pvalue_calculation(self):
        """Test p-value calculation methods"""
        test = EnrichmentTest()
        
        # Create null distribution
        null_dist = np.random.normal(0, 1, 1000)
        
        # Test extreme positive value
        p_value_extreme = test.compute_pvalue(3.0, null_dist)
        self.assertLess(p_value_extreme, 0.05)
        
        # Test moderate value
        p_value_moderate = test.compute_pvalue(0.5, null_dist)
        self.assertGreater(p_value_moderate, 0.05)
        
        # Test with continuity correction
        test_no_correction = EnrichmentTest(continuity_correction=False)
        p_value_no_corr = test_no_correction.compute_pvalue(3.0, null_dist)
        
        # With correction should be slightly higher
        self.assertGreaterEqual(p_value_extreme, p_value_no_corr)


class TestMultipleTestingCorrection(unittest.TestCase):
    """Test multiple testing correction methods"""
    
    def setUp(self):
        """Set up test p-values"""
        self.pvalues = np.array([0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.8, 0.99])
    
    def test_bonferroni_correction(self):
        """Test Bonferroni correction"""
        corrected = MultipleTestingCorrection.bonferroni(self.pvalues)
        
        # Should be more conservative
        self.assertTrue(np.all(corrected >= self.pvalues))
        
        # Should not exceed 1.0
        self.assertTrue(np.all(corrected <= 1.0))
        
        # First p-value should be multiplied by n
        expected_first = min(self.pvalues[0] * len(self.pvalues), 1.0)
        self.assertAlmostEqual(corrected[0], expected_first, places=10)
    
    def test_holm_correction(self):
        """Test Holm-Bonferroni correction"""
        corrected = MultipleTestingCorrection.holm(self.pvalues)
        
        # Should be more conservative than original
        self.assertTrue(np.all(corrected >= self.pvalues))
        
        # Should not exceed 1.0
        self.assertTrue(np.all(corrected <= 1.0))
        
        # Should be less conservative than Bonferroni
        bonferroni = MultipleTestingCorrection.bonferroni(self.pvalues)
        self.assertTrue(np.all(corrected <= bonferroni))
    
    def test_fdr_bh_correction(self):
        """Test Benjamini-Hochberg FDR correction"""
        corrected = MultipleTestingCorrection.fdr_bh(self.pvalues)
        
        # Should be more conservative than original
        self.assertTrue(np.all(corrected >= self.pvalues))
        
        # Should not exceed 1.0
        self.assertTrue(np.all(corrected <= 1.0))
        
        # Should be less conservative than Bonferroni
        bonferroni = MultipleTestingCorrection.bonferroni(self.pvalues)
        self.assertTrue(np.all(corrected <= bonferroni))
    
    def test_fdr_by_correction(self):
        """Test Benjamini-Yekutieli FDR correction"""
        corrected = MultipleTestingCorrection.fdr_by(self.pvalues)
        
        # Should be more conservative than BH
        bh_corrected = MultipleTestingCorrection.fdr_bh(self.pvalues)
        self.assertTrue(np.all(corrected >= bh_corrected))
        
        # Should not exceed 1.0
        self.assertTrue(np.all(corrected <= 1.0))


class TestBootstrapConfidenceInterval(unittest.TestCase):
    """Test bootstrap confidence interval calculation"""
    
    def setUp(self):
        """Set up test data"""
        np.random.seed(42)
        self.data1 = np.random.normal(3.0, 1.0, 50)
        self.data2 = np.random.normal(2.0, 1.0, 50)
        
        def mean_diff(x, y):
            return np.mean(x) - np.mean(y)
        
        self.statistic_func = mean_diff
    
    def test_bootstrap_ci_basic(self):
        """Test basic bootstrap confidence interval"""
        bootstrap_ci = BootstrapConfidenceInterval(n_bootstrap=100, confidence_level=0.95)
        
        ci_lower, ci_upper = bootstrap_ci.compute_ci(
            self.data1, self.data2, self.statistic_func
        )
        
        # Should be valid interval
        self.assertLess(ci_lower, ci_upper)
        self.assertIsInstance(ci_lower, float)
        self.assertIsInstance(ci_upper, float)
        
        # Should contain the true difference approximately
        true_diff = np.mean(self.data1) - np.mean(self.data2)
        self.assertLessEqual(ci_lower, true_diff)
        self.assertGreaterEqual(ci_upper, true_diff)
    
    def test_bootstrap_ci_different_confidence_levels(self):
        """Test bootstrap CI with different confidence levels"""
        bootstrap_90 = BootstrapConfidenceInterval(n_bootstrap=100, confidence_level=0.90)
        bootstrap_95 = BootstrapConfidenceInterval(n_bootstrap=100, confidence_level=0.95)
        bootstrap_99 = BootstrapConfidenceInterval(n_bootstrap=100, confidence_level=0.99)
        
        ci_90 = bootstrap_90.compute_ci(self.data1, self.data2, self.statistic_func)
        ci_95 = bootstrap_95.compute_ci(self.data1, self.data2, self.statistic_func)
        ci_99 = bootstrap_99.compute_ci(self.data1, self.data2, self.statistic_func)
        
        # Higher confidence should give wider intervals
        width_90 = ci_90[1] - ci_90[0]
        width_95 = ci_95[1] - ci_95[0]
        width_99 = ci_99[1] - ci_99[0]
        
        self.assertLess(width_90, width_95)
        self.assertLess(width_95, width_99)


class TestPowerAnalysis(unittest.TestCase):
    """Test statistical power analysis"""
    
    def test_power_calculation_basic(self):
        """Test basic power calculation"""
        # Large effect size should have high power
        power_large = PowerAnalysis.compute_power(
            effect_size=1.0, n_obs=50, n_bg=50, alpha=0.05
        )
        
        # Small effect size should have low power
        power_small = PowerAnalysis.compute_power(
            effect_size=0.1, n_obs=50, n_bg=50, alpha=0.05
        )
        
        # Large effect should have higher power
        self.assertGreater(power_large, power_small)
        
        # Power should be between 0 and 1
        self.assertGreaterEqual(power_large, 0.0)
        self.assertLessEqual(power_large, 1.0)
        self.assertGreaterEqual(power_small, 0.0)
        self.assertLessEqual(power_small, 1.0)
    
    def test_power_with_sample_size(self):
        """Test power changes with sample size"""
        effect_size = 0.5
        
        # Small sample
        power_small = PowerAnalysis.compute_power(
            effect_size=effect_size, n_obs=10, n_bg=10, alpha=0.05
        )
        
        # Large sample
        power_large = PowerAnalysis.compute_power(
            effect_size=effect_size, n_obs=100, n_bg=100, alpha=0.05
        )
        
        # Larger sample should have higher power
        self.assertGreater(power_large, power_small)
    
    def test_power_approximation_fallback(self):
        """Test power approximation fallback"""
        # Test with very small samples (should use approximation)
        power = PowerAnalysis._approximate_power(
            effect_size=0.8, n_obs=5, n_bg=5, alpha=0.05
        )
        
        self.assertGreaterEqual(power, 0.0)
        self.assertLessEqual(power, 1.0)


class TestRobustStatisticalFramework(unittest.TestCase):
    """Test the complete robust statistical framework"""
    
    def setUp(self):
        """Set up test data"""
        np.random.seed(42)
        
        # Create test data
        self.observed_data = {
            'test1': np.random.normal(3.0, 1.0, 30),  # Enriched
            'test2': np.random.normal(2.0, 1.0, 25),  # Neutral
            'test3': np.random.normal(1.0, 1.0, 35),  # Depleted
        }
        
        self.background_data = {
            'test1': np.random.normal(2.0, 1.0, 50),
            'test2': np.random.normal(2.0, 1.0, 50),
            'test3': np.random.normal(2.0, 1.0, 50),
        }
        
        # Create null distributions
        self.null_distributions = {
            'test1': np.random.normal(0, 1, 500),
            'test2': np.random.normal(0, 1, 500),
            'test3': np.random.normal(0, 1, 500),
        }
    
    def test_framework_initialization(self):
        """Test framework initialization"""
        framework = RobustStatisticalFramework(
            test_method='enrichment',
            correction_method='fdr_bh',
            confidence_level=0.95,
            n_bootstrap=100
        )
        
        self.assertEqual(framework.test_method, 'enrichment')
        self.assertEqual(framework.correction_method, 'fdr_bh')
        self.assertEqual(framework.confidence_level, 0.95)
        self.assertEqual(framework.n_bootstrap, 100)
    
    def test_sequential_analysis(self):
        """Test sequential statistical analysis"""
        framework = RobustStatisticalFramework(
            test_method='enrichment',
            correction_method='fdr_bh',
            n_bootstrap=100,
            n_jobs=1  # Sequential
        )
        
        results = framework.run_statistical_analysis(
            self.observed_data, self.null_distributions, self.background_data
        )
        
        # Should have results for all tests
        self.assertEqual(len(results), 3)
        
        # Check required columns
        required_columns = ['test_id', 'statistic', 'pvalue', 'effect_size', 
                          'ci_lower', 'ci_upper', 'power', 'adjusted_pvalue']
        for col in required_columns:
            self.assertIn(col, results.columns)
        
        # Check significance indicators
        self.assertIn('significant', results.columns)
        self.assertIn('effect_size_magnitude', results.columns)
        
        # Results should be sorted by adjusted p-value
        adj_pvals = results['adjusted_pvalue'].values
        self.assertTrue(np.all(adj_pvals[:-1] <= adj_pvals[1:]))
    
    def test_parallel_analysis(self):
        """Test parallel statistical analysis"""
        framework = RobustStatisticalFramework(
            test_method='enrichment',
            correction_method='fdr_bh',
            n_bootstrap=100,
            n_jobs=2  # Parallel
        )
        
        results = framework.run_statistical_analysis(
            self.observed_data, self.null_distributions, self.background_data
        )
        
        # Should have same results as sequential
        self.assertEqual(len(results), 3)
        
        # Check that all tests have valid results
        self.assertTrue(np.all(results['pvalue'] >= 0.0))
        self.assertTrue(np.all(results['pvalue'] <= 1.0))
        self.assertTrue(np.all(results['adjusted_pvalue'] >= 0.0))
        self.assertTrue(np.all(results['adjusted_pvalue'] <= 1.0))
    
    def test_different_correction_methods(self):
        """Test different multiple testing correction methods"""
        correction_methods = ['bonferroni', 'holm', 'fdr_bh', 'fdr_by']
        
        results_dict = {}
        
        for method in correction_methods:
            framework = RobustStatisticalFramework(
                test_method='enrichment',
                correction_method=method,
                n_bootstrap=100,
                n_jobs=1
            )
            
            results = framework.run_statistical_analysis(
                self.observed_data, self.null_distributions, self.background_data
            )
            
            results_dict[method] = results
        
        # Check that different methods give different adjusted p-values
        for method1 in correction_methods:
            for method2 in correction_methods:
                if method1 != method2:
                    adj_pvals1 = results_dict[method1]['adjusted_pvalue'].values
                    adj_pvals2 = results_dict[method2]['adjusted_pvalue'].values
                    
                    # Should be different (unless by coincidence)
                    self.assertFalse(np.array_equal(adj_pvals1, adj_pvals2))
    
    def test_different_test_methods(self):
        """Test different statistical test methods"""
        test_methods = ['enrichment', 'wilcoxon']
        
        results_dict = {}
        
        for method in test_methods:
            framework = RobustStatisticalFramework(
                test_method=method,
                correction_method='fdr_bh',
                n_bootstrap=100,
                n_jobs=1
            )
            
            results = framework.run_statistical_analysis(
                self.observed_data, self.null_distributions, self.background_data
            )
            
            results_dict[method] = results
        
        # Both methods should produce results
        self.assertEqual(len(results_dict['enrichment']), 3)
        self.assertEqual(len(results_dict['wilcoxon']), 3)
        
        # Methods should be recorded correctly
        self.assertTrue(np.all(results_dict['enrichment']['method'] == 'enrichment'))
        self.assertTrue(np.all(results_dict['wilcoxon']['method'] == 'wilcoxon'))
    
    def test_edge_cases(self):
        """Test edge cases and error conditions"""
        framework = RobustStatisticalFramework(
            test_method='enrichment',
            correction_method='fdr_bh',
            n_bootstrap=100,
            min_observations=10
        )
        
        # Test with insufficient data
        small_data = {
            'test1': np.array([1.0, 2.0])  # Too small
        }
        
        small_background = {
            'test1': np.array([1.5, 2.5])  # Too small
        }
        
        small_null = {
            'test1': np.array([0.0, 0.1])
        }
        
        results = framework.run_statistical_analysis(
            small_data, small_null, small_background
        )
        
        # Should handle gracefully (return empty or skip)
        self.assertTrue(len(results) == 0)
        
        # Test with empty data
        empty_results = framework.run_statistical_analysis({}, {}, {})
        self.assertTrue(len(empty_results) == 0)


class TestStatisticalIntegration(unittest.TestCase):
    """Test integration of statistical framework with MEDALT"""
    
    def test_statistical_result_dataclass(self):
        """Test StatisticalResult dataclass"""
        result = StatisticalResult(
            statistic=2.5,
            pvalue=0.01,
            effect_size=0.8,
            confidence_interval=(0.2, 1.4),
            method='enrichment',
            n_observations=50,
            n_permutations=1000,
            power=0.85
        )
        
        self.assertEqual(result.statistic, 2.5)
        self.assertEqual(result.pvalue, 0.01)
        self.assertEqual(result.effect_size, 0.8)
        self.assertEqual(result.confidence_interval, (0.2, 1.4))
        self.assertEqual(result.method, 'enrichment')
        self.assertEqual(result.n_observations, 50)
        self.assertEqual(result.n_permutations, 1000)
        self.assertEqual(result.power, 0.85)
    
    def test_comprehensive_workflow(self):
        """Test comprehensive statistical workflow"""
        np.random.seed(42)
        
        # Create realistic test scenario
        n_tests = 100
        
        # Generate mixed data: some significant, some not
        observed_data = {}
        background_data = {}
        null_distributions = {}
        
        for i in range(n_tests):
            if i < 20:  # 20% truly significant
                obs_mean = 3.0 + np.random.normal(0, 0.5)
            else:  # 80% null
                obs_mean = 2.0 + np.random.normal(0, 0.2)
            
            observed_data[f'test_{i}'] = np.random.normal(obs_mean, 1.0, 30)
            background_data[f'test_{i}'] = np.random.normal(2.0, 1.0, 50)
            null_distributions[f'test_{i}'] = np.random.normal(0, 1, 500)
        
        # Run analysis
        framework = RobustStatisticalFramework(
            test_method='enrichment',
            correction_method='fdr_bh',
            confidence_level=0.95,
            n_bootstrap=100
        )
        
        results = framework.run_statistical_analysis(
            observed_data, null_distributions, background_data
        )
        
        # Should have results for all tests
        self.assertEqual(len(results), n_tests)
        
        # Should have reasonable number of significant results
        significant_count = np.sum(results['significant'])
        self.assertGreater(significant_count, 0)
        self.assertLess(significant_count, n_tests)  # Not all should be significant
        
        # FDR should control false discovery rate
        fdr_controlled = np.sum(results['adjusted_pvalue'] < 0.05)
        raw_significant = np.sum(results['pvalue'] < 0.05)
        
        # FDR should be more conservative
        self.assertLessEqual(fdr_controlled, raw_significant)


def run_performance_benchmark():
    """Run performance benchmark for statistical methods"""
    print("\n" + "="*50)
    print("STATISTICAL ROBUSTNESS PERFORMANCE BENCHMARK")
    print("="*50)
    
    import time
    
    # Create large dataset
    np.random.seed(42)
    n_tests = 1000
    
    observed_data = {}
    background_data = {}
    null_distributions = {}
    
    for i in range(n_tests):
        observed_data[f'test_{i}'] = np.random.normal(2.5, 1.0, 50)
        background_data[f'test_{i}'] = np.random.normal(2.0, 1.0, 100)
        null_distributions[f'test_{i}'] = np.random.normal(0, 1, 1000)
    
    # Test different configurations
    configs = [
        {'method': 'enrichment', 'correction': 'fdr_bh', 'bootstrap': 100},
        {'method': 'enrichment', 'correction': 'bonferroni', 'bootstrap': 100},
        {'method': 'wilcoxon', 'correction': 'fdr_bh', 'bootstrap': 100},
        {'method': 'enrichment', 'correction': 'fdr_bh', 'bootstrap': 1000},
    ]
    
    for config in configs:
        framework = RobustStatisticalFramework(
            test_method=config['method'],
            correction_method=config['correction'],
            n_bootstrap=config['bootstrap'],
            n_jobs=1
        )
        
        start_time = time.time()
        results = framework.run_statistical_analysis(
            observed_data, null_distributions, background_data
        )
        end_time = time.time()
        
        duration = end_time - start_time
        
        print(f"\nConfig: {config['method']}, {config['correction']}, "
              f"bootstrap={config['bootstrap']}")
        print(f"Time: {duration:.2f}s")
        print(f"Results: {len(results)} tests")
        print(f"Significant: {np.sum(results['significant'])}")
        print(f"Rate: {len(results)/duration:.1f} tests/sec")


if __name__ == '__main__':
    # Run unit tests
    print("Running Statistical Robustness Tests...")
    unittest.main(verbosity=2, exit=False)
    
    # Run performance benchmark
    run_performance_benchmark()