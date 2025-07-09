#!/usr/bin/env python
"""
MEDALT Statistical Robustness Framework
Advanced statistical methods for p-value calculations and significance testing
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, Union
import logging
from dataclasses import dataclass
from abc import ABC, abstractmethod
import warnings
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp

# Optional imports with fallbacks
try:
    from scipy import stats
    from scipy.stats import bootstrap
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

try:
    from statsmodels.stats.multitest import multipletests
    HAS_STATSMODELS = True
except ImportError:
    HAS_STATSMODELS = False

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class StatisticalResult:
    """Container for statistical test results"""
    statistic: float
    pvalue: float
    effect_size: float
    confidence_interval: Tuple[float, float]
    method: str
    n_observations: int
    n_permutations: int
    power: Optional[float] = None
    adjusted_pvalue: Optional[float] = None


class StatisticalTest(ABC):
    """Abstract base class for statistical tests"""
    
    @abstractmethod
    def compute_statistic(self, observed_data: np.ndarray, 
                         background_data: np.ndarray) -> float:
        """Compute the test statistic"""
        pass
    
    @abstractmethod
    def compute_pvalue(self, observed_statistic: float, 
                      null_distribution: np.ndarray) -> float:
        """Compute p-value from null distribution"""
        pass
    
    @abstractmethod
    def compute_effect_size(self, observed_data: np.ndarray,
                           background_data: np.ndarray) -> float:
        """Compute effect size"""
        pass


class EnrichmentTest(StatisticalTest):
    """Enrichment test with robust statistics"""
    
    def __init__(self, continuity_correction: bool = True):
        self.continuity_correction = continuity_correction
    
    def compute_statistic(self, observed_data: np.ndarray, 
                         background_data: np.ndarray) -> float:
        """Compute enrichment statistic"""
        obs_mean = np.mean(observed_data)
        bg_mean = np.mean(background_data)
        
        # Use standardized difference
        pooled_std = np.sqrt(
            (np.var(observed_data) + np.var(background_data)) / 2
        )
        
        if pooled_std == 0:
            return 0.0
        
        return (obs_mean - bg_mean) / pooled_std
    
    def compute_pvalue(self, observed_statistic: float, 
                      null_distribution: np.ndarray) -> float:
        """Compute p-value with continuity correction"""
        if len(null_distribution) == 0:
            return 1.0
        
        # Two-sided test
        n_extreme = np.sum(np.abs(null_distribution) >= np.abs(observed_statistic))
        n_total = len(null_distribution)
        
        if self.continuity_correction:
            # Add continuity correction
            p_value = (n_extreme + 1) / (n_total + 1)
        else:
            p_value = n_extreme / n_total
        
        return p_value
    
    def compute_effect_size(self, observed_data: np.ndarray,
                           background_data: np.ndarray) -> float:
        """Compute Cohen's d effect size"""
        obs_mean = np.mean(observed_data)
        bg_mean = np.mean(background_data)
        
        # Pooled standard deviation
        n1, n2 = len(observed_data), len(background_data)
        pooled_std = np.sqrt(
            ((n1 - 1) * np.var(observed_data, ddof=1) + 
             (n2 - 1) * np.var(background_data, ddof=1)) / (n1 + n2 - 2)
        )
        
        if pooled_std == 0:
            return 0.0
        
        return (obs_mean - bg_mean) / pooled_std


class WilcoxonTest(StatisticalTest):
    """Wilcoxon rank-sum test (non-parametric)"""
    
    def compute_statistic(self, observed_data: np.ndarray, 
                         background_data: np.ndarray) -> float:
        """Compute Wilcoxon statistic"""
        if not HAS_SCIPY:
            # Fallback to simple mean difference
            return np.mean(observed_data) - np.mean(background_data)
        
        try:
            statistic, _ = stats.ranksums(observed_data, background_data)
            return statistic
        except Exception:
            return np.mean(observed_data) - np.mean(background_data)
    
    def compute_pvalue(self, observed_statistic: float, 
                      null_distribution: np.ndarray) -> float:
        """Compute p-value using normal approximation"""
        if len(null_distribution) == 0:
            return 1.0
        
        # Use normal approximation for large samples
        mean_null = np.mean(null_distribution)
        std_null = np.std(null_distribution)
        
        if std_null == 0:
            return 1.0 if observed_statistic == mean_null else 0.0
        
        z_score = (observed_statistic - mean_null) / std_null
        
        if HAS_SCIPY:
            p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))
        else:
            # Fallback approximation
            p_value = 2 * (1 - self._normal_cdf(abs(z_score)))
        
        return p_value
    
    def compute_effect_size(self, observed_data: np.ndarray,
                           background_data: np.ndarray) -> float:
        """Compute rank-biserial correlation as effect size"""
        if not HAS_SCIPY:
            return self._cohens_d_fallback(observed_data, background_data)
        
        try:
            statistic, _ = stats.ranksums(observed_data, background_data)
            n1, n2 = len(observed_data), len(background_data)
            
            # Rank-biserial correlation
            r = statistic / np.sqrt(n1 * n2 * (n1 + n2) / 3)
            return r
        except Exception:
            return self._cohens_d_fallback(observed_data, background_data)
    
    def _cohens_d_fallback(self, obs: np.ndarray, bg: np.ndarray) -> float:
        """Fallback Cohen's d calculation"""
        pooled_std = np.sqrt((np.var(obs) + np.var(bg)) / 2)
        if pooled_std == 0:
            return 0.0
        return (np.mean(obs) - np.mean(bg)) / pooled_std
    
    def _normal_cdf(self, x: float) -> float:
        """Fallback normal CDF approximation"""
        return 0.5 * (1 + np.sign(x) * np.sqrt(1 - np.exp(-2 * x * x / np.pi)))


class BootstrapConfidenceInterval:
    """Bootstrap confidence interval calculation"""
    
    def __init__(self, n_bootstrap: int = 1000, confidence_level: float = 0.95):
        self.n_bootstrap = n_bootstrap
        self.confidence_level = confidence_level
    
    def compute_ci(self, observed_data: np.ndarray, 
                   background_data: np.ndarray,
                   statistic_func) -> Tuple[float, float]:
        """Compute bootstrap confidence interval"""
        
        if HAS_SCIPY and len(observed_data) > 10 and len(background_data) > 10:
            try:
                return self._scipy_bootstrap(observed_data, background_data, statistic_func)
            except Exception:
                return self._manual_bootstrap(observed_data, background_data, statistic_func)
        else:
            return self._manual_bootstrap(observed_data, background_data, statistic_func)
    
    def _scipy_bootstrap(self, obs: np.ndarray, bg: np.ndarray, 
                        statistic_func) -> Tuple[float, float]:
        """Use scipy bootstrap if available"""
        def stat_func(x, y):
            return statistic_func(x, y)
        
        # Bootstrap resampling
        res = bootstrap((obs, bg), stat_func, n_resamples=self.n_bootstrap, 
                       confidence_level=self.confidence_level)
        
        return res.confidence_interval.low, res.confidence_interval.high
    
    def _manual_bootstrap(self, obs: np.ndarray, bg: np.ndarray,
                         statistic_func) -> Tuple[float, float]:
        """Manual bootstrap implementation"""
        bootstrap_stats = []
        
        for _ in range(self.n_bootstrap):
            # Resample with replacement
            obs_boot = np.random.choice(obs, size=len(obs), replace=True)
            bg_boot = np.random.choice(bg, size=len(bg), replace=True)
            
            # Compute statistic
            stat = statistic_func(obs_boot, bg_boot)
            bootstrap_stats.append(stat)
        
        # Calculate confidence interval
        alpha = 1 - self.confidence_level
        lower_percentile = 100 * (alpha / 2)
        upper_percentile = 100 * (1 - alpha / 2)
        
        ci_lower = np.percentile(bootstrap_stats, lower_percentile)
        ci_upper = np.percentile(bootstrap_stats, upper_percentile)
        
        return ci_lower, ci_upper


class MultipleTestingCorrection:
    """Multiple testing correction methods"""
    
    @staticmethod
    def bonferroni(pvalues: np.ndarray) -> np.ndarray:
        """Bonferroni correction"""
        return np.minimum(pvalues * len(pvalues), 1.0)
    
    @staticmethod
    def holm(pvalues: np.ndarray) -> np.ndarray:
        """Holm-Bonferroni correction"""
        n = len(pvalues)
        sorted_indices = np.argsort(pvalues)
        sorted_pvalues = pvalues[sorted_indices]
        
        # Calculate Holm correction
        holm_corrections = np.zeros(n)
        for i in range(n):
            holm_corrections[i] = min(1.0, sorted_pvalues[i] * (n - i))
        
        # Ensure monotonicity
        for i in range(1, n):
            holm_corrections[i] = max(holm_corrections[i], holm_corrections[i-1])
        
        # Restore original order
        corrected_pvalues = np.zeros(n)
        corrected_pvalues[sorted_indices] = holm_corrections
        
        return corrected_pvalues
    
    @staticmethod
    def fdr_bh(pvalues: np.ndarray) -> np.ndarray:
        """Benjamini-Hochberg FDR correction"""
        if HAS_STATSMODELS:
            try:
                _, corrected, _, _ = multipletests(pvalues, method='fdr_bh')
                return corrected
            except Exception:
                pass
        
        # Manual implementation
        n = len(pvalues)
        sorted_indices = np.argsort(pvalues)
        sorted_pvalues = pvalues[sorted_indices]
        
        # Calculate BH correction
        bh_corrections = np.zeros(n)
        for i in range(n):
            bh_corrections[i] = min(1.0, sorted_pvalues[i] * n / (i + 1))
        
        # Ensure monotonicity (reverse)
        for i in range(n-2, -1, -1):
            bh_corrections[i] = min(bh_corrections[i], bh_corrections[i+1])
        
        # Restore original order
        corrected_pvalues = np.zeros(n)
        corrected_pvalues[sorted_indices] = bh_corrections
        
        return corrected_pvalues
    
    @staticmethod
    def fdr_by(pvalues: np.ndarray) -> np.ndarray:
        """Benjamini-Yekutieli FDR correction"""
        if HAS_STATSMODELS:
            try:
                _, corrected, _, _ = multipletests(pvalues, method='fdr_by')
                return corrected
            except Exception:
                pass
        
        # Manual implementation
        n = len(pvalues)
        c_n = np.sum(1.0 / np.arange(1, n + 1))  # Harmonic number
        
        # Apply BH correction with c(n) factor
        bh_corrected = MultipleTestingCorrection.fdr_bh(pvalues)
        return np.minimum(bh_corrected * c_n, 1.0)


class PowerAnalysis:
    """Statistical power analysis"""
    
    @staticmethod
    def compute_power(effect_size: float, n_obs: int, n_bg: int,
                     alpha: float = 0.05) -> float:
        """Compute statistical power"""
        if not HAS_SCIPY:
            return PowerAnalysis._approximate_power(effect_size, n_obs, n_bg, alpha)
        
        try:
            # Use normal approximation for power calculation
            # Standard error of the difference
            se_diff = np.sqrt(1/n_obs + 1/n_bg)
            
            # Critical value for two-tailed test
            z_crit = stats.norm.ppf(1 - alpha/2)
            
            # Non-centrality parameter
            ncp = effect_size / se_diff
            
            # Power = P(|Z| > z_crit | H1 is true)
            power = 1 - stats.norm.cdf(z_crit - ncp) + stats.norm.cdf(-z_crit - ncp)
            
            return power
        except Exception:
            return PowerAnalysis._approximate_power(effect_size, n_obs, n_bg, alpha)
    
    @staticmethod
    def _approximate_power(effect_size: float, n_obs: int, n_bg: int,
                          alpha: float) -> float:
        """Approximate power calculation"""
        # Simple approximation based on effect size and sample size
        total_n = n_obs + n_bg
        
        if total_n < 10:
            return 0.1
        elif abs(effect_size) < 0.2:
            return 0.1 + 0.4 * min(total_n / 100, 1)
        elif abs(effect_size) < 0.5:
            return 0.3 + 0.5 * min(total_n / 50, 1)
        else:
            return 0.6 + 0.4 * min(total_n / 30, 1)


class RobustStatisticalFramework:
    """Comprehensive statistical framework for MEDALT"""
    
    def __init__(self, 
                 test_method: str = 'enrichment',
                 correction_method: str = 'fdr_bh',
                 confidence_level: float = 0.95,
                 n_bootstrap: int = 1000,
                 min_observations: int = 5,
                 n_jobs: int = 1):
        
        self.test_method = test_method
        self.correction_method = correction_method
        self.confidence_level = confidence_level
        self.n_bootstrap = n_bootstrap
        self.min_observations = min_observations
        self.n_jobs = n_jobs
        
        # Initialize test object
        if test_method == 'enrichment':
            self.test = EnrichmentTest()
        elif test_method == 'wilcoxon':
            self.test = WilcoxonTest()
        else:
            raise ValueError(f"Unknown test method: {test_method}")
        
        # Initialize bootstrap
        self.bootstrap_ci = BootstrapConfidenceInterval(
            n_bootstrap=n_bootstrap,
            confidence_level=confidence_level
        )
        
        logger.info(f"Initialized RobustStatisticalFramework: {test_method} test, "
                   f"{correction_method} correction, {confidence_level} CI")
    
    def run_statistical_analysis(self, 
                                observed_data: Dict[str, np.ndarray],
                                null_distributions: Dict[str, np.ndarray],
                                background_data: Dict[str, np.ndarray]) -> pd.DataFrame:
        """
        Run comprehensive statistical analysis
        
        Args:
            observed_data: Dict of observed values for each test
            null_distributions: Dict of null distributions for each test
            background_data: Dict of background values for each test
            
        Returns:
            DataFrame with comprehensive statistical results
        """
        logger.info("Running comprehensive statistical analysis...")
        
        # Run tests in parallel if requested
        if self.n_jobs > 1:
            return self._run_parallel_analysis(observed_data, null_distributions, background_data)
        else:
            return self._run_sequential_analysis(observed_data, null_distributions, background_data)
    
    def _run_sequential_analysis(self, observed_data: Dict, null_distributions: Dict,
                               background_data: Dict) -> pd.DataFrame:
        """Run analysis sequentially"""
        results = []
        
        for test_id in observed_data.keys():
            if test_id not in null_distributions or test_id not in background_data:
                continue
                
            obs_data = observed_data[test_id]
            null_dist = null_distributions[test_id]
            bg_data = background_data[test_id]
            
            # Skip if insufficient data
            if len(obs_data) < self.min_observations or len(bg_data) < self.min_observations:
                continue
            
            # Compute statistical result
            stat_result = self._compute_single_test(test_id, obs_data, null_dist, bg_data)
            results.append(stat_result)
        
        return self._format_results(results)
    
    def _run_parallel_analysis(self, observed_data: Dict, null_distributions: Dict,
                             background_data: Dict) -> pd.DataFrame:
        """Run analysis in parallel"""
        # Prepare arguments for parallel processing
        args_list = []
        for test_id in observed_data.keys():
            if test_id not in null_distributions or test_id not in background_data:
                continue
                
            obs_data = observed_data[test_id]
            null_dist = null_distributions[test_id]
            bg_data = background_data[test_id]
            
            if len(obs_data) >= self.min_observations and len(bg_data) >= self.min_observations:
                args_list.append((test_id, obs_data, null_dist, bg_data))
        
        # Run in parallel
        with ProcessPoolExecutor(max_workers=self.n_jobs) as executor:
            results = list(executor.map(self._compute_single_test_wrapper, args_list))
        
        # Filter out None results
        results = [r for r in results if r is not None]
        
        return self._format_results(results)
    
    def _compute_single_test_wrapper(self, args: Tuple) -> Optional[StatisticalResult]:
        """Wrapper for parallel processing"""
        test_id, obs_data, null_dist, bg_data = args
        return self._compute_single_test(test_id, obs_data, null_dist, bg_data)
    
    def _compute_single_test(self, test_id: str, obs_data: np.ndarray,
                           null_dist: np.ndarray, bg_data: np.ndarray) -> StatisticalResult:
        """Compute statistical result for a single test"""
        
        # Compute test statistic
        observed_statistic = self.test.compute_statistic(obs_data, bg_data)
        
        # Compute p-value
        pvalue = self.test.compute_pvalue(observed_statistic, null_dist)
        
        # Compute effect size
        effect_size = self.test.compute_effect_size(obs_data, bg_data)
        
        # Compute confidence interval
        try:
            ci_lower, ci_upper = self.bootstrap_ci.compute_ci(
                obs_data, bg_data, self.test.compute_statistic
            )
        except Exception:
            ci_lower, ci_upper = np.nan, np.nan
        
        # Compute power
        power = PowerAnalysis.compute_power(
            effect_size, len(obs_data), len(bg_data)
        )
        
        return StatisticalResult(
            statistic=observed_statistic,
            pvalue=pvalue,
            effect_size=effect_size,
            confidence_interval=(ci_lower, ci_upper),
            method=self.test_method,
            n_observations=len(obs_data),
            n_permutations=len(null_dist),
            power=power
        )
    
    def _format_results(self, results: List[StatisticalResult]) -> pd.DataFrame:
        """Format results into DataFrame with multiple testing correction"""
        if not results:
            return pd.DataFrame()
        
        # Convert to DataFrame
        data = []
        for i, result in enumerate(results):
            data.append({
                'test_id': i,
                'statistic': result.statistic,
                'pvalue': result.pvalue,
                'effect_size': result.effect_size,
                'ci_lower': result.confidence_interval[0],
                'ci_upper': result.confidence_interval[1],
                'method': result.method,
                'n_observations': result.n_observations,
                'n_permutations': result.n_permutations,
                'power': result.power
            })
        
        df = pd.DataFrame(data)
        
        # Apply multiple testing correction
        if len(df) > 1:
            pvalues = df['pvalue'].values
            
            if self.correction_method == 'bonferroni':
                adjusted_pvalues = MultipleTestingCorrection.bonferroni(pvalues)
            elif self.correction_method == 'holm':
                adjusted_pvalues = MultipleTestingCorrection.holm(pvalues)
            elif self.correction_method == 'fdr_bh':
                adjusted_pvalues = MultipleTestingCorrection.fdr_bh(pvalues)
            elif self.correction_method == 'fdr_by':
                adjusted_pvalues = MultipleTestingCorrection.fdr_by(pvalues)
            else:
                adjusted_pvalues = pvalues
            
            df['adjusted_pvalue'] = adjusted_pvalues
            
            # Add significance indicators
            df['significant'] = df['adjusted_pvalue'] < 0.05
            df['highly_significant'] = df['adjusted_pvalue'] < 0.01
            
            # Add effect size interpretation
            df['effect_size_magnitude'] = pd.cut(
                np.abs(df['effect_size']),
                bins=[0, 0.2, 0.5, 0.8, np.inf],
                labels=['negligible', 'small', 'medium', 'large']
            )
        
        return df.sort_values('adjusted_pvalue' if 'adjusted_pvalue' in df.columns else 'pvalue')


def main():
    """Example usage and testing"""
    
    # Create synthetic data
    np.random.seed(42)
    
    # Simulate observed data (enriched)
    observed_data = {
        'test1': np.random.normal(3.0, 1.0, 50),  # Enriched
        'test2': np.random.normal(2.0, 1.0, 30),  # Neutral
        'test3': np.random.normal(1.0, 1.0, 40),  # Depleted
    }
    
    # Simulate background data
    background_data = {
        'test1': np.random.normal(2.0, 1.0, 100),
        'test2': np.random.normal(2.0, 1.0, 100),
        'test3': np.random.normal(2.0, 1.0, 100),
    }
    
    # Simulate null distributions
    null_distributions = {
        'test1': np.random.normal(0, 1, 1000),
        'test2': np.random.normal(0, 1, 1000),
        'test3': np.random.normal(0, 1, 1000),
    }
    
    # Run statistical analysis
    framework = RobustStatisticalFramework(
        test_method='enrichment',
        correction_method='fdr_bh',
        confidence_level=0.95,
        n_bootstrap=1000
    )
    
    results = framework.run_statistical_analysis(
        observed_data, null_distributions, background_data
    )
    
    print("Statistical Analysis Results:")
    print("=" * 50)
    print(results)
    
    # Show interpretation
    if not results.empty:
        print("\nInterpretation:")
        print("-" * 30)
        significant_tests = results[results['significant']]
        print(f"Significant tests: {len(significant_tests)}")
        
        for idx, row in significant_tests.iterrows():
            print(f"Test {row['test_id']}: "
                  f"effect size = {row['effect_size']:.3f} ({row['effect_size_magnitude']}), "
                  f"p = {row['pvalue']:.3e}, "
                  f"adj. p = {row['adjusted_pvalue']:.3e}")


if __name__ == '__main__':
    main()