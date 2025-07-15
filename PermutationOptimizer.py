#!/usr/bin/env python2
"""
Optimization utilities for MEDALT permutation analysis
"""

import os
import sys
from collections import defaultdict

class PermutationOptimizer(object):
    """
    Optimize permutation analysis with early termination and adaptive sampling
    """
    
    def __init__(self, min_permutations=20, max_permutations=100, significance_threshold=0.05):
        self.min_permutations = min_permutations
        self.max_permutations = max_permutations
        self.significance_threshold = significance_threshold
        self.permutation_results = []
        
    def should_continue_permutation(self, current_permutation, test_statistic):
        """
        Determine if we should continue with more permutations based on 
        current results and statistical power
        """
        self.permutation_results.append(test_statistic)
        
        # Always do minimum number of permutations
        if current_permutation < self.min_permutations:
            return True
            
        # If we've done max permutations, stop
        if current_permutation >= self.max_permutations:
            return False
            
        # Calculate current p-value estimate
        num_extreme = sum(1 for x in self.permutation_results if x >= test_statistic)
        current_p_value = float(num_extreme) / len(self.permutation_results)
        
        # Early termination conditions
        if current_p_value > 0.1:  # Very non-significant
            return False
        elif current_p_value < 0.001:  # Very significant
            return False
        elif current_permutation > 50 and current_p_value > 0.05:  # Moderately non-significant
            return False
            
        return True
        
    def estimate_required_permutations(self, dataset_size):
        """
        Estimate the number of permutations needed based on dataset characteristics
        """
        num_cells, num_genes = dataset_size
        
        # For larger datasets, we can use fewer permutations due to better statistical power
        if num_cells > 5000 or num_genes > 2000:
            return 30
        elif num_cells > 1000 or num_genes > 1000:
            return 50
        else:
            return 100