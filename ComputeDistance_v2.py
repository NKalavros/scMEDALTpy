#!/usr/bin/env python2
"""
Version 2: Ultra-optimized distance calculation for massive datasets
"""

import numpy as np
import time
from collections import deque

def matrixbuilder_v2(node):
    """Ultra-optimized matrix builder for large datasets"""
    print "Building distance matrix with v2 optimizations..."
    
    node_keys = list(node.keys())
    n = len(node_keys)
    
    print "Processing %d cells..." % n
    
    if n > 1000:
        print "WARNING: Large dataset detected. Using memory-efficient processing."
        return matrixbuilder_chunked(node)
    
    # Pre-allocate matrix
    matrix = np.zeros((n, n), dtype=np.float32)  # Use float32 to save memory
    
    # Convert all nodes to numpy arrays once
    print "Converting to numpy arrays..."
    node_arrays = {}
    for i, key in enumerate(node_keys):
        # Flatten all chromosomal segments into one array
        flattened = []
        for segment in node[key]:
            flattened.extend(segment)
        node_arrays[key] = np.array(flattened, dtype=np.int16)  # Use int16 to save memory
        
        if i % 100 == 0:
            print "  Converted %d/%d cells" % (i, n)
    
    print "Calculating distances..."
    
    # Calculate distances with progress reporting
    total_calcs = n * (n - 1) / 2
    calc_count = 0
    
    for i, key1 in enumerate(node_keys):
        for j, key2 in enumerate(node_keys):
            if i == j:
                matrix[i, j] = 0.0
            elif i < j:  # Only calculate upper triangle
                dist_val = dist_numpy_optimized_v2(node_arrays[key1], node_arrays[key2])
                matrix[i, j] = dist_val
                matrix[j, i] = dist_val  # Symmetric
                
                calc_count += 1
                if calc_count % 10000 == 0:
                    print "  Progress: %d/%d (%.1f%%)" % (calc_count, total_calcs, 100.0 * calc_count / total_calcs)
    
    return node_keys, matrix.tolist()

def matrixbuilder_chunked(node):
    """Process very large datasets in chunks to avoid memory issues"""
    node_keys = list(node.keys())
    n = len(node_keys)
    
    print "Using chunked processing for %d cells" % n
    
    # Use smaller chunks for very large datasets
    chunk_size = min(500, n / 4)
    
    # Pre-allocate matrix on disk if too large for memory
    if n > 2000:
        print "Dataset too large for memory - using approximate methods"
        return matrixbuilder_approximate(node)
    
    matrix = np.zeros((n, n), dtype=np.float32)
    
    # Process in chunks
    for chunk_start in range(0, n, chunk_size):
        chunk_end = min(chunk_start + chunk_size, n)
        print "Processing chunk %d-%d..." % (chunk_start, chunk_end)
        
        # Convert chunk to numpy
        chunk_arrays = {}
        for i in range(chunk_start, chunk_end):
            key = node_keys[i]
            flattened = []
            for segment in node[key]:
                flattened.extend(segment)
            chunk_arrays[key] = np.array(flattened, dtype=np.int16)
        
        # Calculate distances for this chunk
        for i in range(chunk_start, chunk_end):
            key1 = node_keys[i]
            for j in range(i, n):
                key2 = node_keys[j]
                
                if i == j:
                    matrix[i, j] = 0.0
                else:
                    # Load second array if not in chunk
                    if key2 not in chunk_arrays:
                        flattened = []
                        for segment in node[key2]:
                            flattened.extend(segment)
                        arr2 = np.array(flattened, dtype=np.int16)
                    else:
                        arr2 = chunk_arrays[key2]
                    
                    dist_val = dist_numpy_optimized_v2(chunk_arrays[key1], arr2)
                    matrix[i, j] = dist_val
                    matrix[j, i] = dist_val
    
    return node_keys, matrix.tolist()

def matrixbuilder_approximate(node):
    """Use approximate methods for extremely large datasets"""
    node_keys = list(node.keys())
    n = len(node_keys)
    
    print "Using approximate distance calculation for %d cells" % n
    print "This will be faster but less precise - suitable for large-scale analysis"
    
    # Sample-based approach for very large datasets
    max_samples = 1000
    if n > max_samples:
        print "Sampling %d cells from %d total" % (max_samples, n)
        # Take every n/max_samples cells
        step = n / max_samples
        sampled_keys = [node_keys[i] for i in range(0, n, step)][:max_samples]
    else:
        sampled_keys = node_keys
    
    print "Computing on %d sampled cells..." % len(sampled_keys)
    
    # Compute distance matrix for sample
    sample_matrix = np.zeros((len(sampled_keys), len(sampled_keys)), dtype=np.float32)
    
    # Convert sample to numpy
    sample_arrays = {}
    for key in sampled_keys:
        flattened = []
        for segment in node[key]:
            flattened.extend(segment)
        sample_arrays[key] = np.array(flattened, dtype=np.int16)
    
    # Calculate sample distances
    for i, key1 in enumerate(sampled_keys):
        for j, key2 in enumerate(sampled_keys):
            if i == j:
                sample_matrix[i, j] = 0.0
            elif i < j:
                dist_val = dist_numpy_optimized_v2(sample_arrays[key1], sample_arrays[key2])
                sample_matrix[i, j] = dist_val
                sample_matrix[j, i] = dist_val
    
    # Interpolate for full matrix
    full_matrix = np.zeros((n, n), dtype=np.float32)
    
    print "Interpolating distances for full dataset..."
    for i in range(n):
        for j in range(i, n):
            if i == j:
                full_matrix[i, j] = 0.0
            else:
                # Find nearest samples
                sample_i = min(i * len(sampled_keys) / n, len(sampled_keys) - 1)
                sample_j = min(j * len(sampled_keys) / n, len(sampled_keys) - 1)
                
                dist_val = sample_matrix[sample_i, sample_j]
                full_matrix[i, j] = dist_val
                full_matrix[j, i] = dist_val
    
    return node_keys, full_matrix.tolist()

def dist_numpy_optimized_v2(arr1, arr2):
    """Ultra-fast numpy distance calculation"""
    # Handle different lengths
    min_len = min(len(arr1), len(arr2))
    if min_len == 0:
        return 1000000.0
    
    a1 = arr1[:min_len]
    a2 = arr2[:min_len]
    
    # Check for zero incompatibilities quickly
    zero_mask1 = (a1 == 0)
    zero_mask2 = (a2 == 0)
    
    if np.any(zero_mask1 != zero_mask2):
        return 1000000.0
    
    # Filter out zeros
    valid_mask = ~(zero_mask1 | zero_mask2)
    if not np.any(valid_mask):
        return 0.0
    
    filtered1 = a1[valid_mask]
    filtered2 = a2[valid_mask]
    
    # Fast MED calculation
    return med_distance_fast(filtered1, filtered2)

def med_distance_fast(arr1, arr2):
    """Optimized MED distance using numpy operations where possible"""
    if len(arr1) != len(arr2):
        return 1000000.0
    
    if len(arr1) == 0:
        return 0.0
    
    if len(arr1) == 1:
        return float(abs(arr1[0] - arr2[0]))
    
    # For very short segments, use original algorithm
    if len(arr1) <= 10:
        return med_distance_simple(arr1.tolist(), arr2.tolist())
    
    # For longer segments, use approximation
    diff = arr1.astype(np.float32) - arr2.astype(np.float32)
    
    # Simple approximation: sum of absolute differences
    # This is faster but less accurate than true MED
    return float(np.sum(np.abs(diff)))

def med_distance_simple(list1, list2):
    """Simple MED calculation for short lists"""
    if len(list1) != len(list2):
        return 1000000.0
    
    if len(list1) == 1:
        return abs(list1[0] - list2[0])
    
    diff_list = deque([list1[i] - list2[i] for i in range(len(list1))])
    
    d = 0
    while diff_list:
        if diff_list[0] == 0:
            diff_list.popleft()
        elif diff_list[0] > 0:
            k = 0
            for i in range(len(diff_list)):
                if diff_list[i] > 0:
                    k = i
                else:
                    break
            for i in range(k + 1):
                diff_list[i] -= 1
            d += 1
        else:
            k = 0
            for i in range(len(diff_list)):
                if diff_list[i] < 0:
                    k = i
                else:
                    break
            for i in range(k + 1):
                diff_list[i] += 1
            d += 1
    
    return d

# Wrapper functions for compatibility
def matrixbuilder(node):
    """Main interface - choose best method based on dataset size"""
    n = len(node)
    
    if n <= 100:
        # Use original optimized version for small datasets
        from ComputeDistance import matrixbuilder as matrixbuilder_orig
        return matrixbuilder_orig(node)
    else:
        # Use v2 optimizations for larger datasets
        return matrixbuilder_v2(node)

def dist(node1, node2):
    """Main distance function interface"""
    # Convert to numpy and use optimized calculation
    flattened1 = []
    flattened2 = []
    
    for segment in node1:
        flattened1.extend(segment)
    for segment in node2:
        flattened2.extend(segment)
    
    arr1 = np.array(flattened1, dtype=np.int16)
    arr2 = np.array(flattened2, dtype=np.int16)
    
    return dist_numpy_optimized_v2(arr1, arr2)