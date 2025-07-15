#!/usr/bin/env python2
"""
Version 2 Accurate: Ultra-optimized distance calculation maintaining 100% accuracy
"""

import numpy as np
import time
from collections import deque

def matrixbuilder_v2_accurate(node):
    """Ultra-optimized but accurate matrix builder"""
    print "Building distance matrix with v2 accurate optimizations..."
    
    node_keys = list(node.keys())
    n = len(node_keys)
    
    print "Processing %d cells..." % n
    
    # Pre-allocate matrix with float32 for memory efficiency
    matrix = np.zeros((n, n), dtype=np.float32)
    
    # Convert nodes to more efficient representation
    print "Preprocessing nodes..."
    node_data = {}
    for i, key in enumerate(node_keys):
        node_data[key] = preprocess_node(node[key])
        if i % 100 == 0 and i > 0:
            print "  Preprocessed %d/%d cells" % (i, n)
    
    print "Calculating distances..."
    
    # Calculate distances with vectorized operations where possible
    total_calcs = n * (n - 1) / 2
    calc_count = 0
    
    for i, key1 in enumerate(node_keys):
        for j, key2 in enumerate(node_keys):
            if i == j:
                matrix[i, j] = 0.0
            elif i < j:  # Only calculate upper triangle
                dist_val = dist_accurate_optimized(node_data[key1], node_data[key2])
                matrix[i, j] = dist_val
                matrix[j, i] = dist_val  # Symmetric
                
                calc_count += 1
                if calc_count % 5000 == 0:
                    progress = 100.0 * calc_count / total_calcs
                    print "  Progress: %d/%d (%.1f%%)" % (calc_count, total_calcs, progress)
    
    print "Distance matrix completed!"
    return node_keys, matrix.tolist()

def preprocess_node(node_segments):
    """Preprocess a node for faster distance calculation"""
    # Flatten all segments but keep track of segment boundaries
    all_values = []
    segment_boundaries = [0]
    
    for segment in node_segments:
        all_values.extend(segment)
        segment_boundaries.append(len(all_values))
    
    return {
        'values': np.array(all_values, dtype=np.int16),
        'boundaries': segment_boundaries[:-1]  # Don't need the last boundary
    }

def dist_accurate_optimized(node1_data, node2_data):
    """Optimized but accurate distance calculation"""
    arr1 = node1_data['values']
    arr2 = node2_data['values']
    boundaries1 = node1_data['boundaries']
    boundaries2 = node2_data['boundaries']
    
    # Ensure same number of segments
    if len(boundaries1) != len(boundaries2):
        raise ValueError("Nodes must have same number of segments")
    
    total_dist = 0.0
    
    # Process each segment
    for i in range(len(boundaries1)):
        # Get segment boundaries
        start1 = boundaries1[i]
        end1 = boundaries1[i + 1] if i + 1 < len(boundaries1) else len(arr1)
        
        start2 = boundaries2[i]
        end2 = boundaries2[i + 1] if i + 1 < len(boundaries2) else len(arr2)
        
        # Extract segments
        seg1 = arr1[start1:end1]
        seg2 = arr2[start2:end2]
        
        # Calculate distance for this segment
        if len(seg1) != len(seg2):
            # Pad shorter segment with zeros
            max_len = max(len(seg1), len(seg2))
            if len(seg1) < max_len:
                seg1 = np.pad(seg1, (0, max_len - len(seg1)), 'constant', constant_values=0)
            if len(seg2) < max_len:
                seg2 = np.pad(seg2, (0, max_len - len(seg2)), 'constant', constant_values=0)
        
        segment_dist = med_distance_accurate(seg1, seg2)
        total_dist += segment_dist
    
    return total_dist

def med_distance_accurate(seg1, seg2):
    """Accurate MED distance calculation with numpy optimizations"""
    # Check for zeros
    zero_mask1 = (seg1 == 0)
    zero_mask2 = (seg2 == 0)
    
    # Handle zero incompatibilities
    if np.any(zero_mask1 != zero_mask2):
        return 1000000.0
    
    # Remove zero positions
    valid_mask = ~(zero_mask1 | zero_mask2)
    if not np.any(valid_mask):
        return 0.0
    
    filtered1 = seg1[valid_mask]
    filtered2 = seg2[valid_mask]
    
    if len(filtered1) == 1:
        return float(abs(filtered1[0] - filtered2[0]))
    
    # Use optimized MED calculation
    return med_distance_core_optimized(filtered1, filtered2)

def med_distance_core_optimized(arr1, arr2):
    """Core MED distance with optimizations but maintaining accuracy"""
    if len(arr1) != len(arr2):
        return 1000000.0
    
    if len(arr1) == 0:
        return 0.0
    
    # Convert to Python lists for the algorithm (numpy slicing is slower for the MED algorithm)
    diff_list = (arr1 - arr2).tolist()
    
    d = 0
    i = 0
    
    while i < len(diff_list):
        if diff_list[i] == 0:
            i += 1
        elif diff_list[i] > 0:
            # Find consecutive positive values
            k = i
            while k < len(diff_list) and diff_list[k] > 0:
                k += 1
            
            # Subtract 1 from all positive values in range [i, k)
            for j in range(i, k):
                diff_list[j] -= 1
            d += 1
            
            # Don't advance i - we might have created zeros
        else:  # diff_list[i] < 0
            # Find consecutive negative values
            k = i
            while k < len(diff_list) and diff_list[k] < 0:
                k += 1
            
            # Add 1 to all negative values in range [i, k)
            for j in range(i, k):
                diff_list[j] += 1
            d += 1
            
            # Don't advance i - we might have created zeros
    
    return float(d)

def matrixbuilder_chunked_accurate(node):
    """Chunked processing for very large datasets while maintaining accuracy"""
    node_keys = list(node.keys())
    n = len(node_keys)
    
    print "Using chunked accurate processing for %d cells" % n
    
    chunk_size = 250  # Process 250 cells at a time
    matrix = np.zeros((n, n), dtype=np.float32)
    
    # Preprocess all nodes
    print "Preprocessing all nodes..."
    node_data = {}
    for i, key in enumerate(node_keys):
        node_data[key] = preprocess_node(node[key])
        if i % 100 == 0 and i > 0:
            print "  Preprocessed %d/%d" % (i, n)
    
    # Process in chunks
    for chunk_start in range(0, n, chunk_size):
        chunk_end = min(chunk_start + chunk_size, n)
        print "Processing chunk %d-%d..." % (chunk_start, chunk_end)
        
        for i in range(chunk_start, chunk_end):
            key1 = node_keys[i]
            for j in range(i, n):
                key2 = node_keys[j]
                
                if i == j:
                    matrix[i, j] = 0.0
                else:
                    dist_val = dist_accurate_optimized(node_data[key1], node_data[key2])
                    matrix[i, j] = dist_val
                    matrix[j, i] = dist_val
    
    return node_keys, matrix.tolist()

# Main interface functions
def matrixbuilder(node):
    """Main interface - choose best method based on dataset size"""
    n = len(node)
    
    print "Selecting optimization level for %d cells..." % n
    
    if n <= 100:
        # Use original optimized version for small datasets
        from ComputeDistance import matrixbuilder as matrixbuilder_orig
        return matrixbuilder_orig(node)
    elif n <= 1000:
        # Use v2 accurate for medium datasets
        return matrixbuilder_v2_accurate(node)
    else:
        # Use chunked processing for very large datasets
        return matrixbuilder_chunked_accurate(node)

def dist(node1, node2):
    """Main distance function interface"""
    node1_data = preprocess_node(node1)
    node2_data = preprocess_node(node2)
    return dist_accurate_optimized(node1_data, node2_data)