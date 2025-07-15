#!/usr/bin/env python2
"""
Parallel processing version of distance calculation for large datasets
"""

import numpy as np
import time
import multiprocessing as mp
from collections import deque

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

def calculate_distance_chunk(args):
    """Calculate distances for a chunk of cell pairs - for parallel processing"""
    node_data, node_keys, row_start, row_end, col_start, col_end = args
    
    chunk_results = []
    
    for i in range(row_start, row_end):
        key1 = node_keys[i]
        for j in range(max(col_start, i), col_end):
            key2 = node_keys[j]
            
            if i == j:
                dist_val = 0.0
            else:
                dist_val = dist_accurate_optimized(node_data[key1], node_data[key2])
            
            chunk_results.append((i, j, dist_val))
    
    return chunk_results

def matrixbuilder_parallel(node):
    """Parallel matrix builder using all available cores"""
    print "Building distance matrix with parallel processing..."
    
    node_keys = list(node.keys())
    n = len(node_keys)
    
    print "Processing %d cells with %d cores..." % (n, mp.cpu_count())
    
    # Pre-allocate matrix
    matrix = np.zeros((n, n), dtype=np.float32)
    
    # Preprocess all nodes
    print "Preprocessing nodes..."
    node_data = {}
    for i, key in enumerate(node_keys):
        node_data[key] = preprocess_node(node[key])
        if i % 100 == 0 and i > 0:
            print "  Preprocessed %d/%d cells" % (i, n)
    
    print "Setting up parallel processing..."
    
    # Create work chunks for parallel processing
    num_cores = mp.cpu_count()
    
    # Divide work by rows to balance load
    chunk_size = max(1, n / num_cores)
    chunks = []
    
    for core_id in range(num_cores):
        row_start = core_id * chunk_size
        row_end = min((core_id + 1) * chunk_size, n) if core_id < num_cores - 1 else n
        
        if row_start < n:
            chunks.append((node_data, node_keys, row_start, row_end, 0, n))
    
    print "Created %d work chunks" % len(chunks)
    
    # Process chunks in parallel
    print "Starting parallel distance calculations..."
    start_time = time.time()
    
    pool = mp.Pool(processes=num_cores)
    
    try:
        results = pool.map(calculate_distance_chunk, chunks)
        
        # Collect results
        total_distances = 0
        for chunk_results in results:
            for i, j, dist_val in chunk_results:
                matrix[i, j] = dist_val
                if i != j:
                    matrix[j, i] = dist_val  # Symmetric
                total_distances += 1
        
        elapsed = time.time() - start_time
        print "Parallel processing completed in %.1f seconds" % elapsed
        print "Calculated %d distances" % total_distances
        
    finally:
        pool.close()
        pool.join()
    
    print "Distance matrix completed!"
    return node_keys, matrix.tolist()

def matrixbuilder_parallel_smart(node):
    """Smart parallel processing that only calculates upper triangle"""
    print "Building distance matrix with smart parallel processing..."
    
    node_keys = list(node.keys())
    n = len(node_keys)
    
    print "Processing %d cells with %d cores..." % (n, mp.cpu_count())
    
    # Pre-allocate matrix
    matrix = np.zeros((n, n), dtype=np.float32)
    
    # Preprocess all nodes
    print "Preprocessing nodes..."
    node_data = {}
    for i, key in enumerate(node_keys):
        node_data[key] = preprocess_node(node[key])
        if i % 100 == 0 and i > 0:
            print "  Preprocessed %d/%d cells" % (i, n)
    
    print "Setting up smart parallel processing..."
    
    # Create work chunks for upper triangle only
    num_cores = mp.cpu_count()
    total_calcs = n * (n - 1) / 2
    calcs_per_core = total_calcs / num_cores
    
    chunks = []
    current_calc = 0
    
    for core_id in range(num_cores):
        target_calcs = (core_id + 1) * calcs_per_core
        chunk_pairs = []
        
        # Find pairs for this chunk
        for i in range(n):
            for j in range(i + 1, n):
                if current_calc >= target_calcs and core_id < num_cores - 1:
                    break
                chunk_pairs.append((i, j))
                current_calc += 1
            if current_calc >= target_calcs and core_id < num_cores - 1:
                break
        
        if chunk_pairs:
            chunks.append((node_data, node_keys, chunk_pairs))
    
    print "Created %d work chunks with %d total calculations" % (len(chunks), total_calcs)
    
    # Process chunks in parallel
    print "Starting smart parallel distance calculations..."
    start_time = time.time()
    
    pool = mp.Pool(processes=num_cores)
    
    try:
        results = pool.map(calculate_distance_pairs, chunks)
        
        # Collect results
        total_distances = 0
        for chunk_results in results:
            for i, j, dist_val in chunk_results:
                matrix[i, j] = dist_val
                matrix[j, i] = dist_val  # Symmetric
                total_distances += 1
        
        elapsed = time.time() - start_time
        print "Smart parallel processing completed in %.1f seconds" % elapsed
        print "Calculated %d unique distances" % total_distances
        
    finally:
        pool.close()
        pool.join()
    
    print "Distance matrix completed!"
    return node_keys, matrix.tolist()

def calculate_distance_pairs(args):
    """Calculate distances for specific pairs - for smart parallel processing"""
    node_data, node_keys, pairs = args
    
    chunk_results = []
    
    for i, j in pairs:
        key1 = node_keys[i]
        key2 = node_keys[j]
        
        dist_val = dist_accurate_optimized(node_data[key1], node_data[key2])
        chunk_results.append((i, j, dist_val))
    
    return chunk_results

# Main interface
def matrixbuilder(node):
    """Main interface - use parallel processing for large datasets"""
    n = len(node)
    
    print "Selecting processing method for %d cells..." % n
    
    if n <= 100:
        # Use simple version for small datasets
        from ComputeDistance import matrixbuilder as matrixbuilder_simple
        return matrixbuilder_simple(node)
    else:
        # Use parallel processing for larger datasets
        return matrixbuilder_parallel_smart(node)

def dist(node1, node2):
    """Main distance function interface"""
    node1_data = preprocess_node(node1)
    node2_data = preprocess_node(node2)
    return dist_accurate_optimized(node1_data, node2_data)