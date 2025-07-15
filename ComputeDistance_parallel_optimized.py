#!/usr/bin/env python2
"""
Optimized parallel processing version with reduced overhead
"""

import numpy as np
import time
import multiprocessing as mp

def preprocess_node(node_segments):
    """Preprocess a node for faster distance calculation"""
    all_values = []
    segment_boundaries = [0]
    
    for segment in node_segments:
        all_values.extend(segment)
        segment_boundaries.append(len(all_values))
    
    return {
        'values': np.array(all_values, dtype=np.int16),
        'boundaries': segment_boundaries[:-1]
    }

def med_distance_core_optimized(arr1, arr2):
    """Core MED distance with optimizations but maintaining accuracy"""
    if len(arr1) != len(arr2):
        return 1000000.0
    
    if len(arr1) == 0:
        return 0.0
    
    diff_list = (arr1 - arr2).tolist()
    
    d = 0
    i = 0
    
    while i < len(diff_list):
        if diff_list[i] == 0:
            i += 1
        elif diff_list[i] > 0:
            k = i
            while k < len(diff_list) and diff_list[k] > 0:
                k += 1
            
            for j in range(i, k):
                diff_list[j] -= 1
            d += 1
        else:
            k = i
            while k < len(diff_list) and diff_list[k] < 0:
                k += 1
            
            for j in range(i, k):
                diff_list[j] += 1
            d += 1
    
    return float(d)

def calculate_distance_batch_optimized(args):
    """Optimized batch distance calculation with minimal data transfer"""
    node_data_subset, pairs = args
    
    results = []
    
    for i, j in pairs:
        # Get preprocessed data
        node1_data = node_data_subset[i]
        node2_data = node_data_subset[j]
        
        arr1 = node1_data['values']
        arr2 = node2_data['values']
        boundaries1 = node1_data['boundaries']
        boundaries2 = node2_data['boundaries']
        
        total_dist = 0.0
        
        # Process each segment
        for seg_idx in range(len(boundaries1)):
            start1 = boundaries1[seg_idx]
            end1 = boundaries1[seg_idx + 1] if seg_idx + 1 < len(boundaries1) else len(arr1)
            
            start2 = boundaries2[seg_idx]
            end2 = boundaries2[seg_idx + 1] if seg_idx + 1 < len(boundaries2) else len(arr2)
            
            seg1 = arr1[start1:end1]
            seg2 = arr2[start2:end2]
            
            # Handle different segment lengths
            if len(seg1) != len(seg2):
                max_len = max(len(seg1), len(seg2))
                if len(seg1) < max_len:
                    seg1 = np.pad(seg1, (0, max_len - len(seg1)), 'constant', constant_values=0)
                if len(seg2) < max_len:
                    seg2 = np.pad(seg2, (0, max_len - len(seg2)), 'constant', constant_values=0)
            
            # Check for zeros
            zero_mask1 = (seg1 == 0)
            zero_mask2 = (seg2 == 0)
            
            if np.any(zero_mask1 != zero_mask2):
                total_dist = 1000000.0
                break
            
            # Remove zero positions
            valid_mask = ~(zero_mask1 | zero_mask2)
            if not np.any(valid_mask):
                continue
            
            filtered1 = seg1[valid_mask]
            filtered2 = seg2[valid_mask]
            
            if len(filtered1) == 1:
                total_dist += float(abs(filtered1[0] - filtered2[0]))
            else:
                total_dist += med_distance_core_optimized(filtered1, filtered2)
        
        results.append((i, j, total_dist))
    
    return results

def matrixbuilder_parallel_optimized(node):
    """Highly optimized parallel matrix builder"""
    print "Building distance matrix with optimized parallel processing..."
    
    node_keys = list(node.keys())
    n = len(node_keys)
    
    num_cores = mp.cpu_count()
    print "Processing %d cells with %d cores..." % (n, num_cores)
    
    # Pre-allocate matrix
    matrix = np.zeros((n, n), dtype=np.float32)
    
    # Preprocess all nodes
    print "Preprocessing nodes..."
    node_data = {}
    for i, key in enumerate(node_keys):
        node_data[key] = preprocess_node(node[key])
        if i % 200 == 0 and i > 0:
            print "  Preprocessed %d/%d cells" % (i, n)
    
    # Create index-based node data for efficient transfer
    node_data_indexed = {}
    for i, key in enumerate(node_keys):
        node_data_indexed[i] = node_data[key]
    
    # Generate all pairs for upper triangle
    print "Generating work pairs..."
    total_pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            total_pairs.append((i, j))
    
    total_calcs = len(total_pairs)
    print "Total unique calculations: %d" % total_calcs
    
    # Distribute pairs across cores with optimal batch sizes
    batch_size = max(100, total_calcs / (num_cores * 4))  # Larger batches for efficiency
    batches = []
    
    for i in range(0, total_calcs, batch_size):
        batch_pairs = total_pairs[i:i+batch_size]
        batches.append((node_data_indexed, batch_pairs))
    
    print "Created %d batches with ~%d pairs each" % (len(batches), batch_size)
    
    # Process batches in parallel
    print "Starting optimized parallel processing..."
    start_time = time.time()
    
    pool = mp.Pool(processes=num_cores)
    
    try:
        all_results = pool.map(calculate_distance_batch_optimized, batches)
        
        # Collect results
        total_processed = 0
        for batch_results in all_results:
            for i, j, dist_val in batch_results:
                matrix[i, j] = dist_val
                matrix[j, i] = dist_val  # Symmetric
                total_processed += 1
        
        elapsed = time.time() - start_time
        print "Optimized parallel processing completed in %.1f seconds" % elapsed
        print "Processed %d distance calculations" % total_processed
        print "Speed: %.1f calculations/second" % (total_processed / elapsed)
        
    finally:
        pool.close()
        pool.join()
    
    print "Distance matrix completed!"
    return node_keys, matrix.tolist()

def matrixbuilder_hybrid(node):
    """Hybrid approach: use parallel for large datasets, sequential for small"""
    n = len(node)
    
    print "Selecting processing method for %d cells..." % n
    
    if n <= 200:
        # Use sequential for small datasets (less overhead)
        from ComputeDistance import matrixbuilder as matrixbuilder_seq
        return matrixbuilder_seq(node)
    else:
        # Use optimized parallel for larger datasets
        return matrixbuilder_parallel_optimized(node)

# Main interface
def matrixbuilder(node):
    """Main interface with hybrid approach"""
    return matrixbuilder_hybrid(node)

def dist(node1, node2):
    """Main distance function interface"""
    node1_data = preprocess_node(node1)
    node2_data = preprocess_node(node2)
    
    # Use the same distance calculation as the parallel version
    arr1 = node1_data['values']
    arr2 = node2_data['values']
    boundaries1 = node1_data['boundaries']
    boundaries2 = node2_data['boundaries']
    
    total_dist = 0.0
    
    for i in range(len(boundaries1)):
        start1 = boundaries1[i]
        end1 = boundaries1[i + 1] if i + 1 < len(boundaries1) else len(arr1)
        
        start2 = boundaries2[i]
        end2 = boundaries2[i + 1] if i + 1 < len(boundaries2) else len(arr2)
        
        seg1 = arr1[start1:end1]
        seg2 = arr2[start2:end2]
        
        if len(seg1) != len(seg2):
            max_len = max(len(seg1), len(seg2))
            if len(seg1) < max_len:
                seg1 = np.pad(seg1, (0, max_len - len(seg1)), 'constant', constant_values=0)
            if len(seg2) < max_len:
                seg2 = np.pad(seg2, (0, max_len - len(seg2)), 'constant', constant_values=0)
        
        zero_mask1 = (seg1 == 0)
        zero_mask2 = (seg2 == 0)
        
        if np.any(zero_mask1 != zero_mask2):
            return 1000000.0
        
        valid_mask = ~(zero_mask1 | zero_mask2)
        if not np.any(valid_mask):
            continue
        
        filtered1 = seg1[valid_mask]
        filtered2 = seg2[valid_mask]
        
        if len(filtered1) == 1:
            total_dist += float(abs(filtered1[0] - filtered2[0]))
        else:
            total_dist += med_distance_core_optimized(filtered1, filtered2)
    
    return total_dist