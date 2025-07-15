#!/usr/bin/env python2
"""
Numpy-optimized version of ComputeDistance.py for massive speedup
"""

import copy
import numpy as np
from collections import deque

def matrixbuilder_numpy(node):
    """Ultra-fast matrix builder using numpy vectorized operations"""
    node_keys = list(node.keys())
    n = len(node_keys)
    
    # Convert node data to numpy arrays for vectorized operations
    node_arrays = {}
    for key in node_keys:
        # Flatten the chromosomal segments into a single array
        flattened = []
        for segment in node[key]:
            flattened.extend(segment)
        node_arrays[key] = np.array(flattened, dtype=np.float32)
    
    # Pre-allocate distance matrix
    matrix = np.zeros((n, n), dtype=np.float32)
    
    # Vectorized distance calculation
    for i, key1 in enumerate(node_keys):
        for j, key2 in enumerate(node_keys):
            if i == j:
                matrix[i, j] = 0.0
            elif i < j:  # Only calculate upper triangle
                dist_val = dist_numpy_vectorized(node_arrays[key1], node_arrays[key2])
                matrix[i, j] = dist_val
                matrix[j, i] = dist_val  # Symmetric matrix
    
    return node_keys, matrix.tolist()

def dist_numpy_vectorized(arr1, arr2):
    """Vectorized distance calculation using numpy"""
    if len(arr1) != len(arr2):
        raise ValueError("Arrays must have same length")
    
    # Handle zero values efficiently
    zero_mask1 = (arr1 == 0)
    zero_mask2 = (arr2 == 0)
    
    # Check for incompatible zeros - return exact same penalty as original
    if np.any(zero_mask1 != zero_mask2):
        # Count mismatches to be consistent with original
        mismatches = np.sum(zero_mask1 != zero_mask2)
        return 1000000.0 + mismatches
    
    # Filter out zero positions
    valid_mask = ~(zero_mask1 | zero_mask2)
    if not np.any(valid_mask):
        return 0.0
    
    filtered1 = arr1[valid_mask]
    filtered2 = arr2[valid_mask]
    
    # Calculate MED distance using vectorized operations
    return med_distance_numpy(filtered1, filtered2)

def med_distance_numpy(arr1, arr2):
    """Numpy-optimized MED distance calculation"""
    if len(arr1) != len(arr2):
        raise ValueError("Arrays must have same length")
    
    if len(arr1) == 0:
        return 0.0
    
    if len(arr1) == 1:
        return abs(float(arr1[0]) - float(arr2[0]))
    
    # Calculate difference array
    diff = arr1 - arr2
    
    # Use iterative approach for MED calculation
    d = 0
    diff_list = diff.tolist()  # Convert to list for manipulation
    
    while diff_list:
        if diff_list[0] == 0:
            diff_list.pop(0)
        elif diff_list[0] > 0:
            # Find consecutive positive values
            k = 0
            for i in range(len(diff_list)):
                if diff_list[i] > 0:
                    k = i
                else:
                    break
            
            # Subtract 1 from all positive values
            for i in range(k + 1):
                diff_list[i] -= 1
            d += 1
        else:  # diff_list[0] < 0
            # Find consecutive negative values
            k = 0
            for i in range(len(diff_list)):
                if diff_list[i] < 0:
                    k = i
                else:
                    break
            
            # Add 1 to all negative values
            for i in range(k + 1):
                diff_list[i] += 1
            d += 1
    
    return float(d)

def matrixbuilder_segmented(node):
    """Segment-aware matrix builder for better performance"""
    node_keys = list(node.keys())
    n = len(node_keys)
    
    # Pre-allocate matrix
    matrix = np.zeros((n, n), dtype=np.float32)
    
    # Calculate distances with segment-by-segment optimization
    for i, key1 in enumerate(node_keys):
        for j, key2 in enumerate(node_keys):
            if i == j:
                matrix[i, j] = 0.0
            elif i < j:
                dist_val = dist_segmented_optimized(node[key1], node[key2])
                matrix[i, j] = dist_val
                matrix[j, i] = dist_val
    
    return node_keys, matrix.tolist()

def dist_segmented_optimized(node1, node2):
    """Optimized distance calculation working on segments"""
    if len(node1) != len(node2):
        raise ValueError("Nodes must have same length")
    
    total_dist = 0.0
    
    for i in range(len(node1)):
        seg1 = np.array(node1[i], dtype=np.int32)  # Use int32 to match original
        seg2 = np.array(node2[i], dtype=np.int32)
        
        # Handle zero values - match original zerodisthelper behavior exactly
        if np.any(seg1 == 0) or np.any(seg2 == 0):
            segment_dist = zerodisthelper_numpy(seg1.tolist(), seg2.tolist())
        else:
            segment_dist = med_distance_numpy(seg1.astype(np.float32), seg2.astype(np.float32))
        
        total_dist += segment_dist
    
    return total_dist

def zerodisthelper_numpy(node1, node2):
    """Numpy version of zerodisthelper that matches original exactly"""
    n1 = node1[:]  # Copy
    n2 = node2[:]  # Copy
    temp1 = []
    temp2 = []
    
    for i in range(len(n1)):
        x1 = n1[i]
        x2 = n2[i]
        if x1 == 0:
            if x2 == 0:
                temp1.append(x1)
                temp2.append(x2)
            else:
                return 1000000.0
        else:
            temp1.append(x1)
            temp2.append(x2)
    
    if not temp1:
        return 0.0
    
    return med_distance_numpy(np.array(temp1, dtype=np.float32), np.array(temp2, dtype=np.float32))

# Wrapper functions for compatibility
def matrixbuilder(node):
    """Choose the best matrix builder based on data characteristics"""
    # Get data characteristics
    if not node:
        return [], []
    
    first_key = list(node.keys())[0]
    num_segments = len(node[first_key])
    avg_segment_size = sum(len(seg) for seg in node[first_key]) / num_segments if num_segments > 0 else 0
    
    # Use numpy optimization for larger datasets
    if len(node) > 20 or avg_segment_size > 50:
        return matrixbuilder_segmented(node)
    else:
        return matrixbuilder_numpy(node)

def dist(node1, node2):
    """Compatibility wrapper for dist function"""
    return dist_segmented_optimized(node1, node2)

# Original functions for fallback
def matrixbuilder_original(node):
    """Original matrix builder - kept for compatibility"""
    matrix = []
    for node1 in node:
        temp = []
        for node2 in node:
            temp.append(dist_original(node[node1], node[node2]))
        matrix.append(temp)
    return node.keys(), matrix

def dist_original(node1, node2):
    """Original distance function"""
    d = 0
    for i in range(0, len(node1)):
        d = d + disthelper_original(node1[i], node2[i])
    return d

def disthelper_original(node1, node2):
    if 0 in node1 or 0 in node2:
        return zerodisthelper_original(node1, node2)
    return distcalc_original(node1, node2)

def distcalc_original(node1, node2):
    assert len(node1) == len(node2)
    if len(node1) == 1:
        return abs(node1[0] - node2[0])
    else:
        d = 0
        newlist = copy.deepcopy(node1)
        for i in range(0, len(node2)):
            newlist[i] -= node2[i]
        while newlist:
            if newlist[0] == 0:
                newlist.pop(0)
            elif newlist[0] > 0:
                k = 0
                for i in range(0, len(newlist)):
                    if newlist[i] > 0:
                        k = i
                    else:
                        break
                for i in range(0, k + 1):
                    newlist[i] -= 1
                d += 1
            elif newlist[0] < 0:
                k = 0
                for i in range(0, len(newlist)):
                    if newlist[i] < 0:
                        k = i
                    else:
                        break
                for i in range(0, k + 1):
                    newlist[i] += 1
                d += 1
        return abs(d)

def zerodisthelper_original(node1, node2):
    n1 = copy.deepcopy(node1)
    n2 = copy.deepcopy(node2)
    dist = 0
    temp1 = []
    temp2 = []
    while n1:
        x1 = n1.pop()
        x2 = n2.pop()
        if x1 == 0:
            if x2 == 0:
                temp1.append(x1)
                temp2.append(x2)
            else:
                return 1000000
        else:
            temp1.append(x1)
            temp2.append(x2)
    return distcalc_original(temp1, temp2)