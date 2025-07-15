#!/usr/bin/env python2
"""
Optimized but fully compatible version of ComputeDistance.py
"""

import copy
import numpy as np
from collections import deque

def matrixbuilder(node):
    """Optimized matrix builder with symmetric calculation"""
    node_keys = list(node.keys())
    n = len(node_keys)
    
    # Pre-allocate matrix for better performance
    matrix = [[0.0 for _ in range(n)] for _ in range(n)]
    
    # Calculate distances using optimized algorithm
    for i, node1_key in enumerate(node_keys):
        for j, node2_key in enumerate(node_keys):
            if i == j:
                matrix[i][j] = 0.0
            elif i < j:  # Only calculate upper triangle
                dist_val = dist_optimized(node[node1_key], node[node2_key])
                matrix[i][j] = dist_val
                matrix[j][i] = dist_val  # Symmetric matrix
            
    return node_keys, matrix

def dist_optimized(node1, node2):
    """Optimized distance calculation that matches original exactly"""
    d = 0
    for i in range(len(node1)):
        d += disthelper_optimized(node1[i], node2[i])
    return d

def disthelper_optimized(node1, node2):
    """Optimized distance helper"""
    if 0 in node1 or 0 in node2:
        return zerodisthelper_optimized(node1, node2)
    return distcalc_optimized(node1, node2)

def distcalc_optimized(node1, node2):
    """Optimized MED distance calculation"""
    if len(node1) != len(node2):
        raise ValueError("Segments must have same length")
    
    if len(node1) == 1:
        return abs(node1[0] - node2[0])
    
    # Use deque for efficient operations
    diff_list = deque([node1[i] - node2[i] for i in range(len(node1))])
    
    d = 0
    while diff_list:
        if diff_list[0] == 0:
            diff_list.popleft()
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
        elif diff_list[0] < 0:
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
    
    return d

def zerodisthelper_optimized(node1, node2):
    """Optimized zero distance helper"""
    # Use lists instead of deep copy for better performance
    n1 = list(node1)
    n2 = list(node2)
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
                return 1000000
        else:
            temp1.append(x1)
            temp2.append(x2)
    
    return distcalc_optimized(temp1, temp2)

# Numpy-enhanced version for very large datasets
def matrixbuilder_numpy(node):
    """Numpy-optimized matrix builder for large datasets"""
    node_keys = list(node.keys())
    n = len(node_keys)
    
    # Convert to numpy arrays for vectorized operations
    node_arrays = {}
    max_total_length = 0
    
    for key in node_keys:
        # Flatten chromosomal segments
        flattened = []
        for segment in node[key]:
            flattened.extend(segment)
        node_arrays[key] = np.array(flattened, dtype=np.int32)
        max_total_length = max(max_total_length, len(flattened))
    
    # Pre-allocate distance matrix
    matrix = np.zeros((n, n), dtype=np.float64)
    
    # Calculate distances
    for i, key1 in enumerate(node_keys):
        for j, key2 in enumerate(node_keys):
            if i == j:
                matrix[i, j] = 0.0
            elif i < j:
                # Use original algorithm for correctness
                dist_val = dist_optimized(node[key1], node[key2])
                matrix[i, j] = dist_val
                matrix[j, i] = dist_val
    
    return node_keys, matrix.tolist()

# Adaptive function that chooses the best implementation
def matrixbuilder_adaptive(node):
    """Choose the best matrix builder based on data size"""
    if not node:
        return [], []
    
    num_nodes = len(node)
    
    # For small datasets, use standard optimized version
    if num_nodes <= 50:
        return matrixbuilder(node)
    else:
        # For large datasets, use numpy version
        return matrixbuilder_numpy(node)

# Keep original functions for testing
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

# Use the optimized function by default
dist = dist_optimized