#!/usr/bin/env python3
"""
Simple text-based tree structure viewer for debugging MEDALT trees
"""

import pandas as pd
import sys
from collections import defaultdict, deque

def build_tree_structure(tree_file):
    """Build parent-child relationships from tree file"""
    df = pd.read_csv(tree_file, sep='\t')
    
    children = defaultdict(list)
    parents = {}
    
    for _, row in df.iterrows():
        parent = row['from']
        child = row['to']
        distance = row['dist']
        
        children[parent].append((child, distance))
        parents[child] = parent
    
    return children, parents

def print_tree(children, root='root', indent=0, max_depth=5):
    """Print tree structure recursively"""
    if indent > max_depth:
        print('  ' * indent + '... (truncated)')
        return
    
    # Sort children by distance for consistent output
    child_list = sorted(children.get(root, []), key=lambda x: x[1])
    
    for child, distance in child_list:
        print('  ' * indent + f'├─ {child} (dist: {distance})')
        print_tree(children, child, indent + 1, max_depth)

def analyze_tree_stats(tree_file):
    """Analyze basic tree statistics"""
    df = pd.read_csv(tree_file, sep='\t')
    
    print("=== Tree Statistics ===")
    print(f"Total edges: {len(df)}")
    print(f"Unique nodes: {len(set(df['from']) | set(df['to']))}")
    
    # Distance distribution
    distances = df['dist'].value_counts().sort_index()
    print(f"\nDistance distribution:")
    for dist, count in distances.items():
        print(f"  Distance {dist}: {count} edges")
    
    # Root analysis
    roots = set(df['from']) - set(df['to'])
    print(f"\nRoot nodes: {roots}")
    
    # Leaf analysis (nodes with no children)
    leaves = set(df['to']) - set(df['from'])
    print(f"Leaf nodes: {len(leaves)} total")
    
    # Branching analysis
    children_count = df['from'].value_counts()
    print(f"\nBranching statistics:")
    print(f"  Nodes with 1 child: {(children_count == 1).sum()}")
    print(f"  Nodes with 2+ children: {(children_count >= 2).sum()}")
    print(f"  Max children per node: {children_count.max()}")
    
    return children_count

def main():
    if len(sys.argv) < 2:
        print("Usage: python tree_viewer.py <tree_file.txt>")
        sys.exit(1)
    
    tree_file = sys.argv[1]
    
    print(f"Analyzing tree: {tree_file}")
    print("=" * 50)
    
    # Basic statistics
    analyze_tree_stats(tree_file)
    
    print("\n" + "=" * 50)
    print("Tree Structure (showing root and first few levels):")
    print("=" * 50)
    
    # Build and display tree structure
    children, parents = build_tree_structure(tree_file)
    
    # Find root
    roots = [node for node in children.keys() if node not in parents]
    
    if not roots:
        print("No root found! This might indicate a cycle in the tree.")
        return
    
    for root in roots:
        print(f"\nFrom root: {root}")
        print_tree(children, root, 0, max_depth=4)
    
    print(f"\n(Tree truncated at depth 4 for readability)")

if __name__ == "__main__":
    main()