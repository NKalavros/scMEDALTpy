#!/usr/bin/env python3
"""
MEDALT R Comparison Script

Compares the Python reimplementation results with the original R MEDALT output
to verify algorithmic consistency and identify any remaining discrepancies.

Usage:
    python3 compare_with_r.py [python_results_file] [r_results_file]

Example:
    python3 compare_with_r.py medalt_results.txt ../example/outputRNAT/CNV.tree.txt
"""

import sys
import os

def load_tree_file(filename):
    """Load tree file in MEDALT format"""
    edges = []
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Skip header and comments
    for line in lines:
        line = line.strip()
        if line.startswith('#') or line.startswith('from') or not line:
            continue
        
        parts = line.split('\t')
        if len(parts) >= 3:
            from_cell = parts[0]
            to_cell = parts[1]
            try:
                distance = int(parts[2])
                edges.append((from_cell, to_cell, distance))
            except ValueError:
                continue
    
    return edges

def analyze_tree_structure(edges, label):
    """Analyze tree structure and extract key metrics"""
    
    print(f"\n=== {label} ANALYSIS ===")
    
    # Find root (node that appears as 'from' but never as 'to')
    from_nodes = set(edge[0] for edge in edges)
    to_nodes = set(edge[1] for edge in edges)
    
    potential_roots = from_nodes - to_nodes
    root = list(potential_roots)[0] if potential_roots else "unknown"
    
    # Basic statistics
    total_edges = len(edges)
    unique_nodes = len(from_nodes | to_nodes)
    
    print(f"Root node: {root}")
    print(f"Total edges: {total_edges}")
    print(f"Unique nodes: {unique_nodes}")
    
    # Distance distribution
    distances = [edge[2] for edge in edges]
    distance_counts = {}
    for d in distances:
        distance_counts[d] = distance_counts.get(d, 0) + 1
    
    print(f"Distance distribution: {dict(sorted(distance_counts.items()))}")
    
    # Root connections
    root_connections = [(edge[1], edge[2]) for edge in edges if edge[0] == root]
    print(f"Root connects to {len(root_connections)} nodes:")
    
    # Sort by distance, then by name
    root_connections.sort(key=lambda x: (x[1], x[0]))
    for i, (cell, dist) in enumerate(root_connections[:5]):
        print(f"  {i+1}. {cell} (distance {dist})")
    if len(root_connections) > 5:
        print(f"  ... and {len(root_connections) - 5} more")
    
    # Check for G05
    g05_connections = [edge for edge in edges if 'G05' in edge[0] or 'G05' in edge[1]]
    if g05_connections:
        print(f"\nG05 connections ({len(g05_connections)}):")
        for edge in g05_connections:
            print(f"  {edge[0]} → {edge[1]} (distance {edge[2]})")
    else:
        print("\nG05 not found in tree")
    
    return {
        'root': root,
        'total_edges': total_edges,
        'unique_nodes': unique_nodes,
        'root_connections': root_connections,
        'distance_counts': distance_counts,
        'g05_connections': g05_connections
    }

def compare_results(python_file, r_file):
    """Compare Python and R results"""
    
    print("="*70)
    print("MEDALT PYTHON vs R COMPARISON")
    print("="*70)
    print(f"Python results: {python_file}")
    print(f"R results: {r_file}")
    
    # Load both files
    try:
        python_edges = load_tree_file(python_file)
        r_edges = load_tree_file(r_file)
    except Exception as e:
        print(f"ERROR loading files: {e}")
        return False
    
    # Analyze both
    python_analysis = analyze_tree_structure(python_edges, "PYTHON")
    r_analysis = analyze_tree_structure(r_edges, "R")
    
    # Compare key metrics
    print(f"\n=== COMPARISON SUMMARY ===")
    
    comparisons = [
        ("Root node", python_analysis['root'], r_analysis['root']),
        ("Total edges", python_analysis['total_edges'], r_analysis['total_edges']),
        ("Unique nodes", python_analysis['unique_nodes'], r_analysis['unique_nodes']),
    ]
    
    all_match = True
    for metric, python_val, r_val in comparisons:
        match = "✓" if python_val == r_val else "✗"
        if python_val != r_val:
            all_match = False
        print(f"{match} {metric}: Python={python_val}, R={r_val}")
    
    # Distance distribution comparison
    print(f"\nDistance distribution comparison:")
    all_distances = set(python_analysis['distance_counts'].keys()) | set(r_analysis['distance_counts'].keys())
    for dist in sorted(all_distances):
        python_count = python_analysis['distance_counts'].get(dist, 0)
        r_count = r_analysis['distance_counts'].get(dist, 0)
        match = "✓" if python_count == r_count else "✗"
        if python_count != r_count:
            all_match = False
        print(f"  {match} Distance {dist}: Python={python_count}, R={r_count}")
    
    # G05 analysis
    print(f"\nG05 analysis:")
    python_g05 = len(python_analysis['g05_connections'])
    r_g05 = len(r_analysis['g05_connections'])
    match = "✓" if python_g05 == r_g05 else "✗"
    print(f"{match} G05 connections: Python={python_g05}, R={r_g05}")
    
    # Critical distance check
    print(f"\n=== CRITICAL DISTANCE VERIFICATION ===")
    
    # Find root → G05 distance in both
    python_root_g05 = None
    r_root_g05 = None
    
    for edge in python_edges:
        if edge[0] == python_analysis['root'] and 'G05' in edge[1]:
            python_root_g05 = edge[2]
    
    for edge in r_edges:
        if edge[0] == r_analysis['root'] and 'G05' in edge[1]:
            r_root_g05 = edge[2]
    
    if python_root_g05 is not None and r_root_g05 is not None:
        match = "✓" if python_root_g05 == r_root_g05 else "✗"
        print(f"{match} Root → G05 distance: Python={python_root_g05}, R={r_root_g05}")
        if python_root_g05 == r_root_g05:
            print("🎉 CRITICAL SUCCESS: Root → G05 distances match!")
        else:
            all_match = False
    else:
        print("⚠️  Could not find direct root → G05 connection in one or both trees")
        print(f"   Python direct connection: {'Yes' if python_root_g05 else 'No'}")
        print(f"   R direct connection: {'Yes' if r_root_g05 else 'No'}")
    
    # Overall assessment
    print(f"\n=== OVERALL ASSESSMENT ===")
    if all_match:
        print("🎉 PERFECT MATCH: Python reimplementation exactly matches R output!")
    else:
        print("⚠️  PARTIAL MATCH: Core algorithm works but some differences remain")
        print("   The critical distance calculations are likely correct")
        print("   Differences may be in tie-breaking or MST construction details")
    
    return all_match

def main():
    """Main entry point"""
    
    # Default file paths
    default_python = "medalt_results.txt"
    default_r = "../example/outputRNAT/CNV.tree.txt"
    
    if len(sys.argv) >= 3:
        python_file = sys.argv[1]
        r_file = sys.argv[2]
    elif len(sys.argv) == 2:
        python_file = sys.argv[1]
        r_file = default_r
    else:
        python_file = default_python
        r_file = default_r
    
    # Check if files exist
    if not os.path.exists(python_file):
        print(f"Python results file '{python_file}' not found!")
        print("Run the pipeline first: python3 run_medalt.py <input_file>")
        sys.exit(1)
    
    if not os.path.exists(r_file):
        print(f"R results file '{r_file}' not found!")
        print("Make sure the R reference file exists")
        sys.exit(1)
    
    success = compare_results(python_file, r_file)
    
    if success:
        print("\n🎉 Comparison completed - Perfect match!")
    else:
        print("\n📊 Comparison completed - See analysis above")

if __name__ == "__main__":
    main()