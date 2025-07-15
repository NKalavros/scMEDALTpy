#!/usr/bin/env python2
"""
Direct test of distance calculation scalability on large datasets
"""

import time
import sys
sys.path.insert(0, '.')

def create_test_nodes_from_file(filename, max_cells=None):
    """Create test nodes from synthetic data file"""
    print "Loading data from %s..." % filename
    
    nodes = {}
    
    with open(filename, 'r') as f:
        # Read header
        header = f.readline().strip().split('\t')
        cell_names = header[1:]  # Skip first empty column
        
        if max_cells:
            cell_names = cell_names[:max_cells]
        
        print "Processing %d cells..." % len(cell_names)
        
        # Initialize nodes
        for cell_name in cell_names:
            nodes[cell_name] = []
        
        # Read gene data
        gene_count = 0
        for line in f:
            if not line.strip():
                continue
                
            parts = line.strip().split('\t')
            gene_name = parts[0]
            values = parts[1:]
            
            # Convert to copy numbers (round to integers)
            for i, cell_name in enumerate(cell_names):
                if i < len(values):
                    cn = int(round(float(values[i]) * 2))  # Convert relative to absolute
                    cn = max(0, min(cn, 6))  # Clamp to reasonable range
                    
                    # Add to appropriate chromosome segment
                    # Simulate chromosomal segments by grouping genes
                    seg_idx = gene_count / 100  # 100 genes per segment
                    
                    # Ensure segment exists
                    while len(nodes[cell_name]) <= seg_idx:
                        nodes[cell_name].append([])
                    
                    nodes[cell_name][seg_idx].append(cn)
            
            gene_count += 1
            
            # Progress indicator
            if gene_count % 500 == 0:
                print "  Processed %d genes..." % gene_count
    
    print "Loaded %d genes into %d cells" % (gene_count, len(cell_names))
    print "Created %d chromosomal segments per cell" % len(nodes[cell_names[0]])
    
    return nodes

def test_distance_calculation_scaling():
    """Test distance calculation with increasing cell counts"""
    
    filename = "synthetic_data_3000g_2000c.txt"
    print "Large-Scale Distance Calculation Test"
    print "=" * 50
    
    # Test with increasing numbers of cells
    cell_counts = [50, 100, 200, 500, 1000]
    
    from ComputeDistance import matrixbuilder
    
    results = []
    
    for num_cells in cell_counts:
        print "\\nTesting with %d cells..." % num_cells
        
        # Load subset of data
        nodes = create_test_nodes_from_file(filename, max_cells=num_cells)
        
        if not nodes:
            print "Failed to load data!"
            continue
        
        # Test distance calculation
        print "Running distance matrix calculation..."
        start_time = time.time()
        
        try:
            keys, matrix = matrixbuilder(nodes)
            end_time = time.time()
            elapsed = end_time - start_time
            
            print "Success! Time: %.3f seconds" % elapsed
            
            # Verify matrix properties
            n = len(keys)
            if len(matrix) == n and len(matrix[0]) == n:
                print "Matrix size: %dx%d" % (n, n)
                
                # Check some properties
                is_symmetric = True
                for i in range(min(10, n)):
                    for j in range(min(10, n)):
                        if abs(matrix[i][j] - matrix[j][i]) > 1e-10:
                            is_symmetric = False
                            break
                    if not is_symmetric:
                        break
                
                print "Matrix is symmetric: %s" % is_symmetric
                print "Sample distances: [%.1f, %.1f, %.1f]" % (matrix[0][1], matrix[0][2], matrix[1][2])
            
            results.append({
                'num_cells': num_cells,
                'time': elapsed,
                'success': True
            })
            
        except Exception as e:
            print "Failed: %s" % str(e)
            results.append({
                'num_cells': num_cells,
                'time': float('inf'),
                'success': False
            })
    
    # Summary
    print "\\n" + "=" * 50
    print "SCALING RESULTS"
    print "=" * 50
    print "%-10s %-12s %-12s" % ("Cells", "Time (s)", "Status")
    print "-" * 35
    
    for result in results:
        if result['success']:
            time_str = "%.3f" % result['time']
            status = "OK"
        else:
            time_str = "FAILED"
            status = "ERROR"
        
        print "%-10d %-12s %-12s" % (result['num_cells'], time_str, status)
    
    # Extrapolation
    successful_results = [r for r in results if r['success']]
    if len(successful_results) >= 2:
        print "\\nExtrapolation to larger scales:"
        
        # Simple quadratic fit
        largest = successful_results[-1]
        base_cells = largest['num_cells']
        base_time = largest['time']
        
        for target_cells in [2000, 5000, 10000]:
            scale_factor = float(target_cells) / base_cells
            estimated_time = base_time * (scale_factor ** 2)
            
            print "  %d cells: ~%.1f seconds (%.1f minutes)" % (
                target_cells, estimated_time, estimated_time / 60.0
            )
    
    return results

if __name__ == "__main__":
    results = test_distance_calculation_scaling()
    
    # Check if we can handle large scales
    success_count = sum(1 for r in results if r['success'])
    if success_count >= 3:
        print "\\nDistance calculation optimization is working well!"
        print "Ready for next-level optimizations."
    else:
        print "\\nDistance calculation needs more optimization."