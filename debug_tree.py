#!/usr/bin/env python3

from src.tree_utils import read_cnv, dist

# Read the segmental CNV data
nodes, root = read_cnv('example/outputRNA_fixed/2_scRNA.CNV_bin_30.csv')

# Get two real cells (not root)
cells = [name for name in nodes.keys() if name != 'root']
cell1, cell2 = cells[0], cells[1]

print('Cell 1:', cell1)
print('Data 1:', nodes[cell1])
print('Cell 2:', cell2)
print('Data 2:', nodes[cell2])

# Calculate distance
distance = dist(nodes[cell1], nodes[cell2])
print('Distance between', cell1, 'and', cell2, ':', distance)

# Test a few more pairs
for i in range(min(5, len(cells)-1)):
    d = dist(nodes[cells[i]], nodes[cells[i+1]])
    print('Distance', cells[i], '->', cells[i+1], ':', d)

# Check if all cells have the same data (explaining 0 distances)
print('\nChecking for identical cells:')
first_cell_data = nodes[cells[0]]
identical_count = 0
for cell in cells:
    if nodes[cell] == first_cell_data:
        identical_count += 1

print(f'Number of cells identical to first cell: {identical_count}/{len(cells)}')

# Show first few cells' data
print('\nFirst 3 cells data:')
for i in range(min(3, len(cells))):
    print(f'{cells[i]}: {nodes[cells[i]]}')