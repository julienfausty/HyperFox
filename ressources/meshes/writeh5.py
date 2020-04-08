import h5py
import numpy as np


dirName = ''

meshName = 'fieldTest.h5'

f = h5py.File(dirName + meshName, 'w')

# meshGrp = f.create_group('Mesh')

# nodes = meshGrp.create_dataset('Nodes', (9, 2), 'f8')
# nodes[0, 0] = 0.0
# nodes[0, 1]  = 0.0
# nodes[1, 0] = 1.0
# nodes[1, 1] = 0.0
# nodes[2, 0] = 1.0
# nodes[2, 1]= 1.0
# nodes[3, 0] = 0.0
# nodes[3, 1]= 1.0
# nodes[4, 0] = 0.5
# nodes[4, 1]= 0.5
# nodes[5, 0] = 0.5
# nodes[5, 1]= 0.0
# nodes[6, 0] = 1.0
# nodes[6, 1]= 0.5
# nodes[7, 0] = 0.5
# nodes[7, 1]= 1.0
# nodes[8, 0] = 0.0
# nodes[8, 1]= 0.5

# Tri
# cells = meshGrp.create_dataset('Cells', (8, 3), dtype='i4')
# cells[0, :] = [0, 5, 4]
# cells[1, :] = [5, 1, 4]
# cells[2, :] = [1, 6, 4]
# cells[3, :] = [6, 2, 4]
# cells[4, :] = [2, 7, 4]
# cells[5, :] = [7, 3, 4]
# cells[6, :] = [3, 8, 4]
# cells[7, :] = [8, 0, 4]

# Quad
# cells = meshGrp.create_dataset('Cells', (4, 4), dtype='i4')
# cells[0, :] = [0, 5, 4, 8]
# cells[1, :] = [5, 1, 6, 4]
# cells[2, :] = [4, 6, 2, 7]
# cells[3, :] = [8, 4, 7, 3]

# Tri 2
# cells = meshGrp.create_dataset('Cells', (2, 6), dtype='i4')
# cells[0, :] = [0, 1, 3, 5, 4, 8]
# cells[1, :] = [3, 1, 2, 4, 6, 7]

# Quad 2
# cells = meshGrp.create_dataset('Cells', (1, 9), 'i4')
# cells[0, :] = [0, 1, 2, 3, 5, 6, 7, 8, 4]

fieldGrp = f.create_group('FieldData')

nodeField = fieldGrp.create_dataset('NodeField', (9,1,1), 'f8')
for i in range(9):
    nodeField[i, 0, 0] = i
nodeField.attrs['ftype'] = 0

cellField = fieldGrp.create_dataset('CellField', (8,2,1), 'f8')
for i in range(8):
    cellField[i, 0, 0] = i
    cellField[i, 1, 0] = -i
cellField.attrs['ftype'] = 3

f.close()
