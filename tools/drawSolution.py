import h5py
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib.cm as cm
import copy

def readGroup(grp):
    meshData = {}
    if(isinstance(grp, h5py.Group)):
        for k in grp.keys():
            if(isinstance(grp[k], h5py.Dataset)):
                meshData[grp[k].name] = np.empty(grp[k].shape)
                grp[k].read_direct(meshData[grp[k].name])
                for p in grp[k].attrs.keys():
                    meshData[grp[k].name + ".Attrs"] = (p, grp[k].attrs[p])
            else:
                meshData[grp[k].name] = readGroup(grp[k])
    else:
        meshData[grp.name] = grp
    return meshData

def analyticalSol(t, p):
    trigPart = np.cos(np.pi*p[0]/2.0)*np.cos(np.pi*p[1]/2.0) - np.sin(np.pi*p[0]/2.0)*np.sin(np.pi*p[1]/2.0)
    timePart = np.exp(-2.0*((np.pi/2.0)**2.0)*t)
    srcPart = 0
    for i in range(2):
        srcPart += (p[i] - 0.5) * math.erf(p[i] - 0.5) + np.exp(-(p[i] - 0.5)**2)/np.sqrt(np.pi)
    return (timePart*trigPart + srcPart)


dirName = '../results/DiffusionConvergence/ImpImp/'

subDirName = 'regression_dim-2_h-7e-2_ord-5_dt-5e-2/'

n = 10
timeStep = 5e-2

fileName = 'res_' + str(n) +'.h5'

f = h5py.File(dirName + subDirName + fileName, 'r')

meshData = readGroup(f)


f.close()

nodes = meshData['/Mesh']['/Mesh/Nodes']
cells = meshData['/Mesh']['/Mesh/Cells']
cellShape = cells.shape

sol = meshData['/FieldData']['/FieldData/Solution']

# sol = np.array([s[0, 0] for s in sol]) #CG

sol = np.array([s[:, 0] for s in sol]) #DG

anasol = np.array([analyticalSol(n*timeStep, node) for node in nodes])

# plt.figure(figsize=(7, 9))
plt.figure(figsize=(14, 5))

plt.subplot(121)

# plt.scatter(nodes[:, 0], nodes[:, 1], color='red')

# for i in range(cellShape[0]):
    # cell = [int(index) for index in cells[i, [0, 1, 2, 0]]]
    # plt.plot(nodes[cell, 0], nodes[cell, 1], color='black')

x = np.zeros(cellShape[0]*cellShape[1])
y = np.zeros(cellShape[0]*cellShape[1])
s = np.zeros(cellShape[0]*cellShape[1])
diff = np.zeros(cellShape[0]*cellShape[1])

for i in range(cellShape[0]):
    cell = [int(ind) for ind in cells[i, :]]
    x[(i*cellShape[1]):((i+1)*cellShape[1])] = nodes[cell, 0]
    y[(i*cellShape[1]):((i+1)*cellShape[1])] = nodes[cell, 1]
    s[(i*cellShape[1]):((i+1)*cellShape[1])] = sol[i] #DG
    # s[(i*cellShape[1]):((i+1)*cellShape[1])] = sol[cell] #GC

plt.tricontourf(x, y, s, 500)

plt.title('Solution')

plt.colorbar()
plt.clim(1.0, 1.5)

plt.subplot(122)

# for i in range(cellShape[0]):
    # cell = [int(index) for index in cells[i, [0, 1, 2, 0]]]
    # plt.plot(nodes[cell, 0], nodes[cell, 1], color='black')

for i in range(cellShape[0]):
    cell = [int(ind) for ind in cells[i, :]]
    diff[i*cellShape[1]:(i+1)*cellShape[1]] = np.abs(anasol[cell] - sol[i]) #DG
    # diff[i*cellShape[1]:(i+1)*cellShape[1]] = np.abs(anasol[cell] - sol[cell]) #CG

plt.tricontourf(x, y, diff, 500, cmap=cm.get_cmap(name='inferno'))

plt.title('Absolute Difference')

plt.colorbar()
plt.clim(0, 1e-2)

plt.show() 
