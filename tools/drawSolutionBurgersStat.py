import h5py
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.tri as mtri
import matplotlib.cm as cm
import scipy.interpolate as sciip
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
def plot_streamlines(P, u_x, u_y):
    grid_x, grid_y = np.mgrid[0:1:1000j, 0:1:1000j]
    grid_u_x = sciip.griddata(P, u_x, (grid_x,grid_y), method='cubic')
    grid_u_y = sciip.griddata(P, u_y, (grid_x,grid_y), method='cubic')
    c = np.hypot(grid_u_x.T, grid_u_y.T)
    plt.streamplot(grid_x[:,0], grid_y[0,:], grid_u_x.T, grid_u_y.T, density=3.0, color=c)

n = 1
order = '3'
meshSize = '7e-2'

dirName = '../results/BurgersStat/HDG/'

subDirName = 'regression_dim-2_h-' + meshSize + '_ord-' + order + '/'

fileName = 'res_' + str(n) +'.h5'

f = h5py.File(dirName + subDirName + fileName, 'r')

meshData = readGroup(f)


f.close()

nodes = meshData['/Mesh']['/Mesh/Nodes']
cells = meshData['/Mesh']['/Mesh/Cells']
cellShape = cells.shape

sol = meshData['/FieldData']['/FieldData/Solution']
anasol = meshData['/FieldData']['/FieldData/Analytical']
residual = meshData['/FieldData']['/FieldData/Residual']

# anasol = np.array([analyticalSol(t, node) for node in nodes])

# plt.figure(figsize=(12, 9))
plt.figure(figsize=(21, 5))

# plt.scatter(nodes[:, 0], nodes[:, 1], color='red')

x = np.zeros(cellShape[0]*cellShape[1])
y = np.zeros(cellShape[0]*cellShape[1])
P = np.zeros((cellShape[0]*cellShape[1], 2))

for i in range(cellShape[0]):
    cell = [int(ind) for ind in cells[i, :]]
    x[(i*cellShape[1]):((i+1)*cellShape[1])] = nodes[cell, 0]
    y[(i*cellShape[1]):((i+1)*cellShape[1])] = nodes[cell, 1]
    P[(i*cellShape[1]):((i+1)*cellShape[1])] = nodes[cell]

s = np.zeros(cellShape[0]*cellShape[1])
diff = np.zeros(cellShape[0]*cellShape[1])

for i in range(cellShape[0]):
    s[i*cellShape[1]:(i+1)*cellShape[1]] = [val[0] for val in sol[i,:]] #DG
    diff[i*cellShape[1]:(i+1)*cellShape[1]] = [val[1] for val in sol[i,:]] #DG
    # cell = [int(ind) for ind in cells[i, :]] #CG
    # s[i*cellShape[1]:(i+1)*cellShape[1]] = [val[0] for val in sol[cell,0]] #CG
    # diff[i*cellShape[1]:(i+1)*cellShape[1]] = [val[1] for val in sol[cell,0]] #CG

plt.subplot(131)

c = np.hypot(s, diff)

# plt.quiver(x, y, s, diff, c)
# plt.colorbar()

plot_streamlines(P, s, diff)
plt.colorbar()

plt.title('Solution')
# plt.clim(0.03, 0.045)

for i in range(cellShape[0]):
    s[i*cellShape[1]:(i+1)*cellShape[1]] = [val[0] for val in anasol[i,:]] #DG
    diff[i*cellShape[1]:(i+1)*cellShape[1]] = [val[1] for val in anasol[i,:]] #DG
    # cell = [int(ind) for ind in cells[i, :]] #CG
    # s[i*cellShape[1]:(i+1)*cellShape[1]] = [val[0] for val in anasol[cell,0]] #CG
    # diff[i*cellShape[1]:(i+1)*cellShape[1]] = [val[1] for val in anasol[cell,0]] #CG

plt.subplot(132)

c = np.hypot(s, diff)

# plt.quiver(x, y, s, diff, c)
# plt.colorbar()

plot_streamlines(P, s, diff)
plt.colorbar()

plt.title('Analytical')


for i in range(cellShape[0]):
    s[i*cellShape[1]:(i+1)*cellShape[1]] = [val[0] for val in residual[i,:]] #DG
    # cell = [int(ind) for ind in cells[i, :]] #CG
    # s[i*cellShape[1]:(i+1)*cellShape[1]] = [val[0] for val in residual[cell,0]] #CG

plt.subplot(133)

plt.tricontourf(x, y, s, 500, cmap=cm.get_cmap(name='inferno'))

plt.colorbar()

plt.title('Residual')

# for i in range(cellShape[0]):
    # cell = [int(index) for index in cells[i, [0, 1, 2, 0]]]
    # plt.plot(nodes[cell, 0], nodes[cell, 1], color='black')

plt.show() 
