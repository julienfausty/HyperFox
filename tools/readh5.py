import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import copy

def readGroup(grp):
    print(grp.name)
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


dirName = '../ressources/meshes/'

meshName = 'lightSquareOrd1.h5'

f = h5py.File(dirName + meshName, 'r')

meshData = readGroup(f)

f.close()

print(meshData)
