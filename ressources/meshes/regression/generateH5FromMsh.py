import subprocess

meshSizes = {}

meshSizes['2'] = ['5e-1', '3e-1', '2e-1', '1e-1', '7e-2', '5e-2', '2e-2', '1e-2']
meshSizes['3'] = ['3e-1', '2e-1', '1e-1', '8e-2']

orders = range(1,6)


for k in meshSizes.keys():
    for h in meshSizes[k]:
        for o in orders:
            command = ["../../../build/bin/convertGmsh2H5HO"]
            command.append("-iregression_dim-" + k + "_h-" + h + ".msh")
            command.append("-d" + k)
            command.append("-r" + str(o))
            command.append("-oregression_dim-" + k + "_h-" + h + "_ord-" + str(o) + ".h5")
            subprocess.call(command)
