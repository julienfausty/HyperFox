import numpy as np
import matplotlib.pyplot as plt
import csv

def readcsv(filename):
    data = []
    f = open(filename)
    myReader = csv.DictReader(f, delimiter=',')
    for line in myReader:
        data.append(line)
    return data

resultsDir = '../results/Burgers/HDG/'
filename = 'Imp/BEuler/Breakdown.csv'

simData = readcsv(resultsDir + filename)

structData = {}

for sim in simData:
    key = ('Imp', sim['dim'], sim['order'])
    if key not in structData:
        structData[key] = {}
        structData[key]['h'] = []
        structData[key]['dt'] = []
        structData[key]['elinAlg'] = []
        structData[key]['el2'] = []
        structData[key]['del2'] = []
        structData[key]['T'] = []
    structData[key]['h'].append(float(sim['h']))
    structData[key]['dt'].append(float(sim['dt']))
    structData[key]['elinAlg'].append(float(sim['linAlgErr']))
    structData[key]['el2'].append(float(sim['l2Err']))
    structData[key]['del2'].append(float(sim['dL2Err']))
    structData[key]['T'].append(float(sim['avgRuntime']))


# These are the "Tableau 20" colors as RGB.    
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]    
  
# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.)

plt.figure()

# Remove the plot frame lines. They are unnecessary chartjunk.
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["bottom"].set_visible(True)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(True)

# Ensure that the axis ticks only show up on the bottom and left of the plot.
# Ticks on the right and top of the plot are generally unnecessary chartjunk.
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

plt.xlabel(r'$\log(\Delta t)$')
# plt.xlabel(r'$\log(h)$')
plt.ylabel(r'$\log(e)$')

variable = 'h'
errType = 'el2'
slicer = 'dt'
valSlice = 1e-4
tol = 0.5


for key in structData:
    structData[key]['sliceVar'] = [structData[key][variable][i] for i in range(len(structData[key][variable])) if (structData[key][slicer][i] == valSlice and structData[key][errType][i] < tol)]
    structData[key]['sliceErr'] = [structData[key][errType][i] for i in range(len(structData[key][errType])) if (structData[key][slicer][i] == valSlice and structData[key][errType][i] < tol)]
    structData[key]['slope'] = [(np.log10(structData[key]['sliceErr'][i]) - np.log10(structData[key]['sliceErr'][i-1]))/(np.log10(structData[key]['sliceVar'][i]) - np.log10(structData[key]['sliceVar'][i-1])) for i in range(1, len(structData[key]['sliceVar']))]


for key in sorted(structData):
    plt.scatter(np.log10(structData[key]['sliceVar']), np.log10(structData[key]['sliceErr']), color=tableau20[int(key[2])])   
    plt.plot(np.log10(structData[key]['sliceVar']), np.log10(structData[key]['sliceErr']), color=tableau20[int(key[2])], label=key[0] + '(p=' + key[2] + ', d=' + key[1]+')')
    for i in range(len(structData[key]['sliceVar'])-1):
        plt.text(np.log10(structData[key]['sliceVar'][i]), np.log10(structData[key]['sliceErr'][i]), '%3.2f' % (structData[key]['slope'][i]), color=tableau20[int(key[2])])

plt.legend()

plt.show()


