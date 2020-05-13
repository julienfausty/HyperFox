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

resultsDir = '../results/DiffusionConvergence/'
filename = 'ImpImp/Breakdown.csv'

simData = readcsv(resultsDir + filename)

structData = {}

for sim in simData:
    key = ('ImpImp', sim['dim'], sim['h'], sim['timeStep'])
    if key not in structData:
        structData[key] = {}
        structData[key]['p'] = []
        structData[key]['elinAlg'] = []
        structData[key]['el2'] = []
        structData[key]['del2'] = []
        structData[key]['T'] = []
        structData[key]['setup'] = []
        structData[key]['assembly'] = []
        structData[key]['resolution'] = []
        structData[key]['post'] = []
    structData[key]['p'].append(float(sim['order']))
    structData[key]['elinAlg'].append(float(sim['linAlgErr']))
    structData[key]['el2'].append(float(sim['l2Err']))
    structData[key]['del2'].append(float(sim['dL2Err']))
    structData[key]['T'].append(float(sim['runtime']))
    structData[key]['setup'].append(float(sim['setup']))
    structData[key]['assembly'].append(float(sim['assembly']))
    structData[key]['resolution'].append(float(sim['resolution']))
    structData[key]['post'].append(float(sim['post']))


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

plotKey = ('ImpImp', '2', '5e-2', '1e-2')

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

plt.xlabel(r'$p$')
plt.ylabel(r'$T$')

plt.bar(structData[plotKey]['p'], structData[plotKey]['T'], width=0.5, color=tableau20[0], label='Total runtime')
plt.bar(structData[plotKey]['p'], structData[plotKey]['resolution'], width=0.5, color=tableau20[2], label='Resolution')
plt.bar(structData[plotKey]['p'], structData[plotKey]['assembly'], bottom=structData[plotKey]['resolution'], width=0.5, color=tableau20[1], label='Assembly')

print(np.array(structData[plotKey]['assembly'])/np.array(structData[plotKey]['T']))
print(np.array(structData[plotKey]['resolution'])/np.array(structData[plotKey]['T']))

plt.legend()

plt.show()


