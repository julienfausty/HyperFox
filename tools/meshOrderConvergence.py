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

resultsDir = '../results/LaplaceConvergence/'
# filename = 'CG/CGLaplace.csv'

# simData = readcsv(resultsDir + filename)

structData = {}

# for sim in simData:
    # key = ('CG', sim['dim'], sim[' order'])
    # if key not in structData:
        # structData[key] = {}
        # structData[key]['h'] = []
        # structData[key]['elinAlg'] = []
        # structData[key]['el2'] = []
        # structData[key]['del2'] = []
        # structData[key]['T'] = []
    # structData[key]['h'].append(float(sim[' h']))
    # structData[key]['elinAlg'].append(float(sim[' linAlgErr']))
    # structData[key]['el2'].append(float(sim[' l2Err']))
    # structData[key]['del2'].append(float(sim[' dL2Err']))
    # structData[key]['T'].append(float(sim[' runtime']))

filename = 'HDG/HDGLaplace.csv'

simData = readcsv(resultsDir + filename)

for sim in simData:
    key = ('HDG', sim['dim'], sim[' order'])
    if key not in structData:
        structData[key] = {}
        structData[key]['h'] = []
        structData[key]['elinAlg'] = []
        structData[key]['el2'] = []
        structData[key]['del2'] = []
        structData[key]['T'] = []
    structData[key]['h'].append(float(sim[' h']))
    structData[key]['elinAlg'].append(float(sim[' linAlgErr']))
    structData[key]['el2'].append(float(sim[' l2Err']))
    structData[key]['del2'].append(float(sim[' dL2Err']))
    structData[key]['T'].append(float(sim[' runtime']))

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

plt.xlabel(r'$\log(h)$')
plt.ylabel(r'$\log(e)$')
plt.xlim((-1.35, -0.0))

errType = 'el2'

for key in structData:
    structData[key]['slope'] = [(np.log10(structData[key][errType][i]) - np.log10(structData[key][errType][i-1]))/(np.log10(structData[key]['h'][i]) - np.log10(structData[key]['h'][i-1])) for i in range(1, len(structData[key][errType]))]

c = 0
for key in sorted(structData):
    plt.scatter(np.log10(structData[key]['h']), np.log10(structData[key][errType]), color=tableau20[c])   
    plt.plot(np.log10(structData[key]['h']), np.log10(structData[key][errType]), color=tableau20[c], label=key[0] + '(p=' + key[2] + ', d=' + key[1]+')')
    for i in range(len(structData[key]['h'])-1):
        plt.text(np.log10(structData[key]['h'][i])-0.05, np.log10(structData[key][errType][i])+0.1, '%3.2f' % (structData[key]['slope'][i]), color=tableau20[c])
    c += 1

plt.legend()

plt.show()

# plt.figure()
# # Remove the plot frame lines. They are unnecessary chartjunk.
# ax = plt.subplot(111)
# ax.spines["top"].set_visible(False)
# ax.spines["bottom"].set_visible(True)
# ax.spines["right"].set_visible(False)
# ax.spines["left"].set_visible(True)

# # Ensure that the axis ticks only show up on the bottom and left of the plot.
# # Ticks on the right and top of the plot are generally unnecessary chartjunk.
# ax.get_xaxis().tick_bottom()
# ax.get_yaxis().tick_left()
# plt.xlabel(r'$\log(h)$')
# plt.ylabel(r'$\log(T)$')
# plt.xlim((-1.35, -0.1))

errType = 'T'

for key in structData:
    structData[key]['slope'] = [(np.log10(structData[key][errType][i]) - np.log10(structData[key][errType][i-1]))/(np.log10(structData[key]['h'][i]) - np.log10(structData[key]['h'][i-1])) for i in range(1, len(structData[key][errType]))]


c = 0
for key in sorted(structData):
    plt.scatter(np.log10(structData[key]['h']), np.log10(structData[key][errType]), color=tableau20[c])   
    plt.plot(np.log10(structData[key]['h']), np.log10(structData[key][errType]), color=tableau20[c], label=key[0]+'(p=' + key[2] + ', d=' + key[1]+')')
    for i in range(len(structData[key]['h'])-1):
        plt.text(np.log10(structData[key]['h'][i])-0.07, np.log10(structData[key][errType][i])-0.05, '%3.2f' % (structData[key]['slope'][i]), color=tableau20[c])
    c += 1

plt.legend()

plt.show()

plt.xlabel(r'$\log(T)$')
plt.ylabel(r'$\log(e)$')


errType = 'el2'

for key in structData:
    structData[key]['slope'] = [(np.log10(structData[key][errType][i]) - np.log10(structData[key][errType][i-1]))/(np.log10(structData[key]['T'][i]) - np.log10(structData[key]['T'][i-1])) for i in range(1, len(structData[key][errType]))]

dim2Data = {key: structData[key] for key in structData if key[1] == '2'}

for key in sorted(dim2Data):
    lineChar = '-'
    if(key[0] == 'HDG'):
        lineChar = '-.'
    plt.scatter(np.log10(structData[key]['T']), np.log10(structData[key][errType]), color=tableau20[int(key[2])])   
    # plt.plot(np.log10(structData[key]['T']), np.log10(structData[key][errType]), label=key[0]+'(p=' + key[2] + ', d=' + key[1]+')')
    plt.plot(np.log10(structData[key]['T']), np.log10(structData[key][errType]), lineChar, label=key[0]+'(p=' + key[2] + ')', color=tableau20[int(key[2])])
    # for i in range(len(structData[key]['h'])-1):
        # plt.text(np.log10(structData[key]['T'][i])+0.1, np.log10(structData[key][errType][i])-0.1, '%3.2f' % (structData[key]['slope'][i]))

plt.legend()

plt.show()
