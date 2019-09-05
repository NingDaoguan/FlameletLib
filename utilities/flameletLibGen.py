# Thu May 16 14:29:15 CST 2019
# xu-zhang
import os
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

p = 101325.0  # constant pressure
ctifile = input(
        'CTI FILE:  '
        '  Ethanol_31.cti'
        '  Heptane0.cti\n')
gas = ct.Solution(ctifile)  # ethanol
data_directory = 'tables'  # output dir
if not os.path.exists(data_directory):
    os.makedirs(data_directory)
#=============================================================================#


speciesNames = gas.species_names  # each species
nsp = len(speciesNames)  # number of species
molW = gas.molecular_weights  # molecular weights of each species

# first line
names = []
xIndex = -1
uIndex = -1
TIndex = -1
COIndex = -1
CO2Index = -1
H2Index = -1
H2OIndex = -1
ARIndex = -1
for filename in os.listdir():
    if filename.endswith('.csv'):
        with open(filename) as fi:
            line0 = fi.readline()
            names = [x.strip() for x in line0.split(',')]
        for i in range(len(names)):
            if names[i] == 'z' or names[i] == 'z (m)' or names[i] == 'x' or names[i] == 'x (m)':
                xIndex = i
            elif names[i] == 'u' or names[i] == 'u (m/s)':
                uIndex = i
            elif names[i] == 'T' or names[i] == 'T (K)':
                TIndex = i
            elif names[i] == 'CO':
                COIndex = i
                continue
            elif names[i] == 'CO2':
                CO2Index = i
                continue
            elif names[i] == 'H2':
                H2Index = i
                continue
            elif names[i] == 'H2O':
                H2OIndex = i
                continue
            elif names[i] == 'AR':
                ARIndex = i
                continue
        break  # only read one file

fig = plt.figure()
ax1 = plt.subplot(221)
ax2 = plt.subplot(222)
ax3 = plt.subplot(223)
ax4 = plt.subplot(224)
n = 0  # count
for filename in os.listdir():
    if filename.endswith('.csv'):
        filename2 = os.path.join(data_directory, 'flameletTable_{:}.csv'.format(n))
        n += 1
        with open(filename2, 'w+') as fo:
            line = 'Z,Yc,h,omegaYc,T'
            for i in speciesNames:
                line += f',{i}'
            fo.write(line+'\n')

        dataReaction = np.loadtxt(filename[:-4]+'-reaction', delimiter=',', skiprows=1).T  # reaction source
        data1orig = np.loadtxt(filename, delimiter=',', skiprows=1)
        data1 = np.transpose(data1orig)

        T = data1[TIndex]
        YAR = data1[ARIndex]
        YARO = max(YAR[-1], YAR[0])
        Z = (YAR - YARO) / (0.0-YARO)
        # Yc = data1[H2OIndex] + data1[H2Index] + data1[CO2Index] + data1[COIndex]
        Yc = data1[H2OIndex] + data1[CO2Index]

        for i in range(len(data1)):
            if names[i] == speciesNames[0]:
                speciesStart = i

        # Calculate omegaYc
        omegaYc = dataReaction[CO2Index - speciesStart + 2] \
                 +dataReaction[H2OIndex - speciesStart + 2]
        Y = []
        h = np.zeros(len(Yc))
        for i in range(len(Yc)):
            Y = data1orig[i][speciesStart::]
            gas.TPY = (T[i], p, Y)
            h[i] = gas.enthalpy_mass
            # omegaYc[i] = gas.net_production_rates[CO2Index - speciesStart] * molW[CO2Index - speciesStart] \
            #           +gas.net_production_rates[H2OIndex - speciesStart] * molW[H2OIndex - speciesStart]

        data2 = []
        data2.append(list(Z))
        data2.append(list(Yc))
        data2.append(list(h))
        data2.append(list(omegaYc))
        ax1.plot(Z,T)
        ax2.plot(Z,h)
        ax3.plot(Z,Yc)
        ax4.plot(Z,omegaYc)
        data2.append(list(T))
        for i in range(len(data1)):
            if i >= speciesStart:
                data2.append(list(data1[i]))
            else:
                pass
        data2 = np.array(data2)
        data2 = np.transpose(data2)
        with open(filename2, 'a') as f:
            np.savetxt(f, data2, delimiter=',',fmt='%f')
plt.savefig('lib.png',dpi=500)
plt.show()
