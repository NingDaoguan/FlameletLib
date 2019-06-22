# Thu May 16 14:29:15 CST 2019
import os
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

p = 101325.0

print("Enter file names:")
filename = []
while True:
    name_input = input('> ')
    if name_input == '':
        break
    else:
        filename.append(name_input)

data_directory = 'hatables/'
if not os.path.exists(data_directory):
    os.makedirs(data_directory)

fuelType = input('1:\tJet-A\n2:\tnc12h26\n3:\tC2H5OH\n> ')
if fuelType=='1':
    gas = ct.Solution("KEROSENE_CRECK231.cti")
    fuel = ['NC12H26','IC16H34','DECALIN','C7H8']
    Y_f = np.array([0.29152, 0.465052, 0.19402, 0.04941])
elif fuelType=='2':
    gas = ct.Solution("nDodecane_CRECK.cti")
    fuel = ['NC12H26']
    Y_f = np.array([1.0])
elif fuelType=='3':
    gas = ct.Solution("Ethanol_31.cti")
    fuel = ['C2H5OH']
    Y_f = np.array([1.0])
else:
    raise Exception('Only support the above fuels')
speciesNames = gas.species_names
nsp = len(speciesNames)
molW = gas.molecular_weights
xIndex = 0
TIndex = 1
N2Index = 2
ARIndex = 3
O2Index = 6
H2Index = 9
H2OIndex = 10
COIndex = 12
CO2Index = 13

ax1 = plt.subplot(221)
ax2 = plt.subplot(222)
ax3 = plt.subplot(223)
ax4 = plt.subplot(224)
for n,file in enumerate(filename):

    filename2 = data_directory + 'flameletTable_{:}.csv'.format(n)
    with open(filename2, 'w+') as fo:
        line = 'Z,ha,Yc,omegaYc,T'
        for i in speciesNames:
            line += f',{i}'
        fo.write(line+'\n')

    data1orig = np.loadtxt(file)
    data1 = np.transpose(data1orig)
    T = data1[TIndex]
    YAR = data1[ARIndex]
    YARO = max(YAR[-1], YAR[0])
    Z = (YAR - YARO) / (0.0-YARO)
    # Yc = data1[H2OIndex] + data1[H2Index] + data1[CO2Index] + data1[COIndex]
    Yc = data1[H2OIndex] + data1[CO2Index]

    # Calculate omegaYc
    Y = []
    ha = np.zeros(len(Yc))
    omegaYc = np.zeros(len(Yc))
    for i in range(len(Yc)):
        Y = data1orig[i][2::]
        gas.TPY = (T[i], p, Y)
        ha[i] = gas.enthalpy_mass
        omegaYc[i] = gas.net_production_rates[gas.species_index('CO2')] \
                    * molW[gas.species_index('CO2')] \
                    +gas.net_production_rates[gas.species_index('H2O')] \
                    * molW[gas.species_index('H2O')]


    data2 = []
    data2.append(list(Z))
    data2.append(list(ha))
    data2.append(list(Yc))
    data2.append(list(omegaYc))
    ax1.plot(Z,T)
    ax2.plot(Z,ha)
    ax3.plot(Z,Yc)
    ax4.plot(Z,omegaYc)
    data2.append(list(T))
    for i in range(len(data1)):
        if i >= 2:
            data2.append(list(data1[i]))
        else:
            pass
    data2 = np.array(data2)
    data2 = np.transpose(data2)
    if Z[0] > 0.2:
        data2 = data2[::-1]
    else:
        pass
    with open(filename2,'a') as f:
        np.savetxt(f, data2, delimiter=',',fmt='%f')
plt.savefig('hatables.png',dpi=500)
plt.show()

