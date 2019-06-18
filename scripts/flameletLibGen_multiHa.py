# Thu May 16 14:29:15 CST 2019
import os
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

p = 101325.0

data_directory_orig = 'multihatables_frozenOmega/ha_orig/'
if not os.path.exists(data_directory_orig):
    os.makedirs(data_directory_orig)

totalHaNum = 6

numLoop = int( input('Enter loop numer\n> ') )
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

# Read first line
with open('initial_solution.csv') as fi:
    line = fi.readline()
    line = line.split(',')
    names = []
    for i in range(len(line)):
        names.append(line[i])
for i in range(len(names)):
    if names[i] == 'z' or names[i] == 'z (m)':
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

ax1 = plt.subplot(221)
ax2 = plt.subplot(222)
ax3 = plt.subplot(223)
ax4 = plt.subplot(224)
for n in range(0,numLoop+1,1):
    if n == 0:
        filename = 'initial_solution.csv'
    else:
        filename = 'strain_loop_{0:02d}.csv'.format(n)

    filename_orig = data_directory_orig + 'flameletTable_{:}.csv'.format(n)


    with open(filename_orig, 'w+') as fo:
        line = 'Z,ha,Yc,omegaYc,T'
        for i in speciesNames:
            line += f',{i}'
        fo.write(line+'\n')

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
    Y = []
    ha = np.zeros(len(Yc))
    omegaYc = np.zeros(len(Yc))
    for i in range(len(Yc)):
        Y = data1orig[i][speciesStart::]
        gas.TPY = (T[i], p, Y)
        ha[i] = gas.enthalpy_mass
        omegaYc[i] = gas.net_production_rates[CO2Index - speciesStart] * molW[CO2Index - speciesStart] \
                    +gas.net_production_rates[H2OIndex - speciesStart] * molW[H2OIndex - speciesStart]

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
        if i >= speciesStart:
            data2.append(list(data1[i]))
        else:
            pass
    data2 = np.array(data2)
    if Z[0] > 0.2:
        dataOut = np.copy(data2.T[::-1])
    else:
        dataOut = np.copy(data2.T)
    with open(filename_orig,'a') as f:
        np.savetxt(f, dataOut, delimiter=',',fmt='%f')



    # low ha tables
    haorig = np.copy(ha)
    dha = [0.1e05, 0.2e05, 0.5e05, 1e05, 1.5e05, 2e05]
    for hi in range(totalHaNum):
        data_directory_low = 'multihatables_frozenOmega/ha_{:}/'.format(totalHaNum-hi-1)
        if not os.path.exists(data_directory_low):
            os.makedirs(data_directory_low)
        filename_low = data_directory_low + 'flameletTable_{:}.csv'.format(n)

        with open(filename_low, 'w+') as fo:
            line = 'Z,ha,Yc,omegaYc,T'
            for i in speciesNames:
                line += f',{i}'
            fo.write(line+'\n')


        ha = haorig - dha[hi]
        # Calculate omegaYc
        Y = []
        T = np.zeros(len(Z))
        # omegaYc = np.zeros(len(Z))
        for i in range(len(Z)):
            Y = data1orig[i][speciesStart::]
            gas.HPY = (ha[i], p, Y)
            T[i] = gas.T
            # omegaYc[i] = gas.net_production_rates[CO2Index - speciesStart] * molW[CO2Index - speciesStart] \
            #             +gas.net_production_rates[H2OIndex - speciesStart] * molW[H2OIndex - speciesStart]

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
            if i >= speciesStart:
                data2.append(list(data1[i]))
            else:
                pass
        data2 = np.array(data2)
        if Z[0] > 0.2:
            dataOut = np.copy(data2.T[::-1])
        else:
            dataOut = np.copy(data2.T)
        with open(filename_low,'a') as f:
            np.savetxt(f, dataOut, delimiter=',',fmt='%f')

plt.savefig('multihatables_frozenOmega.png',dpi=500)
plt.show()
