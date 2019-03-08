# Thu Mar  7 10:15:11 CST 2019
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

numLoop = 25 # Total number of loops
p = 101325.0 # Ambient pressure
Dz0 = 0.2e-4 # Diffusion coefficient of Z @ 273K
filename2 = 'flameletTable.csv'

fuel = ['NC12H26','IC16H34','DECALIN','C7H8']
Y_f = np.array([0.29152, 0.465052, 0.19402, 0.04941])
gas = ct.Solution("KEROSENE_CRECK231.cti")
names = gas.species_names
nsp = len(names)
molW = gas.molecular_weights
muC = np.zeros(nsp)
muH = np.zeros(nsp)
muO = np.zeros(nsp)
mwC = 12.011
mwH = 1.0079
mwO = 15.999
for i, specie in enumerate(names):
    CAtoms = gas.n_atoms(specie,'C')
    HAtoms = gas.n_atoms(specie,'H')
    OAtoms = gas.n_atoms(specie,'O')
    muC[i] = mwC*CAtoms/molW[i]
    muH[i] = mwH*HAtoms/molW[i]
    muO[i] = mwO*OAtoms/molW[i]
# Fuel stream
muC_f = np.zeros(len(fuel))
muH_f = np.zeros(len(fuel))
muO_f = np.zeros(len(fuel))
for i, specie in enumerate(fuel):
    CAtoms = gas.n_atoms(specie,'C')
    HAtoms = gas.n_atoms(specie,'H')
    OAtoms = gas.n_atoms(specie,'O')
    k = gas.species_index(specie)
    muC_f[i] = mwC*CAtoms/molW[k]
    muH_f[i] = mwH*HAtoms/molW[k]
    muO_f[i] = mwO*OAtoms/molW[k]
ZC_f  = np.dot(Y_f, muC_f)
ZH_f  = np.dot(Y_f, muH_f)
ZO_f  = np.dot(Y_f, muO_f)


data2 = []
nLoop = range(0,numLoop+1,1)
for n in nLoop:
    if n == 0:
        filename = 'initial_solution.csv'
        # Read names
        with open(filename) as f:
            line1 = f.readline()
            line1 = line1.split(',')
            names = []
            for i in range(len(line1)):
                names.append(line1[i])
            line2 = 'Z,chi,T'
            for i in range(len(names) - 5):
                line2 += ',' + names[i+5]
            with open(filename2, 'w+') as f:
                f.write(line2)
    else:
        filename = 'strain_loop_{0:02d}.csv'.format(n)

    data1 = np.loadtxt(filename, delimiter=',', skiprows=1)
    data1 = np.transpose(data1)
    x = data1[0]
    T = data1[3]

    # Read Y
    Y = data1[5::]
    Y = np.transpose(Y)
    ZC = np.dot(Y, muC)
    ZH = np.dot(Y, muH)
    ZO = np.dot(Y, muO)
    # Oxidizer stream
    ZC_o = ZC[-1]
    ZH_o = ZH[-1]
    ZO_o = ZO[-1]
    Z = ( 2*(ZC-ZC_o)/mwC + 0.5*(ZH-ZH_o)/mwH - (ZO-ZO_o)/mwO ) \
        / ( 2*(ZC_f-ZC_o)/mwC + 0.5*(ZH_f-ZH_o)/mwH - (ZO_f-ZO_o)/mwO )

    chi = np.zeros(len(x)) # Dissipation rate \chi
    chi[0] = (Z[1] - Z[0]) / (x[1] - x[0]) # Gradient computation
    for i in range(len(x) - 2):
        chi[i+1] = (Z[i+2] - Z[i]) / (x[i+2] - x[i])
    chi[len(x) - 1] = (Z[-1] - Z[-2]) / (x[-1] - x[-2])
    chi = 2 * chi**2 # 2(grad(Z))^2
    Dz = Dz0 * (T/273.0)**1.5 * 101325.0/p # Diffusion coefficient Dz
    chi = Dz*chi # chi = 2D(grad(Z))^2
    if n == 0:
        data2.append(list(Z))
        data2.append(list(chi))
        data2.append(list(T))
        for i in range(len(data1) - 5):
            data2.append(list(data1[i+5]))
    else:
        data2[0].extend(list(Z))
        data2[1].extend(list(chi))
        data2[2].extend(list(T))
        for i in range(len(data1) - 5):
            data2[i+3].extend(list(data1[i+5]))
    # plt.plot(Z,chi)
    plt.plot(Z,T)

data2 = np.array(data2)
data2 = np.transpose(data2)
with open(filename2,'a') as f:
    np.savetxt(f, data2, delimiter=',')
plt.show()
