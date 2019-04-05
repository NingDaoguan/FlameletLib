# Thu Mar  7 10:15:11 CST 2019
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

p = 101325.0 # Ambient pressure
Dz0 = 0.2e-4 # Diffusion coefficient of Z @ 273K

chiYc = input("1:\tchi\n2:\tYc\n")
if chiYc == '1':
    y = input("Enter loop numer\n")
    numLoop = int(y) # Total number of loops
    x = input("1:\tJet-A\n2:\tnc12h26\n")
    if x=='1':
        gas = ct.Solution("KEROSENE_CRECK231.cti")
        fuel = ['NC12H26','IC16H34','DECALIN','C7H8']
        Y_f = np.array([0.29152, 0.465052, 0.19402, 0.04941])
    elif x=='2':
        gas = ct.Solution("nDodecane_CRECK.cti")
        fuel = ['NC12H26']
        Y_f = np.array([1.0])
    else:
        print('Unknown fuel type!')
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


    ax1 = plt.subplot(121)
    ax2 = plt.subplot(122)
    nLoop = range(0,numLoop+1,1)
    for n in nLoop:
        if n == 0:
            filename = 'initial_solution.csv'
        else:
            filename = 'strain_loop_{0:02d}.csv'.format(n)

        data2 = []
        filename2 = './tablesChi/flameletTable_{:}.csv'.format(n)
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

        data2.append(list(Z[::-1]))
        data2.append(list(chi[::-1]))
        data2.append(list(T[::-1]))
        for i in range(len(data1) - 5):
            data2.append(list(data1[i+5][::-1]))
        data2 = np.array(data2)
        data2 = np.transpose(data2)
        with open(filename2,'a') as f:
            np.savetxt(f, data2, delimiter=',')
        ax1.plot(Z,chi)
        ax2.plot(Z,T)
    plt.show()

elif chiYc == '2':
    y = input("Enter loop numer\n")
    numLoop = int(y) # Total number of loops
    x = input("1:\tJet-A\n2:\tnc12h26\n")
    if x=='1':
        gas = ct.Solution("KEROSENE_CRECK231.cti")
        fuel = ['NC12H26','IC16H34','DECALIN','C7H8']
        Y_f = np.array([0.29152, 0.465052, 0.19402, 0.04941])
    elif x=='2':
        gas = ct.Solution("nDodecane_CRECK.cti")
        fuel = ['NC12H26']
        Y_f = np.array([1.0])
    else:
        print('Unknown fuel type!')
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


    ax1 = plt.subplot(121)
    ax2 = plt.subplot(122)
    nLoop = range(0,numLoop+1,1)
    for n in nLoop:
        if n == 0:
            filename = 'initial_solution.csv'
        else:
            filename = 'strain_loop_{0:02d}.csv'.format(n)

        data2 = []
        filename2 = './tablesYc/flameletTable_{:}.csv'.format(n)
        # Read names
        with open(filename) as f:
            line1 = f.readline()
            line1 = line1.split(',')
            names = []
            for i in range(len(line1)):
                names.append(line1[i])
            line2 = 'Z,Yc,T'
            for i in range(len(names) - 5):
                line2 += ',' + names[i+5]
            with open(filename2, 'w+') as f:
                f.write(line2)

        data1 = np.loadtxt(filename, delimiter=',', skiprows=1)
        data1 = np.transpose(data1)
        x = data1[0]
        T = data1[3]
        H2OIndex = 9
        H2Index = 8
        CO2Index = 12
        COIndex = 11
        Yc = data1[H2OIndex] + data1[H2Index] + data1[CO2Index] + data1[COIndex]
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

        data2.append(list(Z[::-1]))
        data2.append(list(Yc[::-1]))
        data2.append(list(T[::-1]))
        for i in range(len(data1) - 5):
            data2.append(list(data1[i+5][::-1]))
        data2 = np.array(data2)
        data2 = np.transpose(data2)
        with open(filename2,'a') as f:
            np.savetxt(f, data2, delimiter=',')
        ax1.plot(Z,Yc)
        ax2.plot(Z,T)
    plt.show()

else:
    print("Input Error!")