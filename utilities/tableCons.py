import numpy as np
import matplotlib.pyplot as plt

numLoop = 13 # Total number of loops
p = 101325.0 # Ambient pressure
Dz0 = 0.2e-4 # Diffusion coefficient of Z @ 273K
filename2 = 'flameletTable.csv'

data2 = []
nLoop = range(0,numLoop+1,1)
for n in nLoop:
    if n == 0:
        filename1 = 'initial_solution.csv'
        # Read names
        with open(filename1) as f:
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
        filename1 = 'strain_loop_{0:02d}.csv'.format(n)

    data1 = np.loadtxt(filename1, delimiter=',', skiprows=1)
    data1 = np.transpose(data1)
    x = data1[0]
    T = data1[3]
    YAR = data1[5]
    YARO = YAR[-1]
    Z = (YAR-YARO) / (0.0-YARO)
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
    plt.plot(Z,chi)

data2 = np.array(data2)
data2 = np.transpose(data2)
with open(filename2,'a') as f:
    np.savetxt(f, data2, delimiter=',')
plt.show()
