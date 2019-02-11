import numpy as np
import matplotlib.pyplot as plt

numLoop = 13 # Total number of loops
p = 101325.0 # Ambient pressure
Dz0 = 0.2e-4 # Diffusion coefficient of Z @ 273K

nLoop = range(0,numLoop+1,1)
for n in nLoop:
    if n == 0:
        filename1 = 'initial_solution.csv'
        filename2 = 'a00.csv'
    else:
        filename1 = 'strain_loop_{0:02d}.csv'.format(n)
        filename2 = 'a{0:02d}.csv'.format(n)

    # Read names
    with open(filename1) as f:
        line1 = f.readline()
        line1 = line1.split(',')
        names = []
        for i in range(len(line1)):
            names.append(line1[i])

    data1 = np.loadtxt(filename1, delimiter=',', skiprows=1)
    data1 = np.transpose(data1)
    x = data1[0]
    T = data1[3]
    YAR = data1[5]
    YARO = YAR[-1]
    Z = (YAR-YARO) / (0.0-YARO)
    chi = np.zeros(len(x)) # Dissipation rate
    chi[0] = (Z[1] - Z[0]) / (x[1] - x[0])
    for i in range(len(x) - 2):
        chi[i+1] = (Z[i+2] - Z[i]) / (x[i+2] - x[i])
    chi[len(x) - 1] = (Z[-1] - Z[-2]) / (x[-1] - x[-2])
    chi = 2 * chi**2 # 2(grad(Z))^2
    Dz = Dz0 * (T/273.0)**1.5 * 101325.0/p # Diffusion coefficient Dz
    chi = Dz*chi # chi = 2D(grad(Z))^2

    data2 = [Z,chi,T]
    line2 = 'Z,chi,T'
    for i in range(len(data1) - 5):
        data2.append(data1[i+5])
        line2 += ',' + names[i+5]
    data2 = np.transpose(data2)
    np.savetxt(filename2, data2, delimiter=',')

    with open(filename2, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line2 + content)

    plt.plot(Z,chi)
plt.show()