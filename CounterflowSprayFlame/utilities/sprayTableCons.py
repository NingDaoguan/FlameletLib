# Thu Mar  7 10:15:11 CST 2019
import numpy as np
import matplotlib.pyplot as plt

numLoop = int(input('Enter loop number:\n')) # Total number of loops
y = input("1:\tchi\n2:\tYc\n")
p = 101325.0 # Ambient pressure
Dz0 = 0.2e-4 # Diffusion coefficient of Z @ 273K

ax1 = plt.subplot(121)
ax2 = plt.subplot(122)
nLoop = range(0,numLoop+1,1)
for n in nLoop:
    data2 = []
    filename = 'a{:}.txt'.format(n)
    if y == '1':
        filename2 = './tablesChi/sprayTable_{:}.csv'.format(n)
    elif y == '2':
        filename2 = './tablesYc/sprayTable_{:}.csv'.format(n)
    else:
        print('Input Error!')
    # Read names
    with open(filename) as f:
        line1 = f.readline()
        line1 = line1.split(',')
        names = []
        for i in range(len(line1)):
            names.append(line1[i])
        if y == '1':
            line2 = 'Z,chi,T'
        elif y == '2':
            line2 = 'Z,Yc,T'
        else:
            print('Input Error!')
        for i in range(len(names) - 7):
            line2 += ',' + names[i+6]
        with open(filename2, 'w+') as f:
            f.write(line2 + '\n')

    data1 = np.loadtxt(filename, delimiter=',', skiprows=1)
    data1 = np.transpose(data1)
    x = data1[0]
    T = data1[3]
    H2OIndex = 10
    H2Index = 9
    CO2Index = 13
    COIndex = 12
    Yc = data1[H2OIndex] + data1[H2Index] + data1[CO2Index] + data1[COIndex]
    YAR = data1[6]
    YARO = YAR[0]
    Z = (YAR - YARO) / (0.0-YARO)
    chi = np.zeros(len(x)) # Dissipation rate \chi
    chi[0] = (Z[1] - Z[0]) / (x[1] - x[0]) # Gradient computation
    for i in range(len(x) - 2):
        chi[i+1] = (Z[i+2] - Z[i]) / (x[i+2] - x[i])
    chi[len(x) - 1] = (Z[-1] - Z[-2]) / (x[-1] - x[-2])
    chi = 2 * chi**2 # 2(grad(Z))^2
    Dz = Dz0 * (T/273.0)**1.5 * 101325.0/p # Diffusion coefficient Dz
    chi = Dz*chi # chi = 2D(grad(Z))^2

    data2.append(list(Z))
    if y == '1':
        data2.append(list(chi))
        ax1.plot(Z,chi)
        ax2.plot(Z,T)
    elif y == '2':
        data2.append(list(Yc))
        ax1.plot(Z,Yc)
        ax2.plot(Z,T)
    else:
        print('Input Error!')
    data2.append(list(T))
    for i in range(len(data1) - 7):
        data2.append(list(data1[i+6]))
    data2 = np.array(data2)
    data2 = np.transpose(data2)
    with open(filename2,'a') as f:
        np.savetxt(f, data2, delimiter=',')
plt.show()

