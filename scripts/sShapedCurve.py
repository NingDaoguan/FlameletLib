# Wed Mar 20 09:02:17 CST 2019
import numpy as np
import matplotlib.pyplot as plt
figs = (13,10)
fonts1 = 23
fonts2 = 25
linew = 3
plt.rcParams['font.family'] = 'serif'
from matplotlib import rc
plt.rcParams['mathtext.fontset'] = 'stix'
# x_ticks = np.arange(0.0, 0.025, 0.005)
Zst = 0
fuel = input('Fuel:\n')
if fuel == 'Jet-A':
    Zst = 0.063388
elif fuel == 'Ethanol':
    Zst = 0.1002
else:
    print('Valid fuels are `Jet-A` and `Ethanol`')
stOrMax = input('`Tmax` or `Tst`?\n')
chiYc = input('`chi` or `Yc`?\n')
numLoop = int(input('Enter loop number:\n')) # Total number of loops
nLoop = range(0,numLoop+1,1)
chiSt = []
Tmax = []
for n in nLoop:
    filename = 'flameletTable_{:}.csv'.format(n)

    data = np.loadtxt(filename, delimiter=',', skiprows=1)
    data = np.transpose(data)
    Z = data[0]
    chi = data[1]
    T = data[2]

    for j,iZ in enumerate(Z):
        if iZ > Zst:
            Zub = iZ
            Zlb = Z[j-1]
            weightU = (Zst - Zlb) / (Zub - Zlb)
            weightL = (Zub - Zst) / (Zub - Zlb)
            pU = j
            pL = j-1
            break
        else:
            continue

    chiSt.append(chi[pU] * weightU + chi[pL] * weightL)
    if stOrMax == 'Tst':
        Tmax.append(T[pU] * weightU + T[pL] * weightL)
    elif stOrMax == 'Tmax':
        Tmax.append(np.max(T))

# Plot T-Z
plt.figure(figsize=figs)
plt.scatter(chiSt,Tmax,c='k',marker='o')
#plt.scatter(np.log10(chiSt),Tmax,c='k',marker='o')
plt.tick_params(labelsize=fonts1)
if chiYc == 'chi':
    plt.xlabel(r'$\chi_{st} \ (1/s)$',fontsize=fonts1)
    #plt.xlabel(r'$log(\chi_{st}) \ (1/s)$',fontsize=fonts1)
elif chiYc == 'Yc':
    plt.xlabel(r'$Y_{c,st} \ (-)$',fontsize=fonts1)
if stOrMax == 'Tst':
    plt.ylabel(r'$T_{st} \ (K)$',fontsize=fonts1)
elif stOrMax == 'Tmax':
    plt.ylabel(r'$T_{max} \ (K)$',fontsize=fonts1)

plt.show()
