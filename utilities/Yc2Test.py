# Wed Mar 20 09:02:17 CST 2019
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
figs = (13,10)
fonts1 = 23
fonts2 = 25
linew = 3
plt.rcParams['font.family'] = 'serif'
from matplotlib import rc
plt.rcParams['mathtext.fontset'] = 'stix'
# x_ticks = np.arange(0.0, 0.025, 0.005)

Zst = 0.0633883
numLoop = int(input('Enter loop number:\n')) # Total number of loops
nLoop = range(0,numLoop+1,1)
Yc1St = []
Yc2St = []
Tmax = []
for n in nLoop:
    filename = 'flameletTable_{:}.csv'.format(n)

    data = np.loadtxt(filename, delimiter=',', skiprows=1)
    data = np.transpose(data)
    Z = data[0]
    Yc1 = data[1]
    Yc2 = data[2]
    T = data[3]

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

    Yc1St.append(Yc1[pU] * weightU + Yc1[pL] * weightL)
    Yc2St.append(Yc2[pU] * weightU + Yc2[pL] * weightL)
    #Tmax.append(T[pU] * weightU + T[pL] * weightL)
    Tmax.append(np.max(T))

# Plot T-Z
fig = plt.figure(figsize=figs)
ax = Axes3D(fig)
ax.scatter(Yc1St,Yc2St,Tmax,c='k',marker='o')
#plt.scatter(np.log10(Yc1St),Tmax,c='k',marker='o')
plt.tick_params(labelsize=fonts1)
plt.xlabel(r'$Yc_{1,st} \ (-)$',fontsize=fonts1)
#plt.xlabel(r'$log(\chi_{st}) \ (1/s)$',fontsize=fonts1)
plt.ylabel(r'$Yc_{2,st} \ (-)$',fontsize=fonts1)

plt.show()
