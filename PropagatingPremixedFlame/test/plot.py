import numpy as np
import matplotlib.pyplot as plt

# x = input("BFER?(y/n)")
x = 'n'
filename = ['output.csv']
figs = (13,10)
fonts1 = 23
fonts2 = 25
linew = 3

plt.figure(figsize=figs)
ax1 = plt.subplot(111)
ax2 = ax1.twinx()
for file in filename:
    data = np.loadtxt(file, delimiter=',', skiprows = 1)
    data = np.transpose(data)
    x = data[0]
    u = data[1]
    T = data[3]
    with open(file) as f:
        line = f.readline()
        line = line.split(',')
        name = []
        for i in range(len(line)):
            name.append(line[i])
        if x is 'n':
            majorIndex = [6,7,9,11,12]
            HCOIndex = 67
            OHIndex = 64
            fuelIndex = 61
            H2OIndex = 9
            H2Index = 8
            CO2Index = 12
            COIndex = 11
        else:
            majorIndex = [6,7,8,9,10]

    YAR = data[6]
    YAR_O = 0.01378956
    YAR_F = 0
    Z = (YAR - YAR_O)/(YAR_F - YAR_O)
    Zave = np.sum(Z) / len(Z)
    print('Z-'+file+f':\t{Zave}')
    ax1.plot(x,T,marker='^',label='T',c='k',ls='-.',lw=linew,ms=6)
    for k in majorIndex:
        ax2.plot(x,data[k],label=name[k],lw=linew)
ax1.set_xlabel(r'x (m)', fontsize=fonts1)
ax1.set_ylabel(r'T (K)',fontsize=fonts1)
ax2.set_ylabel('Mass fractions (-)',fontsize=fonts1)
ax1.tick_params(labelsize=fonts1)
ax2.tick_params(labelsize=fonts1)
ax2.legend(loc=0,fontsize=fonts1)
ax1.legend(loc=0,fontsize=fonts1)
# plt.savefig('x-T.png',dpi=500)
plt.show()