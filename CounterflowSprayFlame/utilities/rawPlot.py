# Wed Mar 20 09:02:17 CST 2019
import numpy as np
import matplotlib.pyplot as plt

filename = ['mf-0.3mo-0.3L-0.02_raw.txt']
xin = 'n'
figs = (13,10)
fonts1 = 23
fonts2 = 25
linew = 3
plt.rcParams['font.family'] = 'serif'
from matplotlib import rc
plt.rcParams['mathtext.fontset'] = 'stix'
x_ticks = np.arange(0.0, 0.025, 0.005)

# Read names
file = filename[0]
with open(file) as f:
    line = f.readline()
    line = line.split(',')
    name = []
    for i in range(len(line)):
        name.append(line[i])
majorIndex = []
for i in range(len(name)):
    if name[i] == 'O2':
        majorIndex.append(i)
        continue
    elif name[i] == 'N2':
        N2Index = i
        majorIndex.append(i)
        continue
    elif name[i] == 'CO':
        COIndex = i
        majorIndex.append(i)
        continue
    elif name[i] == 'CO2':
        CO2Index = i
        majorIndex.append(i)
        continue
    elif name[i] == 'H2':
        H2Index = i
        majorIndex.append(i)
        continue
    elif name[i] == 'H2O':
        H2OIndex = i
        majorIndex.append(i)
        continue
    elif name[i] == 'HCO':
        HCOIndex = i
        continue
    elif name[i] == 'OH':
        OHIndex = i
        continue
    elif name[i] == 'AR':
        ARIndex = i
        continue

# Plot T, Y and Z
for j,file in enumerate(filename):
    plt.figure(figsize=figs)
    ax1 = plt.subplot(111)
    ax2 = ax1.twinx()
    data = np.loadtxt(file, delimiter=',', skiprows = 1)
    data = np.transpose(data)
    x = data[0]
    u = data[1]
    a = (u[0] - u[-1]) / (x[-1] - x[0])
    T = data[3]
    if xin is 'n':
        YIN = data[ARIndex]
        YIN_O = YIN[0]
    else:
        YIN = data[N2Index]
        YIN_O = YIN[0]
    YIN_F = 0
    Z = (YIN - YIN_O)/(YIN_F - YIN_O)
    Zave = np.sum(Z) / len(Z)
    print('Z-'+file+f':\t{Zave}')
    print(f'Strain rate:\t{a}')
    ax1.plot(x,T,marker='^',label='T',c='k',ls='-.',lw=linew,ms=6)
    for k in majorIndex:
        ax2.plot(x,data[k],label=name[k],lw=linew)
    # ax2.plot(x,10000*data[HCOIndex],label=name[HCOIndex] + r'$\times 10^4$',lw=linew)
    # ax2.plot(x,10*data[OHIndex],label=name[OHIndex] + r'$\times 10$',lw=linew)
    ax2.plot(x,Z,label='Z',ls=':',c='r',lw=linew)
    ax1.set_xlabel(r'x (m)', fontsize=fonts1)
    ax1.set_ylabel(r'Temperature (K)',fontsize=fonts1)
    ax2.set_ylabel('Mass fractions / mixture fraction (-)',fontsize=fonts1)
    ax1.tick_params(labelsize=fonts1)
    ax2.tick_params(labelsize=fonts1)
    ax2.legend(loc=0,fontsize=fonts1)
    ax1.legend(loc=0,fontsize=fonts1)
    plt.xticks(x_ticks,color='k')
    plt.savefig(f'{j}'+'-x-T.png',dpi=500)


# Plot T-Z
plt.figure(figsize=figs)
for i,file in enumerate(filename):
    data = np.loadtxt(file, delimiter=',', skiprows = 1)
    data = np.transpose(data)
    x = data[0]
    T = data[3]
    if xin is 'n':
        YIN = data[ARIndex]
        YIN_O = YIN[0]
    else:
        YIN = data[N2Index]
        YIN_O = YIN[0]
    YIN_F = 0
    Z = (YIN - YIN_O)/(YIN_F - YIN_O)
    plt.plot(Z,T,c='k',lw=linew,label = f'{i}')
plt.tick_params(labelsize=fonts1)
plt.xlabel(r'Z (-)',fontsize=fonts1)
plt.ylabel(r'T (K)',fontsize=fonts1)
plt.legend(loc=0,fontsize=fonts1)
plt.savefig('Z-T.png', dpi=500)


# Plot droplet with T_gas
for ic,filename1 in enumerate(filename):
    plt.figure(figsize=figs)
    ax1 = plt.subplot(111)
    ax2 = ax1.twinx()
    filename2 = 'LAG' + filename1
    data1 = np.loadtxt(filename1,delimiter=',',skiprows=1)
    data1 = np.transpose(data1)
    grid1 = data1[0]
    ax1.plot(grid1,data1[3],label=r'$T_{gas}$',c='k',lw=linew)
    data2 = np.loadtxt(filename2,delimiter=',',comments='#',skiprows=0)
    data2 = np.transpose(data2)
    #time
    t1 = data2[0]
    p1 = data2[1]
    d1 = data2[3]
    d2 = d1/d1[0]
    p = []
    d = []
    for i,ip in enumerate(p1):
        if i%3==0:
            p.append(ip)
            #d.append(d2[i]*d2[i])
            d.append(d2[i])
    p.append(p1[-1])
    d.append(0)
    plt.scatter(p,d,label=r'Droplet',marker='o',s=50,c='r')
    
    ax1.set_xlabel(r'x (m)', fontsize=fonts1)
    ax1.set_ylabel(r'Temperature (K)',fontsize=fonts1)
    #ax2.set_ylabel(r'$\left(\frac{d}{d_0}\right)^2$ (-)',fontsize=fonts1)
    ax2.set_ylabel(r'$\frac{d}{d_0}$ (-)',fontsize=fonts1)
    ax1.tick_params(labelsize=fonts1)
    ax2.tick_params(labelsize=fonts1)
    ax2.legend(loc='upper right',fontsize=fonts1)
    ax1.legend(loc='center right',fontsize=fonts1)
    plt.xticks(x_ticks,color='k')
    plt.savefig(f'{ic}'+'-dd0-T.png',dpi=500)

plt.show()
