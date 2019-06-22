# Wed Mar 20 09:02:17 CST 2019
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

p0 = 101325.0 # Ambient pressure
figs = (14,9)
fonts1 = 23
fonts2 = 25
linew = 3
plt.rcParams['font.family'] = 'serif'
from matplotlib import rc
plt.rcParams['mathtext.fontset'] = 'stix'
# x_ticks = np.arange(0.0, 0.025, 0.005)
print("Enter file names:")
filename = []
while True:
    name_input = input('> ')
    if name_input == '':
        break
    else:
        filename.append(name_input)

xIndex = 0
TIndex = 1
N2Index = 2
ARIndex = 3
O2Index = 6
H2Index = 9
H2OIndex = 10
COIndex = 12
CO2Index = 13
majorIndex = [2,3,6,9,10,12,13]
name = ['N2','AR','O2','H2','H2O','CO','CO2']


# Plot T, Y, Z and droplets
for j,file in enumerate(filename):
    fig = plt.figure(figsize=figs)
    ax1 = plt.subplot(111)
    ax2 = ax1.twinx()
    data = np.loadtxt(file)
    data = np.transpose(data)
    x = data[xIndex]
    T = data[TIndex]
    YIN = data[ARIndex]
    YIN_O = max(YIN[-1], YIN[0])
    YIN_F = 0
    Z = (YIN - YIN_O)/(YIN_F - YIN_O)
    ax1.plot(x,T,marker='^',label='T',c='k',ls='-.',lw=linew,ms=6)
    for k in range(len(majorIndex)):
        ax2.plot(x,data[majorIndex[k]],label=name[k],lw=linew)
    # ax2.plot(x,10000*data[HCOIndex],label=name[HCOIndex] + r'$\times 10^4$',lw=linew)
    # ax2.plot(x,10*data[OHIndex],label=name[OHIndex] + r'$\times 10$',lw=linew)
    ax2.plot(x,Z,label='Z',ls=':',c='r',lw=linew)
    ax2.set_ylabel(r'$Y$ / $Z$ (-)',fontsize=fonts1)
    ax1.set_xlabel(r'x (m)', fontsize=fonts1)
    ax1.set_ylabel(r'Temperature (K)',fontsize=fonts1)
    ax1.tick_params(labelsize=fonts1)
    ax2.tick_params(labelsize=fonts1)
    # ax2.legend(loc=0,fontsize=fonts1)
    # ax1.legend(loc=0,fontsize=fonts1)
    fig.legend(fontsize=fonts1,bbox_to_anchor=(1,1),bbox_transform=ax1.transAxes)
    # plt.xticks(x_ticks,color='k')
    # plt.savefig(file+'-x-T.png',dpi=500,bbox_inches='tight')


# Plot Z-Yc:T
plt.figure(figsize=figs)
for i,file in enumerate(filename):
    data = np.loadtxt(file)
    data = np.transpose(data)
    x = data[xIndex]
    T = data[TIndex]
    YIN = data[ARIndex]
    YIN_O = max(YIN[-1], YIN[0])
    YIN_F = 0
    Z = (YIN - YIN_O)/(YIN_F - YIN_O)
    #Yc = data[COIndex] + data[H2Index] + data[CO2Index] + data[H2OIndex]
    Yc = data[CO2Index] + data[H2OIndex]
    norm = matplotlib.colors.Normalize(vmin=300, vmax=2500)
    sc = plt.scatter(Z,Yc,c=T,cmap=plt.cm.rainbow,s=20,norm=norm)
v = [300,500,1000,1500,2000,2500]
cbar = plt.colorbar(sc,ticks=v)
cbar.set_label(r'T (K)', fontsize=fonts2)
for t in cbar.ax.get_yticklabels():
    t.set_fontsize(fonts2)
plt.tick_params(labelsize=fonts1)
plt.xlabel(r'$Z \ (-)$',fontsize=fonts1)
plt.ylabel(r'$Y_c \ (-)$',fontsize=fonts1)
# plt.savefig('Z-Yc-T.png', dpi=500,bbox_inches='tight')

plt.show()
