# Wed Mar 20 09:02:17 CST 2019
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.size'] = 20
mpl.rcParams['font.weight'] = 'medium'
mpl.rcParams['font.style'] = 'normal'
# mpl.rcParams['font.serif'] = 'DejaVu Serif'
# mpl.rcParams['font.serif'] = 'Georgia'
# mpl.rcParams['font.serif'] = 'Times New Roman'
# mpl.rcParams['text.usetex'] = True
mpl.rcParams['mathtext.fontset'] = 'stix'
# mpl.rcParams['mathtext.fallback_to_cm'] = True
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['savefig.bbox'] = 'tight'

p0 = 101325.0 # Ambient pressure
figs = (13,10)
fonts1 = 23
fonts2 = 25
linew = 3
# x_ticks = np.arange(0.0, 0.025, 0.005)
xin = 'n'
print("Enter file names:")
filename = []
while True:
    name_input = input('> ')
    if name_input == '':
        break
    else:
        filename.append(name_input)

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
    if name[i] == 'z' or name[i] == 'z (m)' or name[i] == 'x (m)':
        xIndex = i
    elif name[i] == 'u' or name[i] == 'u (m/s)':
        uIndex = i
    elif name[i] == 'volFrac (-)':
        volFracIndex = i
        continue
    elif name[i] == 'T' or name[i] == 'T (K)':
        TIndex = i
    elif name[i] == 'O2':
        O2Index = i
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
        majorIndex.append(i)
        continue
    elif name[i] == 'OH':
        OHIndex = i
        majorIndex.append(i)
        continue
    elif name[i] == 'AR':
        ARIndex = i
        majorIndex.append(i)
        continue
    elif name[i] == 'NC12H26':
        NC12H26Index = i
        continue
    elif name[i] == 'IC16H34':
        IC16H34Index = i
        continue
    elif name[i] == 'DECALIN':
        DECALINIndex = i
        continue
    elif name[i] == 'C7H8':
        C7H8Index = i
        continue
    elif name[i] == 'C2H5OH':
        majorIndex.append(i)
        C2H5OHIndex = i
        continue
    # elif name[i] == 'CH4':
    #     majorIndex.append(i)
    #     CH4Index = i
    #     continue

# Plot T, Y, Z and droplets
for j,file in enumerate(filename):
    fig = plt.figure(figsize=figs)
    ax1 = plt.subplot(111)
    ax2 = ax1.twinx()
    data = np.loadtxt(file, delimiter=',', skiprows = 1)
    data = np.transpose(data)
    x = data[xIndex]
    u = data[uIndex]
    a = (u[0] - u[-1]) / (x[-1] - x[0])
    T = data[TIndex]
    YIN = data[ARIndex]
    YIN_O = max(YIN[-1], YIN[0])
    YIN_F = 0
    Z = (YIN - YIN_O)/(YIN_F - YIN_O)
    print(f'Mean strain rate:\t{a}')
    ax1.plot(x,T,marker='^',label='T',c='k',ls='-.',lw=linew,ms=6)
    for k in majorIndex:
        ax2.plot(x,data[k],label=name[k],lw=linew)
    # ax2.plot(x,10000*data[HCOIndex],label=name[HCOIndex] + r'$\times 10^4$',lw=linew)
    # ax2.plot(x,10*data[OHIndex],label=name[OHIndex] + r'$\times 10$',lw=linew)
    # ax2.plot(x,data[volFracIndex],label='Volume fraction',ls=':',c='b',lw=linew)
    ax2.plot(x,Z,label='Z',ls=':',c='r',lw=linew)

    filename2 = 'LAG' + file
    try:
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
            if i%10==0:
                p.append(ip)
                #d.append(d2[i]*d2[i])
                d.append(d2[i])
        p.append(p1[-1])
        d.append(0)
        ax2.scatter(p,d,label=r'Droplet',marker='o',s=50,c='r')
        ax2.set_ylabel(r'$Y$ / $Z$ / $d^*$ (-)')
    except IOError:
        ax2.set_ylabel(r'$Y$ / $Z$ (-)')
        pass
    
    filename2 = 'dispersed' + file[6::]
    try:
        data2 = np.loadtxt(filename2,delimiter=',',comments='#',skiprows=1)
        data2 = np.transpose(data2)
        volFrac = data2[1]
        rhod = data2[2]
        Td = data2[3]
        ax2.plot(data[0],volFrac,label='Volume fraction',ls=':',c='b',lw=linew)
        ax2.set_ylabel(r'$Y$ / $Z$ / $\alpha$ (-)')
    except IOError:
        pass

    ax1.set_xlabel(r'$x$ (m)')
    ax1.set_ylabel(r'$T$ (K)')
    ax1.margins(x=0.0)
    ax1.set_ylim(300,2500)
    # ax1.set_ylim(bottom=300)
    ax2.set_ylim(0,1)
    # ax2.set_ylim(bottom=0)
    fig.legend(bbox_to_anchor=(1,1), bbox_transform=ax1.transAxes)
    x_ticks = np.linspace(x[0], x[-1], 5)
    plt.xticks(x_ticks,color='k')
    fig.tight_layout()
    # plt.savefig(file+'-x-T.png',dpi=500,bbox_inches='tight')


# Plot T-Z
fig, ax = plt.subplots(figsize=figs)
for i,file in enumerate(filename):
    data = np.loadtxt(file, delimiter=',', skiprows = 1)
    data = np.transpose(data)
    x = data[xIndex]
    T = data[TIndex]
    if xin is 'n':
        YIN = data[ARIndex]
        YIN_O = max(YIN[-1], YIN[0])
    else:
        YIN = data[N2Index]
        YIN_O = max(YIN[-1], YIN[0])
    YIN_F = 0
    Z = (YIN - YIN_O)/(YIN_F - YIN_O)
    ax.plot(Z,T,label=file)
ax.set_xlabel(r'$Z$ (-)')
ax.legend()
ax.set_ylim(bottom=300)
# ax.margins(x=0.0)
# ax.set_xlim(0,1.0)
# ax.set_ylim(300,2200)
ax.set_ylabel(r'$T$ (K)')
fig.tight_layout()

# Plot Z-Yc:T
fig, ax = plt.subplots(figsize=figs)
for i,file in enumerate(filename):
    data = np.loadtxt(file, delimiter=',', skiprows = 1)
    data = np.transpose(data)
    x = data[xIndex]
    T = data[TIndex]
    YIN = data[ARIndex]
    YIN_O = max(YIN[-1], YIN[0])
    YIN_F = 0
    Z = (YIN - YIN_O)/(YIN_F - YIN_O)
    #Yc = data[COIndex] + data[H2Index] + data[CO2Index] + data[H2OIndex]
    Yc = data[CO2Index] + data[H2OIndex]
    sc = ax.scatter(Z,Yc,c=T,s=20,vmin=300, vmax=2200)
v = [300,500,1000,1500,2000,2200]
cbar = fig.colorbar(sc,ticks=v)
cbar.set_label(r'$T$ (K)')

# ax.margins(x=0.0)
ax.set_ylim(bottom=0)
ax.set_xlabel(r'$Z \ (-)$')
ax.set_ylabel(r'$Y_c \ (-)$')
fig.tight_layout()
# plt.savefig('Z-Yc-T.png', dpi=500,bbox_inches='tight')

plt.show()