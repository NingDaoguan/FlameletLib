# Wed Mar 20 09:02:17 CST 2019
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Configurations
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
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['savefig.bbox'] = 'tight'
figs = (13,9)  # figure size
# x_ticks = np.arange(0.0, 0.025, 0.005)

# Input file names
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


# ------Plot 1------
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
#    ax1.plot(x,T,marker='^',label='T',c='k',ls='-',ms=6)
    ax1.plot(x,T,label='T',c='k',ls='-')
    for k in majorIndex:
        ax2.plot(x,data[k],label=name[k])
    # ax2.plot(x,10000*data[HCOIndex],label=name[HCOIndex] + r'$\times 10^4$')
    # ax2.plot(x,10*data[OHIndex],label=name[OHIndex] + r'$\times 10$')
    # ax2.plot(x,data[volFracIndex],label='Volume fraction',ls=':',c='b')
    ax2.plot(x,Z,label='Z',ls=':',c='r')

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
        nd = data2[2]
        rhod = data2[3]
        Td = data2[4]
        dd = data2[5]
        ax2.plot(data[0],1000*volFrac,label=r'Volume fraction$\times 10^3$',ls=':',c='b')
        ax2.scatter(data[0], dd/dd[0], label=r'Diameter-Norm',marker='o',s=50,c='b')
        ax2.set_ylabel(r'$Y$ / $Z$ / $\alpha$ / $d^*$ (-)')
    except IOError:
        pass

    ax1.set_xlabel(r'$x$ (m)')
    ax1.set_ylabel(r'$T$ (K)')
    ax1.margins(x=0.0)
#    ax1.set_ylim(300,2500)
#    ax1.set_ylim(300,2500)
#    ax1.set_ylim(bottom=300)
    ax1.set_ylim(280, 2200)
    ax2.set_ylim(0,1)
    # ax2.set_ylim(bottom=0)
    fig.legend(bbox_to_anchor=(1,1), bbox_transform=ax1.transAxes)
    x_ticks = np.linspace(x[0], x[-1], 5)
    plt.xticks(x_ticks,color='k')
    fig.tight_layout()
    # plt.savefig(file+'-x-T.png',dpi=500,bbox_inches='tight')


# ------Plot 2------
# Plot T-Z
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
    ax.plot(Z,T,label=file)
ax.set_xlabel(r'$Z$ (-)')
ax.legend()
#ax.set_ylim(bottom=300)
# ax.margins(x=0.0)
# ax.set_xlim(0,1.0)
# ax.set_ylim(300,2200)
ax.set_ylabel(r'$T$ (K)')
fig.tight_layout()

# ------Plot 3------
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

## ------Plot 4------
## Plot Td
#fig, ax = plt.subplots(figsize=figs)
#for i, file in enumerate(filename):
#    filename2 = 'dispersed' + file[6::]
#    try:
#        data2 = np.loadtxt(filename2,delimiter=',',comments='#',skiprows=1)
#        data2 = np.transpose(data2)
#        x = []
#        Td = []
#        for j in range(len(data2[0])):
#            if data2[4][j] > 100:
#                x.append(data2[0][j])
#                Td.append(data2[4][j])
#        ax.scatter(x, Td, marker='o',s=50)
#    except IOError:
#        pass
#ax.set_xlim(0,0.02)
#ax.set_xlabel(r'$x$ (m)')
#ax.set_ylabel(r'$T_d$ (K)')
#ax.legend()
#fig.tight_layout()

# ------Plot 5------
# Plot u
fig, ax = plt.subplots(figsize=figs)
for i, file in enumerate(filename):
    data2 = np.loadtxt(file,delimiter=',',comments='#',skiprows=1)
    data2 = np.transpose(data2)
    x = data2[0]
    u = data2[1]
    ax.plot(x,u,label=file)
ax.set_xlim(0,0.02)
ax.set_xlabel(r'$x$ (m)')
ax.set_ylabel(r'$u$ (m/s)')
ax.legend()
fig.tight_layout()

plt.show()
