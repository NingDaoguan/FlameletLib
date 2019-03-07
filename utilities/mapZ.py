# Thu Mar  7 10:15:11 CST 2019
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
fonts1 = 20
fonts2 = 25
linew = 3
figs = (13,10)
Tmax = 2700
x_ticks = np.arange(0.0, 0.025, 0.005)
cs = ['r','c','b','m']
lss = ['-','--','-.',':']
mks = ['o','s','^','p']

fuel = ['NC12H26','IC16H34','DECALIN','C7H8']
Y_f = np.array([0.322, 0.446, 0.185, 0.047])
gas = ct.Solution("KEROSENE_CRECK231.cti")


names = gas.species_names
nsp = len(names)
molW = gas.molecular_weights
muC = np.zeros(nsp)
muH = np.zeros(nsp)
muO = np.zeros(nsp)
mwC = 12.017
mwH = 1.0079
mwO = 15.999
for i, specie in enumerate(names):
    CAtoms = gas.n_atoms(specie,'C')
    HAtoms = gas.n_atoms(specie,'H')
    OAtoms = gas.n_atoms(specie,'O')
    muC[i] = mwC*CAtoms/molW[i]
    muH[i] = mwH*HAtoms/molW[i]
    muO[i] = mwO*OAtoms/molW[i]

# Fuel stream
muC_f = np.zeros(len(fuel))
muH_f = np.zeros(len(fuel))
muO_f = np.zeros(len(fuel))
for i, specie in enumerate(fuel):
    CAtoms = gas.n_atoms(specie,'C')
    HAtoms = gas.n_atoms(specie,'H')
    OAtoms = gas.n_atoms(specie,'O')
    k = gas.species_index(specie)
    muC_f[i] = mwC*CAtoms/molW[k]
    muH_f[i] = mwH*HAtoms/molW[k]
    muO_f[i] = mwO*OAtoms/molW[k]
ZC_f  = np.dot(Y_f, muC_f)
ZH_f  = np.dot(Y_f, muH_f)
ZO_f  = np.dot(Y_f, muO_f)

# Read Y
filename = 'output.csv'
data = np.loadtxt(filename,delimiter=',',skiprows=1)
data = np.transpose(data)
Y = data[6:-1]
Y = np.transpose(Y)
ZC = np.dot(Y, muC)
ZH = np.dot(Y, muH)
ZO = np.dot(Y, muO)

# Oxidizer stream
ZC_o = ZC[-1]
ZH_o = ZH[-1]
ZO_o = ZO[-1]

Z = ( 2*(ZC-ZC_o)/mwC + 0.5*(ZH-ZH_o)/mwH - (ZO-ZO_o)/mwO ) \
    / ( 2*(ZC_f-ZC_o)/mwC + 0.5*(ZH_f-ZH_o)/mwH - (ZO_f-ZO_o)/mwO )

x = data[0]
YAR = data[6]
YARO = YAR[-1]
Z1 = (YAR - YARO) / (0.0-YARO)
T = data[3]

plt.figure(figsize=figs)
plt.plot(Z1,T,label = 'Z based on AR',c='k',lw=linew)
plt.plot(Z,T,label = 'Bilger',c='r',lw=linew)
plt.tick_params(labelsize=fonts1)
plt.xlabel(r'Z (-)',fontsize=fonts1)
plt.ylabel(r'T (K)',fontsize=fonts1)
plt.legend(loc=0,fontsize=fonts1)
plt.savefig('Z-T.png', dpi=500)
plt.show()

'''
filename = ['output.csv']
fonts1 = 20
fonts2 = 25
linew = 3
figs = (13,10)
Tmax = 2700
x_ticks = np.arange(0.0, 0.025, 0.005)
cs = ['r','c','b','m']
lss = ['-','--','-.',':']
mks = ['o','s','^','p']

# plt.figure(figsize=figs)
for filename1 in filename:
    data1 = np.loadtxt(filename1,delimiter=',',skiprows=1)
    data1 = np.transpose(data1)
    #f.grid
    grid1 = data1[0]
    #f.T
    T1 = data1[1]
    #plot T
    # plt.plot(grid1,T1,label = fr'${filename1}$',lw=linew)
    print(filename1)
    a = (T1[0] - T1[-1]) / (grid1[-1]-grid1[0])
    print(f'Mean Strain Rate:\t{a}\n')
# plt.tick_params(labelsize=fonts1)
# plt.title('Velocity Field',fontsize=fonts1)
# plt.xlim(0.0,0.02)
# plt.xticks(x_ticks,color='k')
# plt.xlabel(r'x (m)',fontsize=fonts1)
# plt.ylabel(r'u (m/s)',fontsize=fonts1)
# plt.legend(loc=0,fontsize=fonts1)
# plt.savefig('VelocityField.png', dpi=500)


for ic,filename1 in enumerate(filename):
    plt.figure(figsize=figs)
    data1 = np.loadtxt(filename1,delimiter=',',skiprows=1)
    data1 = np.transpose(data1)
    #f.grid
    grid1 = data1[0]
    plt.plot(grid1,data1[3],label=r'$T$',ls='-',c='k',lw=linew)
    plt.tick_params(labelsize=fonts1)
    plt.xlabel(r'x (m)',fontsize=fonts1)
    plt.ylabel(r'T (K)',fontsize=fonts1)
    plt.ylim(500,Tmax)
    plt.legend(loc='upper right',fontsize=fonts1)
    plt.twinx()
    with open(filename1) as f:
        line = f.readline()
        line = line.split(',')
        name = []
        for i in range(len(line)):
            name.append(line[i])
        majorIndex = [7,8,10,12,13]
        HCOIndex = 68
        OHIndex = 65
        fuelIndex = 62
        H2OIndex = 10
        H2Index = 9
        CO2Index = 13
        COIndex = 12
    data1 = np.loadtxt(filename1,delimiter=',',skiprows=1)
    data1 = np.transpose(data1)
    nsp = data1.shape[0]-1
    Yc = data1[H2OIndex] + data1[H2Index] + data1[CO2Index] + data1[COIndex]
    # Yc = data1[H2OIndex]/0.018 + data1[H2Index]/0.002 + data1[CO2Index]/0.044 + data1[COIndex]/0.028
    # Ycb = 0.081/0.018 + 0.0006/0.002 + 0.1533/0.044 + 0.04/0.028
    # Yc = Yc/Ycb
    plt.plot(grid1,Yc,label = r'$Y_c$',ls='--',lw=linew, c='k')
    plt.tick_params(labelsize=fonts1)
    # plt.title('Temperature and Progress Variable Field @ ' + fr'${filename1}$',fontsize=fonts1)
    plt.xlim(0.0,0.02)
    # plt.ylim(0,1)
    plt.xticks(x_ticks,color='k')
    plt.xlabel(r'x (m)',fontsize=fonts1)
    plt.ylabel(r'$Y_c$ (-)',fontsize=fonts1)
    plt.legend(loc='center right',fontsize=fonts1)
    plt.savefig('SpeciesTemperature'+filename1+'.png', dpi=500)
    # plt.show()



plt.figure(figsize=(13,10))
for ic,filename1 in enumerate(filename):
    filename2 = filenameLAGR[ic]
    data1 = np.loadtxt(filename1,delimiter=',',skiprows=1)
    data1 = np.transpose(data1)

    grid1 = data1[0]

    plt.plot(grid1,data1[3],label=r'$T_{gas}$ @ '+fr'${filename1}$',c=cs[ic],lw=linew)
plt.tick_params(labelsize=fonts2)
plt.xlabel(r'x (m)',fontsize=fonts2)
plt.ylabel(r'T (K)',fontsize=fonts2)
plt.ylim(500,Tmax)
plt.legend(loc='upper right',fontsize=fonts2)
plt.twinx()
for ic,filename1 in enumerate(filename):
    filename2 = filenameLAGR[ic]

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
            d.append(d2[i])
    p.append(p1[-1])
    d.append(0)
    plt.scatter(p,d,label=r'$d/d_0$ @ '+fr'${filename1}$',marker=mks[ic],s=60,c=cs[ic])
    
plt.tick_params(labelsize=fonts2)
# plt.title('Temperature Field and Droplet Diameter @ ' + fr'${filename1}$',fontsize=fonts2)
plt.xlim(0.0,0.02)
plt.ylim(0,1.02)
plt.xticks(x_ticks,color='k')
plt.xlabel(r'x (m)',fontsize=fonts2)
plt.ylabel(r'$d/d_0$ (-)',fontsize=fonts2)
plt.legend(loc='center right',fontsize=fonts2)
plt.savefig('DropletDiameterTemperature.png', dpi=500)
# plt.show()



plt.figure(figsize=figs)
for ic,filename1 in enumerate(filename):
    data1 = np.loadtxt(filename1,delimiter=',',skiprows=1)
    data1 = np.transpose(data1)
    # f.grid
    grid1 = data1[0]
    # YAR
    YAR = data1[6]
    YARO = YAR[-1]

    Z1 = (YAR - YARO) / (0.0-YARO)
    plt.plot(grid1,Z1,label = fr'$Z$ @ ${filename1}$',c='k',ls=lss[ic],lw=linew)
plt.tick_params(labelsize=fonts1)
# plt.title('Mixture Fraction',fontsize=fonts1)
plt.xticks(x_ticks,color='k')
plt.xlabel(r'x (m)',fontsize=fonts1)
plt.ylabel(r'Z (-)',fontsize=fonts1)
plt.xlim(0,0.02)
# plt.ylim(0,1)
plt.legend(loc=0,fontsize=fonts1)
plt.savefig('MixtureFraction.png', dpi=500)
# plt.show()



plt.figure(figsize=figs)
for ic,filename1 in enumerate(filename):
    data1 = np.loadtxt(filename1,delimiter=',',skiprows=1)
    data1 = np.transpose(data1)
    # f.grid
    grid1 = data1[0]
    # YAR
    YAR = data1[6]
    YARO = YAR[-1]
    Z1 = (YAR - YARO) / (0.0-YARO)
    plt.plot(Z1,data1[3],label=fr'${filename1}$',ls=lss[ic],c=cs[ic],lw=linew)
    # plt.scatter(Z1,data1[3],marker=mks[ic],label=fr'${filename1}$',c=cs[ic],s=40)
plt.tick_params(labelsize=fonts1)
plt.xlabel(r'Z (-)',fontsize=fonts1)
plt.ylabel(r'T (K)',fontsize=fonts1)
plt.ylim(500,Tmax)
plt.legend(loc='lower right',fontsize=fonts1)
plt.savefig('ZT.png', dpi=500)



plt.figure(figsize=figs)
for ic,filename1 in enumerate(filename):
    with open(filename1) as f:
        line = f.readline()
        line = line.split(',')
        name = []
        for i in range(len(line)):
            name.append(line[i])
    H2OIndex = 10
    H2Index = 9
    CO2Index = 13
    COIndex = 12

    data1 = np.loadtxt(filename1,delimiter=',',skiprows=1)
    data1 = np.transpose(data1)
    # f.grid
    grid1 = data1[0]
    # YAR
    YAR = data1[6]
    YARO = YAR[-1]
    Z1 = (YAR - YARO) / (0.0-YARO)
    Yc = data1[H2OIndex] + data1[H2Index] + data1[CO2Index] + data1[COIndex]
    plt.plot(Z1,Yc,label = fr'${filename1}$',ls=lss[ic],c=cs[ic],lw=linew)

    plt.tick_params(labelsize=fonts1)
    plt.xlabel(r'Z (-)',fontsize=fonts1)
    plt.ylabel(r'$Y_c$ (-)',fontsize=fonts1)
    plt.legend(loc='upper left',fontsize=fonts1)
plt.savefig('ZC.png', dpi=500)
# plt.show()
'''