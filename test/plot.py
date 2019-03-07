#Sat Dec 15 15:50:09 CST 2018
import numpy as np
import matplotlib.pyplot as plt

filename = ['./diffusion_flame_batch_data/initial_solution.csv']
fonts1 = 20
fonts2 = 25
linew = 3
figs = (13,10)
cs = ['r','c','b','m']
lss = ['-','--','-.',':']
mks = ['o','s','^','p']

plt.figure(figsize=figs)
for filename1 in filename:
    data1 = np.loadtxt(filename1,delimiter=',',skiprows=1)
    data1 = np.transpose(data1)
    #f.grid
    grid1 = data1[0]
    #f.T
    T1 = data1[1]
    #plot T
    plt.plot(grid1,T1,label = fr'${filename1}$',lw=linew)
    print(filename1)
    a = (T1[0] - T1[-1]) / (grid1[-1]-grid1[0])
    print(f'Mean Strain Rate:\t{a}\n')
plt.tick_params(labelsize=fonts1)
plt.title('Velocity Field',fontsize=fonts1)
plt.xlabel(r'x (m)',fontsize=fonts1)
plt.ylabel(r'u (m/s)',fontsize=fonts1)
plt.legend(loc=0,fontsize=fonts1)
# plt.savefig('VelocityField.png', dpi=500)


plt.figure(figsize=figs)
for ic,filename1 in enumerate(filename):
    data1 = np.loadtxt(filename1,delimiter=',',skiprows=1)
    data1 = np.transpose(data1)
    grid1 = data1[0]
    plt.plot(grid1,data1[3],label=r'$T_{gas}$ @ '+fr'${filename1}$',c=cs[ic],lw=linew)
plt.tick_params(labelsize=fonts2)
plt.xlabel(r'x (m)',fontsize=fonts2)
plt.ylabel(r'T (K)',fontsize=fonts2)



# plt.figure(figsize=figs)
# for ic,filename1 in enumerate(filename):
#     data1 = np.loadtxt(filename1,delimiter=',',skiprows=1)
#     data1 = np.transpose(data1)
#     # f.grid
#     grid1 = data1[0]
#     # YAR
#     YAR = data1[6]
#     YARO = YAR[-1]

#     Z1 = (YAR - YARO) / (0.0-YARO)
#     plt.plot(grid1,Z1,label = fr'$Z$ @ ${filename1}$',c='k',ls=lss[ic],lw=linew)
# plt.tick_params(labelsize=fonts1)
# # plt.title('Mixture Fraction',fontsize=fonts1)
# plt.xticks(x_ticks,color='k')
# plt.xlabel(r'x (m)',fontsize=fonts1)
# plt.ylabel(r'Z (-)',fontsize=fonts1)
# plt.xlim(0,0.02)
# # plt.ylim(0,1)
# plt.legend(loc=0,fontsize=fonts1)
# plt.savefig('MixtureFraction.png', dpi=500)
plt.show()