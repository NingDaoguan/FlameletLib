#Sat Dec 15 15:50:09 CST 2018
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

x = input("Enter loop number:\n")
x = int(x)
y = input("T(0) or omega(1)?\n")
filename = []
if x > 1:
    for n in range(x+1):
       filename.append('flameletTable_{:}.csv'.format(n))
else:
    filename = ['lookupResult.csv']

fonts1 = 23
fonts2 = 25
linew = 3
figs = (13,10)

plt.figure(figsize=figs)
plt.rcParams['font.family'] = "Serif"
from matplotlib import rc
plt.rcParams['mathtext.fontset'] = 'stix'
for ic,filename1 in enumerate(filename):
    data1 = np.loadtxt(filename1,delimiter=',',skiprows=1)
    data1 = np.transpose(data1)
    YAR = data1[3]
    YARO = YAR[0]
    Z1 = data1[0]
    Yc = data1[2]

    if y == '0':
        v = [300,500,1000,1500,2000]
        norm = matplotlib.colors.Normalize(vmin=300, vmax=2300)
        sc = plt.scatter(Z1,Yc,c=data1[4],cmap=plt.cm.rainbow,s=2,norm=norm)
    elif y == '1':
        #v = [0.0,0.5e9,1e9,1.5e9,2e9,2.5e9]
        #norm = matplotlib.colors.Normalize(vmin=0, vmax=2.2e9)
        v = [0,100,200,300,350]
        norm = matplotlib.colors.Normalize(vmin=0, vmax=350)
        sc = plt.scatter(Z1,Yc,c=data1[3],cmap=plt.cm.rainbow,s=2,norm=norm)

cbar = plt.colorbar(sc,ticks=v)
if y == '0':
    cbar.set_label(r'T (K)', fontsize=fonts2)
elif y == '1':
    cbar.set_label(r'$\dot \omega_{Y_c} \ (kg/m^3 s)$', fontsize=fonts2)
else:
    print("Input Error")
for t in cbar.ax.get_yticklabels():
    t.set_fontsize(fonts2)
# plt.xlim(0,1)
# plt.ylim(0,0.5)
plt.tick_params(labelsize=fonts2)
plt.xlabel(r'Z (-)',fontsize=fonts2)
plt.ylabel(r'$Y_c$ (-)',fontsize=fonts2)
# plt.savefig('SLF_Z'+y+'.png',dpi=500,bbox_inches='tight')
plt.show()
