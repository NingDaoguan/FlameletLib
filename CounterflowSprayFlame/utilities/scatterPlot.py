#Sat Dec 15 15:50:09 CST 2018
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

x = input("Enter loop number:\n")
x = int(x)
y = input("Yc(0) or chi(1)?\n")
#filename = []
#for n in range(x+1):
#    filename.append('flameletTable_{:}.csv'.format(n))
filename = ['lookupResult.csv']

fonts1 = 23
fonts2 = 25
linew = 3
figs = (13,10)

plt.figure(figsize=figs)
plt.rcParams['font.family'] = "Serif"
for ic,filename1 in enumerate(filename):
    data1 = np.loadtxt(filename1,delimiter=',',skiprows=1)
    data1 = np.transpose(data1)
    YAR = data1[3]
    YARO = YAR[0]
    Z1 = data1[0]
    Yc = data1[1]
    norm = matplotlib.colors.Normalize(vmin=300, vmax=2500)
    sc = plt.scatter(Z1,Yc,c=data1[2],cmap=plt.cm.rainbow,s=3,norm=norm)

v = [300,500,1000,1500,2000,2500]
cbar = plt.colorbar(sc,ticks=v)
cbar.set_label(r'T (K)', fontsize=fonts2)
for t in cbar.ax.get_yticklabels():
    t.set_fontsize(fonts2)
plt.tick_params(labelsize=fonts2)
plt.xlabel(r'$Z (-)$',fontsize=fonts2)
if y == '0':
    plt.ylabel(r'$Y_c (-)$',fontsize=fonts2)
elif y == '1':
    plt.ylabel(r'$\chi (1/s)$',fontsize=fonts2)
else:
    print("Input Error")

# plt.xlim(0,1)
# plt.ylim(0,0.35)
plt.savefig('SLF_Z'+y+'.png',dpi=500)
#plt.show()
