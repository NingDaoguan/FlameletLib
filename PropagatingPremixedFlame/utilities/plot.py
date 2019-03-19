#Sat Dec 15 15:50:09 CST 2018
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

x = input("Enter loop number:\n")
x = int(x)
y = input("Yc(0) or chi(1)?\n")
filename = []
for n in range(x):
    filename.append('flameletTable_{:}.csv'.format(n))

fonts1 = 23
fonts2 = 25
linew = 3
figs = (13,10)

plt.figure(figsize=figs)
plt.rcParams['font.family'] = "Serif"
for ic,filename1 in enumerate(filename):
    data1 = np.loadtxt(filename1,delimiter=',',skiprows=1)
    data1 = np.transpose(data1)
    grid1 = data1[0]
    YAR = data1[3]
    YARO = YAR[0]
    Z1 = (YAR - YARO) / (0.0-YARO)
    Yc = data1[1]
    norm = matplotlib.colors.Normalize(vmin=300, vmax=2300)
    sc = plt.scatter(Z1,Yc,c=data1[2],cmap=plt.cm.rainbow,s=10,norm=norm)

v = [300,500,1000,1500,2000]
cbar = plt.colorbar(sc,ticks=v)
cbar.set_label('T (K)', fontsize=fonts2)
for t in cbar.ax.get_yticklabels():
    t.set_fontsize(fonts2)
plt.tick_params(labelsize=fonts2)
plt.xlabel(r'Z (-)',fontsize=fonts2)
if y == '0':
    plt.ylabel(r'Yc (-)',fontsize=fonts2)
elif y == '1':
    plt.ylabel(r'chi (-)',fontsize=fonts2)
else:
    print("Input Error")

plt.show()
