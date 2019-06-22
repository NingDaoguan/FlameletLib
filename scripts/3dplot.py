import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import FormatStrFormatter


plt.rcParams['font.family'] = 'serif'
from matplotlib import rc
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.labelpad'] = 20
fs = 20

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

num = input('number:\n>')
for i in range(int(num)):
    fname = 'flameletTable_{:}.csv'.format(i)
    data = np.loadtxt(fname,delimiter=',',skiprows=1)
    Z = data.T[0]
    ha = data.T[1] / 1000.0
    Yc = data.T[2]
    omega = data.T[3]
    T = data.T[4]
    norm = matplotlib.colors.Normalize(vmin=300,vmax=2200)

    sc = ax.scatter(Z,ha,Yc,c=T,cmap=plt.cm.rainbow,norm=norm)

v = [300,500,1000,1500,2000]
cbar = plt.colorbar(sc,ticks=v)
cbar.set_label(r'T (K)',fontsize=fs)
for t in cbar.ax.get_yticklabels():
    t.set_fontsize(fs)
ax.set_xlabel(r'$Z$ (-)',fontsize=fs)
ax.set_xlim(0,1)
ax.set_ylabel(r'$h_a (kJ/kg)$',fontsize=fs)
ax.set_ylim(-5000,500)
ax.set_zlabel(r'$Y_c$ (-)',fontsize=fs)
ax.tick_params(labelsize=fs-6)
#ax.yaxis.set_major_formatter(FormatStrFormatter('%f'))
plt.show()
