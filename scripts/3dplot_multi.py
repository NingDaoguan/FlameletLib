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

dim = input('dim:\n>')
for idm in range(int(dim)):
    for i in range(100):
        fname = 'ha_{}/flameletTable_{:}.csv'.format(idm,i)
        try:
            data = np.loadtxt(fname,delimiter=',',skiprows=1)
            Z = data.T[0]
            ha = data.T[1] / 1000.0
            Yc = data.T[2]
            omega = data.T[3]
            T = data.T[4]
            norm = matplotlib.colors.Normalize(vmin=250,vmax=2200)

            sc = ax.scatter(Z,ha,Yc,c=T,cmap=plt.cm.rainbow,norm=norm)
        except IOError:
            break

v = [250,500,1000,1500,2000]
cbar = plt.colorbar(sc,ticks=v)
cbar.set_label(r'T (K)',fontsize=fs)
for t in cbar.ax.get_yticklabels():
    t.set_fontsize(fs)
ax.set_xlabel(r'$Z$ (-)',fontsize=fs)
ax.set_ylabel(r'$h_a (kJ/kg)$',fontsize=fs)
ax.set_zlabel(r'$Y_c$ (-)',fontsize=fs)
ax.tick_params(labelsize=fs-6)
#ax.yaxis.set_major_formatter(FormatStrFormatter('%f'))
plt.show()
