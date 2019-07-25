import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import FormatStrFormatter
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
plt.rcParams['axes.labelpad'] = 18
mks= ["o", "D", "d", "s", "p", "H", 0, 4, "<", "3",
      1, 5, ">", "4", 2, 6, "^", "2", 3, 7, "v", "1", "None", None, " ", ""]

fig = plt.figure(figsize=(13,8))
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
    sc = ax.scatter(Z,ha,Yc,c=T,cmap='rainbow',vmin=300,vmax=2600,s=60)

v = [300,500,1000,1500,2000,2500]
cbar = plt.colorbar(sc,ticks=v)
cbar.set_label(r'T (K)')
# for t in cbar.ax.get_yticklabels():
#     t.set_fontsize(fs)
ax.set_xlabel(r'$Z$ (-)')
ax.set_ylabel(r'$h (kJ/kg)$')
ax.set_zlabel(r'$Y_c$ (-)')
ax.tick_params(labelsize=14)
#ax.yaxis.set_major_formatter(FormatStrFormatter('%f'))
plt.savefig('3d.png',dpi=500)
plt.show()