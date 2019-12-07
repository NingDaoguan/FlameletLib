#Sat Dec 15 15:50:09 CST 2018
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Configurations
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.size'] = 20
mpl.rcParams['font.weight'] = 'medium'
mpl.rcParams['font.style'] = 'normal'
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['mathtext.fallback_to_cm'] = True
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['savefig.bbox'] = 'tight'
figs = (13,9)

# Input file names
N = 201
print("Enter file name:")
filename = input('> ')
data = np.loadtxt(filename, delimiter=',').T
z = data[0]
c = data[1]
t = data[3]

points1 = np.linspace(np.min(z), np.max(z), N)
points2 = np.linspace(np.min(c), np.max(c), N)
xx, yy = np.meshgrid(points1, points2)
zz = xx + yy

for i in range(N):
    for j in range(N):
        zz[j][i] = t[N*i + j]

plt.figure(figsize=figs)
ct = plt.contourf(xx, yy, zz, vmin = 300, vmax = 2200, cmap = 'rainbow')
# v = [300,400,500,1000,1500,2000,2300]
cbar = plt.colorbar(ct)
cbar.set_label(r'$T$ $\mathrm{(K)}$')
plt.xlim(0,0.2)
plt.ylim(0,0.3)
plt.xlabel(r'$Z$ (-)')
plt.ylabel(r'$Y_c$ (-)')
plt.tight_layout()
plt.show()

