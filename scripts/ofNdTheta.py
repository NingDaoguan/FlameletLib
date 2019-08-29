import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
# Configurations
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.size'] = 22
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

# data loading 
dataNd = np.loadtxt('pNd.csv',delimiter=',').T
dataTheta = np.loadtxt('pTheta.csv',delimiter=',').T


# Spatial variation
x = np.linspace(0.0, 0.02, 41)
nds = np.zeros(len(x))
thetas = np.zeros(len(x))
for i in range(len(x)):
    nds[i] = np.mean(dataNd[i+1])
    thetas[i] = np.mean(dataTheta[i+1])

fig, axes = plt.subplots(1,2,figsize=figs)
axes[0].scatter(x, nds)
axes[1].scatter(x, thetas)
axes[0].set_ylabel('Number density $\mathrm{(m^{-3})}$')
axes[0].set_xlabel('x (m)')
axes[0].set_xlim(0,0.02)
axes[1].set_ylabel('Volume fraction (-)')
axes[1].set_xlabel('x (m)')
axes[1].set_xlim(0,0.02)
axes[1].set_ylim(0,0.0002)
fig.tight_layout()
plt.show() 

# Time evolution
fig, axes = plt.subplots(1,2,figsize=figs)
for i in range(len(dataNd[0])):
    axes[0].scatter(x, dataNd.T[i][1::], s=5, c='k')
    axes[1].scatter(x, dataTheta.T[i][1::], s=5, c='k')
axes[0].set_ylabel('Number density $\mathrm{(m^{-3})}$')
axes[0].set_xlabel('x (m)')
axes[0].set_xlim(0,0.02)
axes[0].set_ylim(0,3e10)
axes[1].set_ylabel('Volume fraction (-)')
axes[1].set_xlabel('x (m)')
axes[1].set_xlim(0,0.02)
axes[1].set_ylim(0, 0.001)
fig.tight_layout()
plt.show()
