# Mon May 13 10:25:04 CST 2019
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.size'] = 18
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

# Files
figname = '1'
print("Enter file names:")
filename = []
while True:
    name_input = input('> ')
    if name_input == '':
        break
    else:
        filename.append(name_input)

# Plot data
file = filename[0]
with open(file) as f:
    line = f.readline()
    line = line.split(',')
    name = []
    for i in range(len(line)):
        name.append(line[i])

fig1 = plt.figure(figsize=(15,9))
ax1 = fig1.add_subplot(211)
ax2 = fig1.add_subplot(212)
for i,file in enumerate(filename):
    data = np.loadtxt(file,delimiter=',',comments='#',skiprows=1)
    data = np.transpose(data)
    x = data[0]
    for k in range(len(name)):
        if k in [2,7,9,10,18,19,20]:
            ax1.plot(x,data[k],label=name[k],ls='-')
    ax2.plot(x,data[1],label=name[1],ls='-',c='k')


# Fig
ax1.legend(loc=0)
ax2.legend(loc=0)
x_ticks = np.linspace(x[0], x[-1], 5)
ax1.set_xticks(x_ticks)
# ax1.set_ylim(-20, 60)
ax1.set_xlabel(r'$x$ (m)')
ax1.set_ylabel(r'$\dot{\omega}$ (kg/m3 s)')

ax2.set_xticks(x_ticks)
# ax2.set_ylim(-0.5e8, 6e8)
ax2.set_xlabel(r'$x$ (m)')
ax2.set_ylabel(r'$\dot{Q}$ (J/m3 s)')
fig1.tight_layout()
# plt.savefig(figname+'.png',dpi=500,bbox_inches='tight')
plt.show()
