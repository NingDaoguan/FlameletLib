import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

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

print("Enter file names:")
filename = []
while True:
    name_input = input('> ')
    if name_input == '':
        break
    else:
        filename.append(name_input)

fig, ax = plt.subplots(figsize=(12,8))
for ifile in filename:
    data = np.loadtxt(ifile,delimiter=',',skiprows=1)
    data = np.transpose(data)
    ax.plot(data[0], data[1], label=ifile)
ax.set_title(r'Temperature evolution with increasing strain rate')
ax.legend()
ax.set_xlabel(r'$t$ $\mathrm{(s)}$')
ax.set_ylabel(r'$T$ $\mathrm{(K)}$')
plt.show()
