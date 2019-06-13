import numpy as np
import matplotlib.pyplot as plt

fname = input('file name:\n>')
data = np.loadtxt(fname, delimiter=',', comments='#')
data = np.transpose(data)
x = data[0]
h = data[1]
m = data[2]
# total = 0.0
# for i in range(len(x)):
#     if i == 0:
#         total += 0.5*(x[1] - x[0]) * m[0]
#     elif i == len(m)-1:
#         total += 0.5*(x[-1] - x[-2]) * m[-1]
#     else:
#         total += 0.5*(x[i+1] - x[i-1]) * m[i]
ax1 = plt.subplot(121)
ax2 = plt.subplot(122)
ax1.plot(x,h)
ax2.plot(x,m)
plt.show()

