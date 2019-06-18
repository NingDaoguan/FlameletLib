import numpy as np
import matplotlib.pyplot as plt

#fname = input('file name:\n>')
fname = 'trans.csv'
data = np.loadtxt(fname, delimiter=',', comments='#')
data = np.transpose(data)
x = data[0]
n = data[1]
rda = data[2]
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
ax1.plot(x,n)
ax2.plot(x,rda)
plt.show()

