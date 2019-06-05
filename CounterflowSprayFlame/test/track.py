import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('LAGml-1.1768mr-1.1768L-0.02_raw.csv',delimiter=',',comments='#')
data = np.transpose(data)
t = data[0]
p = data[1]
ax1 = plt.subplot(121)
ax2 = plt.subplot(122)
ax1.plot(t,data[1])
ax2.plot(t,data[3])
plt.show()
