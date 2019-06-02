import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('LAGml-1.1768mr-1.1768L-0.02_raw.csv',delimiter=',',comments='#')
data = np.transpose(data)
t = data[0]
p = data[1]
plt.plot(t,p)
plt.show()
plt.plot(t,data[4])
plt.show()
