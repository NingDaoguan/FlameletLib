import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('transLAGml-1.1768mr-1.1768L-0.02_raw.csv',delimiter=',',comments='#')
data = np.transpose(data)
x = data[0]
h = data[1]
m = data[2]
total = 0.0
for i in range(len(x)):
    if i == 0:
        total += 0.5*(x[1] - x[0]) * m[0]
    elif i == len(m)-1:
        total += 0.5*(x[-1] - x[-2]) * m[-1]
    else:
        total += 0.5*(x[i+1] - x[i-1]) * m[i]
#plt.plot(x,m)
plt.plot(x,m)
plt.show()
print(total)
