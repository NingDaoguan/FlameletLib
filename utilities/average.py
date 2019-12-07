import numpy as np
import matplotlib.pyplot as plt

time = [0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3]

vec = []

for it in range(len(time)):
    t = time[it]
    f = f'{t}/line_rho_T_N2_AR_HE_H_O2_OH_O_H2_H2O_HO2_CO_CO2_HCO_CH3_CH4_CH2O_T-CH2_S-CH2_C2H4_CH3O_C2H5_C2H6_H2O2_C2H2_C2H3_CH2CHO_CH2CO_CH2OH_CH3CHO_CH3CO_C2H5OH_CH2CH2OH_CH3CHOH_CH3CH2O.xy'
    data = np.loadtxt(f).T
    vec.append(data)

output = vec[0].copy()
for ix in range(len(output[0])):
    for iv in range(len(output)):
        avg = 0.0
        for it in range(len(time)):
            avg += vec[it][iv][ix]
        avg /= len(time)
        output[iv][ix] = avg

with open('singleGraph.csv', 'w+') as f:
    np.savetxt(f, output.T, delimiter=',', fmt='%f')

plt.plot(output[0], output[2])
plt.show()

