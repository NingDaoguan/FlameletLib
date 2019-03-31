import cantera as ct
import numpy as np
# import os

# Create directory for output data files
data_directory = './patch_data/'
# if not os.path.exists(data_directory):
    # os.makedirs(data_directory)

print('1:\tJet-A-300K\n2:\tJet-A-300KLe1\n3:\tJet-A-800K\n4:\tJet-A-800KLe1\n5:\tnc12h26-300K\n6:\tnc12h26-300KLe1')
x = input()
x = int(x)
width = 0.015
# PART 1: INITIALIZATION
if x==1:
    reaction_mechanism = 'KEROSENE_CRECK231.cti'
    gas = ct.Solution(reaction_mechanism)
    f = ct.CounterflowDiffusionFlame(gas, width=width)
    # Define the operating pressure and boundary conditions
    f.P = 1.e5  # 1 bar
    f.fuel_inlet.mdot = 5.5 # kg/m^2/s
    f.fuel_inlet.X = 'NC12H26:0.3, IC16H34:0.36, DECALIN:0.246, C7H8:0.094'
    f.fuel_inlet.T = 480.0  # K
    f.oxidizer_inlet.mdot = 3 # kg/m^2/s
    f.oxidizer_inlet.X = 'O2:0.21, N2:0.78, AR:0.01'
    f.oxidizer_inlet.T = 300.0  # K
    temperature_limit_extinction = 480  # K
elif x==2:
    reaction_mechanism = 'KEROSENE_CRECK231.cti'
    gas = ct.Solution(reaction_mechanism)
    gas.transport_model = 'UnityLewis'
    f = ct.CounterflowDiffusionFlame(gas, width=width)
    # Define the operating pressure and boundary conditions
    f.P = 1.e5  # 1 bar
    f.fuel_inlet.mdot = 0.6 # kg/m^2/s
    f.fuel_inlet.X = 'NC12H26:0.3, IC16H34:0.36, DECALIN:0.246, C7H8:0.094'
    f.fuel_inlet.T = 480.0  # K
    f.oxidizer_inlet.mdot = 0.35 # kg/m^2/s
    f.oxidizer_inlet.X = 'O2:0.21, N2:0.78, AR:0.01'
    f.oxidizer_inlet.T = 300.0  # K
    temperature_limit_extinction = 480  # K
elif x==3:
    reaction_mechanism = 'KEROSENE_CRECK231.cti'
    gas = ct.Solution(reaction_mechanism)
    f = ct.CounterflowDiffusionFlame(gas, width=width)
    # Define the operating pressure and boundary conditions
    f.P = 1.e5  # 1 bar
    f.fuel_inlet.mdot = 0.9 # kg/m^2/s
    f.fuel_inlet.X = 'NC12H26:0.3, IC16H34:0.36, DECALIN:0.246, C7H8:0.094'
    f.fuel_inlet.T = 480.0  # K
    f.oxidizer_inlet.mdot = 0.4 # kg/m^2/s
    f.oxidizer_inlet.X = 'O2:0.21, N2:0.78, AR:0.01'
    f.oxidizer_inlet.T = 800.0  # K
    temperature_limit_extinction = 800  # K
elif x==4:
    reaction_mechanism = 'KEROSENE_CRECK231.cti'
    gas = ct.Solution(reaction_mechanism)
    gas.transport_model = 'UnityLewis'
    f = ct.CounterflowDiffusionFlame(gas, width=width)
    # Define the operating pressure and boundary conditions
    f.P = 1.e5  # 1 bar
    f.fuel_inlet.mdot = 0.9 # kg/m^2/s
    f.fuel_inlet.X = 'NC12H26:0.3, IC16H34:0.36, DECALIN:0.246, C7H8:0.094'
    f.fuel_inlet.T = 480.0  # K
    f.oxidizer_inlet.mdot = 0.4 # kg/m^2/s
    f.oxidizer_inlet.X = 'O2:0.21, N2:0.78, AR:0.01'
    f.oxidizer_inlet.T = 800.0  # K
    temperature_limit_extinction = 800  # K
elif x==5:
    reaction_mechanism = 'nDodecane_CRECK.cti'
    gas = ct.Solution(reaction_mechanism)
    f = ct.CounterflowDiffusionFlame(gas, width=width)
    # Define the operating pressure and boundary conditions
    f.P = 1.e5  # 1 bar
    f.fuel_inlet.mdot = 0.6 # kg/m^2/s
    f.fuel_inlet.X = 'NC12H26:1.0'
    f.fuel_inlet.T = 480.0  # K
    f.oxidizer_inlet.mdot = 0.35 # kg/m^2/s
    f.oxidizer_inlet.X = 'O2:0.21, N2:0.78, AR:0.01'
    f.oxidizer_inlet.T = 300.0  # K
    temperature_limit_extinction = 480  # K
elif x==6:
    reaction_mechanism = 'nDodecane_CRECK.cti'
    gas = ct.Solution(reaction_mechanism)
    gas.transport_model = 'UnityLewis'
    f = ct.CounterflowDiffusionFlame(gas, width=width)
    # Define the operating pressure and boundary conditions
    f.P = 1.e5  # 1 bar
    f.fuel_inlet.mdot = 0.6 # kg/m^2/s
    f.fuel_inlet.X = 'NC12H26:1.0'
    f.fuel_inlet.T = 480.0  # K
    f.oxidizer_inlet.mdot = 0.35 # kg/m^2/s
    f.oxidizer_inlet.X = 'O2:0.21, N2:0.78, AR:0.01'
    f.oxidizer_inlet.T = 300.0  # K
    temperature_limit_extinction = 480  # K

# Set refinement parameters, if used
f.set_refine_criteria(ratio=4.0, slope=0.3, curve=0.3, prune=0.04)


# Initialize and solve
print('Restore an initial solution')
file_name = 'patchInit.xml'
f.restore(filename=data_directory + file_name, name='solution', loglevel=0)

f.solve(loglevel=0, auto=True)

f.save(data_directory + 'patch1' + '.xml', name='solution',
       description='Cantera version ' + ct.__version__ +
       ', reaction mechanism ' + reaction_mechanism)
f.write_csv(data_directory + 'patch1' + '.csv', species='Y', quiet=False)
