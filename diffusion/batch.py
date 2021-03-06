import cantera as ct
import numpy as np
import os

# Create directory for output data files
data_directory = './diffusion_flame_batch_data/'
if not os.path.exists(data_directory):
    os.makedirs(data_directory)

# PART 1: INITIALIZATION
# reaction_mechanism = 'Ethanol_31.cti'
reaction_mechanism = 'Heptane0.cti'
gas = ct.Solution(reaction_mechanism)
gas.transport_model = 'UnityLewis'
width = 0.02
f = ct.CounterflowDiffusionFlame(gas, width=width)
# Define the operating pressure and boundary conditions
f.P = 5.0e6  # 1 bar
f.fuel_inlet.mdot = 24 # kg/m^2/s
f.fuel_inlet.X = 'C7H16:1'
f.fuel_inlet.T = 300.0  # K
f.oxidizer_inlet.mdot = 12 # kg/m^2/s
f.oxidizer_inlet.X = 'O2:0.21, N2:0.78, AR:0.01'
f.oxidizer_inlet.T = 800.0  # K
temperature_limit_extinction = np.maximum(f.fuel_inlet.T, f.oxidizer_inlet.T)  # K

# Set refinement parameters, if used
f.set_refine_criteria(ratio=4.0, slope=0.2, curve=0.2, prune=0.04)

# Define a limit for the maximum temperature below which the flame is
# considered as extinguished and the computation is aborted
# This increases the speed of refinement is enabled
def interrupt_extinction(t):
    if np.max(f.T) < temperature_limit_extinction:
        raise Exception('Flame extinguished')
    return 0.
f.set_interrupt(interrupt_extinction)

# Initialize and solve
print('Creating the initial solution')
f.solve(loglevel=0, auto=True)

# Save to data directory
file_name = 'initial_solution'
f.save(data_directory + file_name + '.xml', name='solution',
       description='Cantera version ' + ct.__version__ +
       ', reaction mechanism ' + reaction_mechanism)
f.write_csv(data_directory + file_name + '.csv', species='Y', quiet=False)


# PART 2: STRAIN RATE LOOP
# Compute counterflow diffusion flames at increasing strain rates at 1 bar
# The strain rate is assumed to increase by 25% in each step until the flame is
# extinguished
strain_factor = 1.25

# Exponents for the initial solution variation with changes in strain rate
# Taken from Fiala and Sattelmayer (2014)
exp_d_a = - 1. / 2.
exp_u_a = 1. / 2.
exp_V_a = 1.
exp_lam_a = 2.
exp_mdot_a = 1. / 2.

# Restore initial solution
file_name = 'initial_solution.xml'
f.restore(filename=data_directory + file_name, name='solution', loglevel=0)

# Counter to identify the loop
n = 0
# Do the strain rate loop
while np.max(f.T) > temperature_limit_extinction:
    n += 1
    print('strain rate iteration', n)
    # Create an initial guess based on the previous solution
    # Update grid
    f.flame.grid *= strain_factor ** exp_d_a
    normalized_grid = f.grid / (f.grid[-1] - f.grid[0])
    # Update mass fluxes
    f.fuel_inlet.mdot *= strain_factor ** exp_mdot_a
    f.oxidizer_inlet.mdot *= strain_factor ** exp_mdot_a
    # Update velocities
    f.set_profile('u', normalized_grid, f.u * strain_factor ** exp_u_a)
    f.set_profile('V', normalized_grid, f.V * strain_factor ** exp_V_a)
    # Update pressure curvature
    f.set_profile('lambda', normalized_grid, f.L * strain_factor ** exp_lam_a)
    try:
        # Try solving the flame
        f.solve(loglevel=0)
        file_name = 'strain_loop_' + format(n, '02d')
        f.save(data_directory + file_name + '.xml', name='solution', loglevel=1,
               description='Cantera version ' + ct.__version__ +
               ', reaction mechanism ' + reaction_mechanism)
        f.write_csv(data_directory + file_name + '.csv', species='Y', quiet=False)
    except Exception as e:
        if e.args[0] == 'Flame extinguished':
            print('Flame extinguished')
        else:
            print('Error occurred while solving:', e)
        break

