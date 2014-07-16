from src.sim_double_well_SDE import sim_double_well, load_parameters
from src.metrics_double_well import average_activation_energy

import json
import numpy as np


'''
SAMPLING METHOD: Standard
This simulation is the test case, e.g. 
the standard dynamics with no enchanced sampling.

Computes the energy barrier between a double well using overdamped 
Langevin dynamics with an Euler-Maruyama SDE.
'''

# Load the simulation parameters

f_json = "simulation_setups/example_EnergyBarrier.json"
params = load_parameters(f_json)

S = sim_double_well(**params)
S.run()

# Compute the exact value for error measurements
Um0, Ub, Um1 = map(S["potential"], [-1,0,1])
exact_avg_activation_energy = np.array([Ub-Um0, Ub-Um1]).mean()

# Compute the error
err_T = S.traj_metric_t
err   = np.abs((np.array(S.traj_metric) - exact_avg_activation_energy))

# Save the results
np.savetxt(params["f_results"],[err_T, err])

'''
Plot the results
'''

from src.plots_double_well import plot_simulation
plot_simulation(S, err)
    
    
    

