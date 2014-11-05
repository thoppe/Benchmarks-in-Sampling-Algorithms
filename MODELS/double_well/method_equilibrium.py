from src.sim_double_well_SDE import sim_double_well
from src.helper_functions import startup_simulation, finalize_simulation
from src.metrics_double_well import compute_activation_error

import argparse, logging
import numpy as np

desc = '''
SAMPLING METHOD: Equilibrium
This simulation is the test case, e.g. 
standard dynamics with no enchanced sampling.

Computes the energy barrier between a double well using overdamped 
Langevin dynamics with an Euler-Maruyama SDE.
'''

params = startup_simulation(desc)
S = sim_double_well(**params)

# First equilibrate the system
S.run(fixed_time=params["warmup_time"], record=False)

# Run the simulation
S.run()

# Finish the simulation, make plots, etc...
finalize_simulation(S, metric_function=compute_activation_error)    
    

