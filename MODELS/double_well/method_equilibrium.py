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

cargs, params = startup_simulation(desc)
S = sim_double_well(**params)

# First equilibrate the system
S.run(fixed_time=params["warmup_time"], record=False)

# Run the simulation
S.run()

# Close any open file handlers
S.close()

# Compute and save the errors
compute_activation_error(S, **params)

# Plot the results if requested
if params["show_plot"]:
    import src.plots_double_well as plot
    plot.plot_simulation(S,**params)

    

