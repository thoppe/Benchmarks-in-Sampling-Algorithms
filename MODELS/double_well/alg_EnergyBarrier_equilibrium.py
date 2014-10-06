from src.sim_double_well_SDE import sim_double_well
from src.helper_functions import load_parameters, compute_activation_error
from src.helper_functions import save_results

import argparse, logging
import numpy as np

desc = '''
SAMPLING METHOD: Equilibrium
This simulation is the test case, e.g. 
standard dynamics with no enchanced sampling.

Computes the energy barrier between a double well using overdamped 
Langevin dynamics with an Euler-Maruyama SDE.
'''

parser = argparse.ArgumentParser(description=desc)
parser.add_argument('parameter_file_json')
parser.add_argument('--replica_n', type=int, default=0,
                    help="Appends this number onto the results file")
cargs = vars(parser.parse_args())

# Start the logger
logging.root.setLevel(logging.INFO)

# Load the simulation parameters
params = load_parameters(cargs["parameter_file_json"])

# Set the filenames
params["f_trajectory"] = params["f_trajectory"].format(**cargs)
params["f_results"] = params["f_results"].format(**cargs)

S = sim_double_well(**params)

# First equlibrate the system
S.run(fixed_time=params["warmup_time"], record=False)

# Run the simulation
S.run()
S.close()

# Compute the error, summed over both wells
epsilon,time_steps = compute_activation_error(S)

# Save the results to file
save_results(params["f_results"], time_steps, epsilon)

# Plot the results if asked
if "show_plot" in params and params["show_plot"]:
    from src.plots_double_well import plot_simulation
    plot_simulation(S)

    
    

