from src.sim_double_well_SDE import sim_double_well, load_parameters
from src.metrics_double_well import activation_energy

import json, argparse, logging
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

with sim_double_well(**params) as S:

    # First equlibrate the system
    S.run(fixed_time=params["warmup_time"], record=False)

    # Run the simulation
    S.run()

    # Compute the exact value for error measurements
    Um0, Ub, Um1 = map(S["potential"], [-1,0,1])
    exact_activation_energy = np.array([Ub-Um0, Ub-Um1])

    # Measure the barrier height from the trajectory
    estimated_activation_energy, time_steps = activation_energy(**S)

    # Compute the error, summed over both wells
    epsilon = np.abs(estimated_activation_energy-exact_activation_energy)
    epsilon = epsilon.sum(axis=1)

    # Save the results
    np.savetxt(S["f_results"],np.array([time_steps, epsilon]).T)


# Plot the results if asked
if "show_plot" in params and params["show_plot"]:
    from src.plots_double_well import plot_simulation
    plot_simulation(S)

    
    

