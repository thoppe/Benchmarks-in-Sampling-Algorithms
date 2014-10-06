from src.sim_double_well_SDE import sim_double_well
from src.helper_functions import load_parameters, compute_activation_error
from src.helper_functions import save_results

import argparse, logging, random
import numpy as np

desc = '''
SAMPLING METHOD: Umbrella Sampling with Harmonic Bias Potentials
To be written...

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

# Create a simulation for every temperature, set the filenames
REPLICAS = []
for replica_n in range(params["umbrella_windows"]):
    p = params.copy()
    p["replica_n"] = replica_n
    p["f_trajectory"] = p["f_trajectory"].format(**p)
    p["f_results"] = p["f_results"].format(**p)
    REPLICAS.append( sim_double_well(**p) )

# Modify the potential of each simulation
ubounds  = REPLICAS[0]["umbrella_bounds"]
uwindows = REPLICAS[0]["umbrella_windows"]
U_X = np.linspace(ubounds[0],ubounds[1],uwindows)
for n,S in enumerate(REPLICAS):
        
    S["umbrella_strength"] = S["umbrella_strength"]
    S["umbrella_center"]   = U_X[n]

    def bias_potential(x,t,**kw) : 
        return (kw["umbrella_strength"]/2)*(x-kw["umbrella_center"])**2
    def bias_force(x,t,**kw)   : 
        return -kw["umbrella_strength"]*(x-kw["umbrella_center"])
    S["bias_potential"] = bias_potential
    S["bias_force"]     = bias_force


# Let the systems equlibrate on their own
for S in REPLICAS: 
    S.run(params["warmup_time"], record=False)

# Run each simulation

for S in REPLICAS:
    S.run()
    S.close()


for S in REPLICAS:

    # Compute the error, summed over both wells
    epsilon,time_steps = compute_activation_error(S)

    # Save the results to file
    save_results(S["f_results"], time_steps, epsilon)
  
# Plot the results if asked
if "show_plot" in params and params["show_plot"]:
    from src.plots_double_well import plot_simulation
    plot_simulation(S)
    
    
    

