from src.sim_double_well_SDE import sim_double_well
from src.helper_functions import load_parameters, compute_activation_error
from src.helper_functions import save_results

import argparse, logging, random
import numpy as np

desc = '''
SAMPLING METHOD: Replica Exchange

[kT_list]       defines the temperatures used.
[exchange_time] simulation time before a replica exchanged is attempted

Computes the energy barrier between a double well using overdamped 
Langevin dynamics with an Euler-Maruyama SDE.
'''

parser = argparse.ArgumentParser(description=desc)
parser.add_argument('parameter_file_json')
cargs = vars(parser.parse_args())

# Start the logger
logging.root.setLevel(logging.INFO)

# Load the simulation parameters
params = load_parameters(cargs["parameter_file_json"])

# Create a simulation for every temperature, set the filenames
REPLICAS = []
for replica_n,kT in enumerate(params["kT_list"]):
    p = params.copy()
    p["kT"] = kT
    p["replica_n"] = replica_n
    p["f_trajectory"] = p["f_trajectory"].format(**p)
    p["f_results"] = p["f_results"].format(**p)
    REPLICAS.append( sim_double_well(**p) )

# Let the systems equlibrate on their own
for S in REPLICAS: 
    S.run(params["warmup_time"], record=False)

while not REPLICAS[0].is_complete():
    for S in REPLICAS: 
        S.run(fixed_time = params["exchange_time"])

    # Choose two replicas and compute probability to exchange
    s0, s1 = random.sample(REPLICAS, 2)
    u0 = s0["potential"](s0["xi"])
    u1 = s1["potential"](s1["xi"])
    beta0 = 1.0/s0["kT"]
    beta1 = 1.0/s1["kT"]
    p = np.exp((beta0-beta1)*(u0-u1))

    def exchange_replicas(s0,s1):
        s0["xi"], s1["xi"] = s1["xi"], s0["xi"]

    # Exchange replica coordinates if Metropolis-Hastings condition is met
    if np.random.random() < p:
        exchange_replicas(s0,s1)

for S in REPLICAS: 
    
    # Close any files open
    S.close()

    # Compute the error, summed over both wells
    epsilon,time_steps = compute_activation_error(S)

    # Save the results to file
    save_results(S["f_results"], time_steps, epsilon)

# Plot the results if asked
if "show_plot" in params and params["show_plot"]:
    from src.plots_double_well import plot_simulation
    for S in REPLICAS: 
        plot_simulation(S)
    
    
    

