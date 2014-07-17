from src.sim_double_well_SDE import sim_double_well, load_parameters
from src.metrics_double_well import average_activation_energy

import json, argparse, logging, random
import numpy as np

desc = '''
SAMPLING METHOD: Replica Exchange

[kT_list]         defines the temperatures used.
[exchange_steps]  number of simulation steps before a replica exchanged is attempted

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

# Create a simulation for every temperature

REPLICAS = []
for kT in params["kT_list"]:
    p = params.copy()
    p["kT"] = kT
    REPLICAS.append( sim_double_well(**p) )


def exchange_replicas(s0,s1):
    s0["kT"], s1["kT"] = s1["kT"], s0["kT"]

exchange_steps = params["exchange_steps"]

while not REPLICAS[0].is_complete():
    for S in REPLICAS: S.run(exchange_steps)

    # Choose two replicas and compute probability to exchange
    s0, s1 = random.sample(REPLICAS, 2)
    u0 = s0["potential"](s0["xi"])
    u1 = s1["potential"](s1["xi"])
    delta = (1/s0["kT"])-(1/s1["kT"]) * (u1-u0)
    p = np.exp(-delta)

    # Exchange replica coordinates if Metropolis-Hastings condition is met
    if np.random.random() < p:
        exchange_replicas(s0,s1)

# Compute the exact value for error measurements
Um0, Ub, Um1 = map(REPLICAS[0]["potential"], [-1,0,1])
exact_avg_activation_energy = np.array([Ub-Um0, Ub-Um1]).mean()

# Compute the error
for k,S in enumerate(REPLICAS): 
    err_T = S.traj_metric_t
    err   = np.abs((np.array(S.traj_metric) - exact_avg_activation_energy))

    # Format the output filename
    f_out = S["f_results"].format(replica_n = k)
    np.savetxt(f_out,np.array([err_T, err]).T)

# Plot the results if asked
if "show_plot" in params and params["show_plot"]:
    from src.plots_double_well import plot_simulation

    err_T = S.traj_metric_t
    err   = np.abs((np.array(S.traj_metric) - exact_avg_activation_energy))

    for S in REPLICAS: 
        plot_simulation(S, err)
    
    
    

