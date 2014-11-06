from src.sim_double_well_SDE import sim_double_well
from src.helper_functions import startup_simulation, finalize_simulation
from src.metrics_double_well import compute_activation_error

import logging, random
import numpy as np

desc = '''
SAMPLING METHOD: Replica Exchange

[kT_list]       defines the temperatures used.
[exchange_time] simulation time before a replica exchanged is attempted

Computes the energy barrier between a double well using overdamped 
Langevin dynamics with an Euler-Maruyama SDE.
'''

cargs, params = startup_simulation(desc, format_filenames=False)

# Load multiprocessing library if requested
if cargs["parallel"]: 
    #import multiprocessing as mp
    logging.warning("Multiprocessing not implemented yet.")

# Set the parameters for every simulation and the filenames
param_set = []
for replica_n,kT in enumerate(params["kT_list"]):
    p = params.copy()
    p["kT"] = kT
    p["replica_n"] = replica_n
    p["f_trajectory"] = p["f_trajectory"].format(**p)
    p["f_results"] = p["f_results"].format(**p)
    param_set.append(p)

# Create a simulation for every temperature
REPLICAS = []
for p in param_set:
    REPLICAS.append( sim_double_well(**p) )

# Let the systems equilibrate on their own
for S in REPLICAS: 
    S.run(params["warmup_time"], record=False)

while not REPLICAS[0].is_complete():
    for S in REPLICAS: 
        S.run(fixed_time = params["exchange_time"])

    # Choose two replicas and compute probability to exchange
    s0, s1 = random.sample(REPLICAS, 2)
    u0 = s0.potential(s0.xi)
    u1 = s1.potential(s1.xi)
    beta0 = 1.0/s0.kT
    beta1 = 1.0/s1.kT
    p = np.exp((beta0-beta1)*(u0-u1))

    # Exchange replica coordinates if Metropolis-Hastings condition is met
    if np.random.random() < p:
        s0.xi, s1.xi = s1.xi, s0.xi

for S in REPLICAS:
    S.close()

# Compute and save the errors
for S,p in zip(REPLICAS,param_set):    
    print S.f_trajectory
    compute_activation_error(S, **p)

# Plot the results if requested
if params["show_plot"]:
    import src.plots_double_well as plot
    for S,p in zip(REPLICAS,param_set):
        plot.plot_simulation(S,**p)



