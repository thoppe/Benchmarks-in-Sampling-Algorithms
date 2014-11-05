from src.sim_double_well_SDE import sim_double_well
from src.helper_functions import startup_simulation, finalize_simulation
from src.metrics_double_well import compute_activation_error

import logging, random
import numpy as np

desc = '''
SAMPLING METHOD: Umbrella Sampling with Harmonic Bias Potentials
To be written...

Computes the energy barrier between a double well using overdamped 
Langevin dynamics with an Euler-Maruyama SDE.
'''

params = startup_simulation(desc, format_filenames=False)

# Create a simulation for every temperature, set the filenames
REPLICAS = []
for replica_n in range(params["umbrella_windows"]):
    p = params.copy()
    p["replica_n"] = replica_n
    p["f_trajectory"] = p["f_trajectory"].format(**p)
    p["f_results"] = p["f_results"].format(**p)
    REPLICAS.append( sim_double_well(**p) )

# Modify the potential of each simulation
U_X = np.linspace(*params["umbrella_bounds"], 
                  num=params["umbrella_windows"])

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

# Finish the simulation, make plots, etc...
for S in REPLICAS:    
    finalize_simulation(S, metric_function=compute_activation_error)
    
    
    

