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

cargs, params = startup_simulation(desc, format_filenames=False)

# Set the parameters for every simulation and the filenames
param_set = []
for replica_n in range(params["umbrella_windows"]):
    p = params.copy()
    p["replica_n"] = replica_n
    p["f_trajectory"] = p["f_trajectory"].format(**p)
    p["f_results"] = p["f_results"].format(**p)
    param_set.append(p)

# Create a simulation for every bias
REPLICAS = []
for p in param_set:
    REPLICAS.append( sim_double_well(**p) )

# Modify the potential of each simulation
U_X = np.linspace(*params["umbrella_bounds"], 
                  num=params["umbrella_windows"])

for n,S in enumerate(REPLICAS):
        
    S.umbrella_strength = params["umbrella_strength"]
    S.umbrella_center   = U_X[n]

    class umbrella_potential:
        def __init__(self,x0,A):
            self.strength = A
            self.center   = x0
        def __call__(self,x,t):
            return (self.strength/2)*(x-self.center)**2

    class umbrella_force:
        def __init__(self,x0,A):
            self.strength = A
            self.center   = x0
        def __call__(self,x,t):
            return -self.strength*(x-self.center)

    S.bias_potential = umbrella_potential(U_X[n], params["umbrella_strength"])
    S.bias_force = umbrella_force(U_X[n], params["umbrella_strength"])

# Let the systems equlibrate on their own
for S in REPLICAS:
    S.run(params["warmup_time"], record=False)

# Run each simulation
for S in REPLICAS:
    S.run()



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
    
    
    

