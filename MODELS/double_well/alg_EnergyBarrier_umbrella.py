from src.sim_double_well_SDE import sim_double_well, load_parameters
from src.metrics_double_well import activation_energy

import json, argparse, logging
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

    def bias_potential(x,**kw) : 
        return (kw["umbrella_strength"]/2)*(x-kw["umbrella_center"])**2
    def bias_force(x,t,**kw)   : 
        return -kw["umbrella_strength"]*(x-kw["umbrella_center"])
    S["bias_potential"] = bias_potential
    S["bias_force"]     = bias_force

import multiprocessing
def equilibrate(i): 
    REPLICAS[i].run(params["warmup_steps"], record=False)
    return i

def run_system(i):
    REPLICAS[i].run()
    REPLICAS[i].close()
    return i

P = multiprocessing.Pool()
sol = P.imap(equilibrate, range(len(REPLICAS)))
for n in sol: print "Equilibrated ", n

sol = P.imap(run_system, range(len(REPLICAS)))
for n in sol: print "Completed ", n
P.close(); P.join()


exit()
    

# Compute the exact value for error measurements
Um0, Ub, Um1 = map(S["potential"], [-1,0,1])
exact_activation_energy = np.array([Ub-Um0, Ub-Um1])

# Measure the barrier height from the trajectory
estimated_activation_energy, time_steps = activation_energy(**S)

# Compute the error
epsilon = np.abs(estimated_activation_energy-exact_activation_energy)

# Average over the two wells
epsilon = epsilon.sum(axis=1)

# Save the results
np.savetxt(S["f_results"],np.array([time_steps, epsilon]).T)

# Plot the results if asked
if "show_plot" in params and params["show_plot"]:
    from src.plots_double_well import plot_simulation
    plot_simulation(S)
    
    
    

