import pylab as plt
import seaborn as sns
import numpy as np
import glob,os

from src.sim_double_well_SDE import load_parameters
import multiprocessing

cutoff = 100000

# Convert the trajectories
def convert_txt_to_numpy_traj(f):
    f_np = os.path.basename(f).replace('.txt','.npy')  
    f_np = os.path.join("trajectory", f_np)
    if not os.path.exists(f_np):
        print f_np
        time, pos = np.loadtxt(f).T
        np.save(f_np,pos)
    return f

P = multiprocessing.Pool()
F_TRAJ = glob.glob("trajectory/umbrella_EnergyBarrier_r*")
sol = P.imap(convert_txt_to_numpy_traj, F_TRAJ)
for x in sol: pass

from pymbar import MBAR

params = load_parameters("simulation_setups/umbrella_EnergyBarrier.json")
F_TRAJ = sorted(glob.glob("umbrella_EnergyBarrier_r*.npy"))

#### Build the replicas here for templating
from src.sim_double_well_SDE import sim_double_well, load_parameters
from src.metrics_double_well import activation_energy

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


F_TRAJ = sorted(glob.glob("trajectory/umbrella_EnergyBarrier_r*.npy"))
Y = np.array([np.load(f)[:cutoff] for f in F_TRAJ])
print "Cutoff/datapoints ", cutoff, np.load(f).shape

k_states = len(Y)
l_states = len(Y)
N_k = np.array([Y[k].shape[0] for k in range(k_states)])
N_max = max(N_k)
u_kln = np.zeros((k_states,l_states,N_max))

print "Building u_kln"

for k in xrange(k_states):
    for l in xrange(l_states):
        u_kln[k,l,:] =  REPLICAS[l].U(Y[k])


# Estimate free energies from simulation using MBAR.
print "Estimating relative free energies using MBAR, this may take a while..."

import pymbar
M = pymbar.MBAR(u_kln, N_k)
dF, dF_un = M.getFreeEnergyDifferences()
print pymbar.timeseries.statisticalInefficiency(Y[0],Y[1])

x_idx = np.argmin(np.abs(U_X))
F = dF[x_idx]
F -= F[x_idx]
F += 1

plt.errorbar(U_X, F,yerr=dF_un[x_idx],label="MBAR {} samples".format(cutoff))

'''
def weights(X, mu, var):
    sigma = np.sqrt(var)
    w = np.exp(-0.5*(((X-mu)/sigma)**2))
    return w/(sigma*np.sqrt(2*np.pi))

def diff_free_energy(X, mu, var, idx):
    beta = 1.0/params["kT"]
    x0 = centers[idx]
    k  = params["umbrella_strength"]
    A = (1/beta)*((X-mu)/var) + k*(X-x0)
    return A

observed_mean = [y.mean() for y in Y]
observed_var  = [y.var()  for y in Y]

P = [weights(X,mu,var) for (mu,var) in zip(observed_mean, observed_var)]
P = np.array(P)
P /= P.sum(axis=0)

dA = [diff_free_energy(X,mu,var,idx) for idx,(mu,var) in 
      enumerate(zip(observed_mean, observed_var))]

dA_avg = (dA*P).sum(axis=0)

from scipy.integrate import cumtrapz, simps
A = cumtrapz(dA_avg,X)

# Center for visual appeal
#A -= A.min()
#A -= A[5000] -1
'''
def U(x): return (x**2-1.0)**2 

#plt.plot(X,dA_avg,label="Umbrella integration")
plt.plot(U_X,U(U_X),'--',label="Exact")
#plt.plot(X[1:],A,label="Umbrella integration")
#plt.xlim(-2,2)
#plt.ylim(0,2)
plt.legend(loc="best")
plt.show()
