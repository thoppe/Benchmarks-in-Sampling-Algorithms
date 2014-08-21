import pylab as plt
import seaborn as sns
import numpy as np
import glob,os

from src.sim_double_well_SDE import load_parameters

F_TRAJ = glob.glob("trajectory/umbrella_EnergyBarrier_r*")

for f in F_TRAJ:
    f_np = os.path.basename(f).replace('.txt','.npy')  
    if not os.path.exists(f_np):
        print f_np
        time, pos = np.loadtxt(f).T
        np.save(f_np,pos)

from pymbar import MBAR

params = load_parameters("simulation_setups/umbrella_EnergyBarrier.json")
F_TRAJ = sorted(glob.glob("umbrella_EnergyBarrier_r*.npy"))

def U(x,**kw): return (x**2-1.0)**2 
samples = 1000
simulations = 8
X_k = np.array([np.load(f)[:][:samples] for f in F_TRAJ])
N_k = np.array([np.shape(x)[0] for x in X_k])
#U_k = U(X_k)
U_k = np.zeros((simulations,simulations,samples))

for k in xrange(simulations):
    u = U(X_k[k])
    for l in xrange(simulations):
        U_k[l,k,:] = u


# Estimate free energies from simulation using MBAR.
print "Estimating relative free energies using MBAR, this may take a while..."

# Initialize the MBAR class, determining the free energies.
mbar = MBAR(U_k, N_k, 
            method = 'adaptive',relative_tolerance=1.0e-10,verbose=True)


(Delta_f_ij_estimated, dDelta_f_ij_estimated) = mbar.getFreeEnergyDifferences()

A_kn = U_k
(A_ij, dA_ij) = mbar.computeExpectations(A_kn, output='differences')
print dA_ij
exit()

#print mbar.getFreeEnergyDifferences()
#print Delta_f_ij_estimated.shape
X = np.zeros((simulations,200))
for n in range(simulations):
    X[n] = np.linspace(-2.5,2.5,200)
print mbar.computeExpectations(X)
help(mbar)
exit()

exit()
# Compute error from analytical free energy differences.
Delta_f_ij_error = Delta_f_ij_estimated - Delta_f_ij_analytical

print "Error in free energies is:"
print Delta_f_ij_error

print "Standard deviations away is:"
# mathematical manipulation to avoid dividing by zero errors; we don't care
# about the diagnonals, since they are identically zero.
df_ij_mod = dDelta_f_ij_estimated + numpy.identity(K)
stdevs = numpy.abs(Delta_f_ij_error/df_ij_mod)
for k in range(K):
  stdevs[k,k] = 0;
print stdevs
exit()

# Number of integration bins
integration_bins = 10000


X = np.linspace(-1.5,1.5,integration_bins)


ubounds  = params["umbrella_bounds"]
uwindows = params["umbrella_windows"]
centers  = np.linspace(ubounds[0],ubounds[1],uwindows)
umbrella_strength = params["umbrella_strength"]

beta = 1.0/params["kT"]

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

'''
Y = [np.load(f)[:] for f in F_TRAJ]
for k,y in enumerate(Y):
    print k
    label = "center = %s"%centers[k]
    sns.distplot(y,label=label,hist=False)
plt.legend()
plt.show()
exit()
#for y in Y:
#    print len(y)
#exit()
'''

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




#plt.plot(X,dA_avg,label="Umbrella integration")
#plt.plot(X,U(X),'--',label="Exact")
plt.plot(X[1:],A,label="Umbrella integration")
#plt.xlim(-2,2)
#plt.ylim(0,2)
plt.legend(loc="best")
plt.show()
