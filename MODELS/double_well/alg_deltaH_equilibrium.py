from sim_double_well_SDE import sim_double_well
from metrics_double_well import average_activation_energy

import numpy as np
from scipy.stats.kde import gaussian_kde as KDE
from scipy.integrate import quad

'''
SAMPLING METHOD: Standard
This simulation is the test case, e.g. the standard dynamics with no enchanced sampling.

Computes delta H (the energy barrier) between a double well using overdamped 
Langevin dynamics with an Euler-Maruyama SDE.
'''

args = {"kT":1.5, 
        "friction_coeff":.1,
        "simulation_time": 30,
        "SIM_metric_func":average_activation_energy}

S = sim_double_well(**args)
S.run()

# Compute the exact value for error measurements
Um0, Ub, Um1 = map(S["potential"], [-1,0,1])
exact_avg_activation_energy = np.array([Ub-Um0, Ub-Um1]).mean()

# Compute the error
err = np.abs((np.array(S.traj_metric) - exact_avg_activation_energy))

'''
Plot the results
'''

import pylab as plt
import seaborn as sns

xbounds = 2.5
xp = np.linspace(-xbounds,xbounds,1000)
fig, axes = plt.subplots(2, 2, figsize=(12, 7))

ax = axes[0,0]
ax.set_title(r"Particle trajectory")

T,X = S.traj_t, S.traj_x
ax.plot(T,X)
ax.set_xlim(min(T),max(T))
ax.set_xlabel(r"$t$")
ax.set_ylabel(r"$x$")

ax = axes[0,1]
ax.set_title(r"Potential")
ax.plot(xp,S["potential"](xp))
ax.set_ylim(-1, 4)
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$U(x)$")

# Plot the estimated potential
H = KDE(X)
U = -np.log(H(xp))*S["kT"]
ax.plot(xp,U,color='r',alpha=.75)

ax = axes[1,1]
ax.set_title(r"Observed vs Expected measure")
sns.distplot(X,ax=ax,label=r"$\mu(x)$, invariant measure $e^{-\beta U(x)}$")

invariant_norm = quad(S.invariant_measure, -xbounds, xbounds)[0]
target = S.invariant_measure(xp) / invariant_norm
print invariant_norm

ax.plot(xp,target,'r',alpha=.75,label=r"$H(x)$")
ax.set_ylim(0,1)
ax.legend(loc=0)

ax = axes[1,0]
err_T = S.traj_metric_t
ax.set_title(r"$L_1$ error, $\int \ |\Delta U_{ex} - \Delta U_{approx}| dx$")
ax.plot(err_T,err)
ax.set_xlim(min(T),max(T))
ax.semilogy(err_T,np.zeros(len(err_T)),'k--',alpha=.5)
ax.set_ylim(ymin=0,ymax=1)
plt.tight_layout()
plt.show()
    





    
    
    

