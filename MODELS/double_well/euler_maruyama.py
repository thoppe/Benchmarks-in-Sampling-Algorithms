import numpy as np
from numpy.random import normal
from scipy.misc import derivative
from scipy.stats.kde import gaussian_kde as KDE
from scipy.integrate import quad, trapz

'''
DRAFT for:
Overdamped Langevin dynamics in a double well potential
Euler-Maruyama SDE
Author: Travis Hoppe
'''

# Simultation paramters
simulation_time = 1000.0
dt    = .005

# Initial position
x0 = 0.0
t0 = 0.0

# Extent to calculation the invariant measure over (for errors)
xbounds = 2.5
xp = np.linspace(-xbounds,xbounds,1000)
error_check = 1000

def double_well(x, **kwargs):
    return (x**2-1.0)**2
    
args = {"kT": 1.0,
        "friction_xi": 1.0,
        "potential":double_well}

def force(x, **kwargs):
    return -derivative(kwargs["potential"], x, dx=.001)/kwargs["friction_xi"]

def brownian_motion(x, **kwargs):
    diffusion_coeff = kwargs["kT"]/kwargs["friction_xi"]
    return np.sqrt(2*diffusion_coeff)

def EULER_MARUYAMA_SDE(x, dt, f, g,**kw):
    return f(x,**kw)*dt + g(x,**kw)*normal()*np.sqrt(dt)

def invariant_measure(x):
    diffusion_coeff    = args["kT"]/args["friction_xi"]
    return np.exp(-args["potential"](x)/diffusion_coeff)

invariant_norm = quad(invariant_measure, -xbounds, xbounds)[0]
target = invariant_measure(xp) / invariant_norm

SDE_f = force
SDE_g = brownian_motion

X,T = [x0,],[t0,]
err, err_T = [],[]

while T[-1] < simulation_time:

    x,t = X[-1],T[-1]
    dx  = EULER_MARUYAMA_SDE(x, dt, SDE_f, SDE_g, **args)
    X.append(x+dx)
    T.append(t+dt)

    if (len(T)-2)%error_check==0:
        H = KDE(X)
        err_T.append(T[-1])

        def estimated_pot(x): return -np.log(H(x))*args["kT"]
        Em0, Eb, Em1 = map(estimated_pot, [-1,0,1])
        Um0, Ub, Um1 = map(args["potential"], [-1,0,1])

        dE = np.array([Eb-Em0, Eb-Em1])
        dU = np.array([Ub-Um0, Ub-Um1])
        err_term = np.abs(dE - dU).mean()
        err.append(err_term)

        print t, err[-1]

# Plot the results

import pylab as plt
import seaborn as sns
fig, axes = plt.subplots(2, 2, figsize=(12, 7))

ax = axes[0,0]
ax.set_title(r"Particle trajectory")
ax.plot(T,X)
ax.set_xlim(min(T),max(T))
ax.set_xlabel(r"$t$")
ax.set_ylabel(r"$x$")

ax = axes[0,1]
ax.set_title(r"Potential")
ax.plot(xp,args["potential"](xp))
ax.set_ylim(-1, 4)
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$U(x)$")

# Plot the estimated potential
H = KDE(X)
U = -np.log(H(xp))*args["kT"]
ax.plot(xp,U,color='r',alpha=.75)

ax = axes[1,1]
ax.set_title(r"Observed vs Expected measure")
sns.distplot(X,ax=ax,label=r"$\mu(x)$, invariant measure $e^{-\beta U(x)}$")
ax.plot(xp,target,'r',alpha=.75,label=r"$H(x)$")
ax.set_ylim(0,1)
ax.legend(loc=0)

ax = axes[1,0]
ax.set_title(r"$L_1$ error, $\int \ |\Delta U_{ex} - \Delta U_{approx}| dx$")
ax.plot(err_T,err)
ax.set_xlim(min(T),max(T))
ax.semilogy(err_T,np.zeros(len(err_T)),'k--',alpha=.5)
ax.set_ylim(ymin=0,ymax=1)
plt.tight_layout()
plt.show()
    





    
    
    

