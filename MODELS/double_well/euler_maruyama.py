import numpy as np
from numpy.random import normal

'''
DRAFT for:
Overdamped Langevin dynamics in a double well potential
Euler-Maruyama SDE
Author: Travis Hoppe
'''

class SDE_euler_maruyama(dict):
    ''' 
    Integrates a stochastic differential equation of the form:
        dx(t) = f(x,t)*dt + g(x,t)*dW
    here W is a Wigner process. Euler-Maruyama advances the SDE by the scheme:
        x(t+dt) = f(x,t)*dt + g(x,t)*normal()*np.sqrt(dt)
    '''
   

    def __init__(self, **simulation_args):

        # Set the defaults
        self["x0"] = 0.0
        self["t0"] = 0.0
        
        self["kT"] = 1.0
        self["dt"] = 0.001
        
        self["f"] = self["g"] = None

        # Update the system parameters if passed
        self.update( simulation_args )

        self.x = self["x0"]
        self.t = self["t0"]

        self.n = 0 # Number of function calls

    def step(self,**kw):

        # Make sure functions have been defined
        if not (self["f"] or self["g"]):
            msg = "Must define functions f and g"
            raise SyntaxError(msg)

        x,t = self.x, self.t
        dt  = self["dt"]

        dx  = self["f"](x,t,**self)*dt
        dx += self["g"](x,t,**self)*normal()*np.sqrt(dt)

        self.x += dx
        self.t += dt
        self.n += 1


''' 
Define the Langevin overdamped double-well potential system:
'''

from scipy.misc import derivative
from scipy.stats.kde import gaussian_kde as KDE
from scipy.integrate import quad, trapz


def double_well(x, **kwargs):
    return (x**2-1.0)**2 

def force(x, t, **kwargs):
    return -derivative(kwargs["potential"], x, dx=.001)/kwargs["friction_xi"]

def brownian_motion(x, t, **kwargs):
    diffusion_coeff = kwargs["kT"]/kwargs["friction_xi"]
    return np.sqrt(2*diffusion_coeff)

args = {
    "f":force, 
    "g":brownian_motion,
    "potential":double_well,
    "friction_xi": 1.0,
    }


S = SDE_euler_maruyama(**args)

'''
Determine the metric and invariant measure (for errors)
'''

def invariant_measure(x,S):
    diffusion_coeff    = S["kT"]/S["friction_xi"]
    return np.exp(-S["potential"](x)/diffusion_coeff)


xbounds = 2.5
xp = np.linspace(-xbounds,xbounds,1000)
invariant_norm = quad(invariant_measure, -xbounds, xbounds, args=(S,))[0]
target = invariant_measure(xp,S) / invariant_norm

'''
Run the simulation:
'''

error_check = 1000
simulation_time = 100.0

X,T = [S.x,],[S.t,]
err, err_T = [],[]

while S.t < simulation_time:
    S.step()
    X.append(S.x)
    T.append(S.t)

    if (S.n%error_check==0):
        H = KDE(X)
        err_T.append(T[-1])

        def estimated_pot(x): return -np.log(H(x))*S["kT"]
        Em0, Eb, Em1 = map(estimated_pot, [-1,0,1])
        Um0, Ub, Um1 = map(args["potential"], [-1,0,1])

        dE = np.array([Eb-Em0, Eb-Em1])
        dU = np.array([Ub-Um0, Ub-Um1])
        err_term = np.abs(dE - dU).mean()
        err.append(err_term)

        print S.t, err[-1]

'''
Plot the results
'''

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
U = -np.log(H(xp))*S["kT"]
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
    





    
    
    

