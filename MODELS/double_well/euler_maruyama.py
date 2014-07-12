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
        self["kT"] = 1.0
        self["dt"] = 0.001

        # The current position and time
        self["xi"] = 0.0
        self["ti"] = 0.0

        # Number of integration calls
        self["n"]  = 0
        
        self["SDE_f"] = self["SDE_g"] = None

        # Update the system parameters if passed
        self.update( simulation_args )

    def check(self):
        # Make sure functions have been defined
        if not (self["SDE_f"] or self["SDE_g"]):
            msg = "Must define functions f and g"
            raise SyntaxError(msg)

    def step(self):
        self.check()

        x,t,dt = self["xi"], self["ti"], self["dt"]

        dx  = self["SDE_f"](x,t,**self)*dt
        dx += self["SDE_g"](x,t,**self)*normal()*np.sqrt(dt)

        self["xi"] += dx
        self["ti"] += dt
        self["n"]  += 1

class overdamped_langevin(SDE_euler_maruyama):
    ''' 
    Defines a generic one-dimensional Langevin overdamped system
    integrated by Euler-Maruyama.
    '''

    def brownian_motion(self, x, t,**kw):
        return np.sqrt(2*self["kT"]/self["friction_coeff"])

    def __init__(self, **simulation_args):
        super(overdamped_langevin, self).__init__(**simulation_args)

        # Set the defaults
        self["friction_coeff"] = 1.0
        self["kT"] = 1.0
       
        # Update the system parameters if passed
        self.update( simulation_args )

        # The overdamped langevin experiences constant brownian motion
        self["SDE_g"] = self.brownian_motion

class double_well(overdamped_langevin):
    ''' 
    Defines a symmetric double well of height 1kT, under a 
    Langevin overdamped system integrated by Euler-Maruyama.
    '''

    def invariant_measure(self,x):
        return np.exp(-self.U(x)/self['kT'])

    def U(self,x,**kw):
        return (x**2-1.0)**2 
    
    def force(self,x,t,**kw):
        dU = 4*x*(x**2-1)
        return -dU/self["friction_coeff"]

    def __init__(self, **simulation_args):
        super(double_well, self).__init__(**simulation_args)

        self["potential"] = self.U
        self["SDE_f"] = self.force   


from scipy.stats.kde import gaussian_kde as KDE
from scipy.integrate import quad

args = {"kT":1.5, "friction_coeff":.1}
S = double_well(**args)

'''
Run the simulation:
'''


error_check = 1000
simulation_time = 300.0

X,T = [S["xi"],],[S["ti"],]
err, err_T = [],[]

while S["ti"] < simulation_time:
    S.step()
    X.append(S["xi"])
    T.append(S["ti"])

    if (S["n"]%error_check==0):
        H = KDE(X)
        err_T.append(T[-1])

        def estimated_pot(x): return -np.log(H(x))*S["kT"]
        Em0, Eb, Em1 = map(estimated_pot, [-1,0,1])
        Um0, Ub, Um1 = map(S["potential"], [-1,0,1])

        dE = np.array([Eb-Em0, Eb-Em1])
        dU = np.array([Ub-Um0, Ub-Um1])
        err_term = np.abs(dE - dU).mean()
        err.append(err_term)

        print S["ti"], err[-1]


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
ax.set_title(r"$L_1$ error, $\int \ |\Delta U_{ex} - \Delta U_{approx}| dx$")
ax.plot(err_T,err)
ax.set_xlim(min(T),max(T))
ax.semilogy(err_T,np.zeros(len(err_T)),'k--',alpha=.5)
ax.set_ylim(ymin=0,ymax=1)
plt.tight_layout()
plt.show()
    





    
    
    

