import numpy as np
from numpy.random import normal
from scipy.misc import derivative
from scipy.stats.kde import gaussian_kde as KDE

'''
DRAFT for:
Euler-Maruyama SDE (high-friction limit) for double well potential
[units may not be correct yet]
Author: Travis Hoppe
'''


total_timesteps = 100000
sigma = .75
dt    = .010
x0 = 0.0
t0 = 0.0

def potential(x, **kwargs):
    return (x**2-1.0)**2

def force(x, **kwargs):
    return -derivative(potential, x, dx=.001)

def constant_noise(x, **kwargs):
    return kwargs["sigma"]

def EULER_MARUYAMA_SDE(x, dt, f, g,**kw):
    return f(x,**kw)*dt + g(x,**kw)*normal()*np.sqrt(dt)

SDE_f = force
SDE_g = constant_noise

X,T = [x0,],[t0,]
xp = np.linspace(-2,2,1000)
err, err_T = [],[]

for k in xrange(total_timesteps):

    x,t = X[-1],T[-1]
    dx  = EULER_MARUYAMA_SDE(x, dt, SDE_f, SDE_g, sigma=sigma)
    X.append(x+dx)
    T.append(t+dt)

    if k and k%1000==0:
        print k
        H = KDE(X)        
        err.append( (H(-1)/H(1))[0] )
        err_T.append(T[-1])

import pylab as plt
import seaborn as sns
fig, axes = plt.subplots(2, 2, figsize=(12, 7))
ax = axes[0,0]
ax.plot(T,X)

ax = axes[0,1]
ax.plot(xp,potential(xp))

ax = axes[1,0]
sns.distplot(X,ax=ax)

#H = KDE(X)(xp)
#ax.plot(xp, H,color='r')

ax = axes[1,1]
ax.plot(err_T,err)
ax.plot(err_T,np.ones(len(err_T)),'k--',alpha=.5)
ax.set_ylim(0.5, 1.5)

plt.show()
    





    
    
    

