import numpy as np
from scipy.stats.kde import gaussian_kde as KDE
from scipy.integrate import quad

import pylab as plt
import seaborn as sns

def plot_simulation(S, err):

    xbounds = 2.5
    xp = np.linspace(-xbounds,xbounds,1000)
    fig, axes = plt.subplots(2, 2, figsize=(12, 7))

    ax = axes[0,0]
    ax.set_title(r"Particle trajectory")

    T,X = np.array(S.traj_t), np.array(S.traj_x)

    # Plot only the last 100 s of simulation time
    idx = T > (T.max()-25)
        
    ax.plot(T[idx],X[idx])
    ax.set_xlim(min(T[idx]),max(T[idx]))
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

    # Shift the estimated potential down for visual aide
    # (only relative differences matter anyways)
    U += S["potential"](0) + np.log(H(0))*S["kT"]

    ax.plot(xp,U,color='r',alpha=.75)

    ax = axes[1,1]
    ax.set_title(r"Observed vs Expected measure")
    sns.distplot(X,ax=ax,label=r"$\mu(x)$, invariant measure $e^{-\beta U(x)}$")

    invariant_norm = quad(S.invariant_measure, -xbounds, xbounds)[0]
    target = S.invariant_measure(xp) / invariant_norm

    ax.plot(xp,target,'r',alpha=.75,label=r"$H(x)$")
    ax.set_ylim(0,1)
    ax.legend(loc=0)

    ax = axes[1,0]
    err_t = S.traj_metric_t
    ax.set_title(r"$L_1$ error, $\int \ |\Delta U_{ex} - \Delta U_{approx}| dx$")
    ax.plot(err_t,err)
    ax.set_xlim(min(T),max(T))
    ax.set_ylim(ymin=0,ymax=1)
    plt.tight_layout()
    plt.show()

