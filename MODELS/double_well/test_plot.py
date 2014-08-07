import pylab as plt
import seaborn as sns
import numpy as np
import glob,os


F_TRAJ = glob.glob("trajectory/umbrella_EnergyBarrier_r*")
for f in F_TRAJ:
    f_np = os.path.basename(f).replace('.txt','.npy')  
    if not os.path.exists(f_np):
        print f_np
        time, pos = np.loadtxt(f).T
        np.save(f_np,pos)

F_TRAJ = sorted(glob.glob("umbrella_EnergyBarrier_r*.npy"))
centers =  [-1,0,1]
umbrella_strength = 2.0
X = np.linspace(-2,2, 500)

def weights(X, Y):
    mu,sigma = Y.mean(), np.sqrt(Y.var())
    w = np.exp(-0.5*(((X-mu)/sigma)**2))
    return w*(sigma*np.sqrt(2*np.pi))

def diff_free_energy(X, Y, k, x0, kT=1.0):
    mu,variance = Y.mean(), Y.var()
    A = kT*((X-mu)/variance) - k*(X-x0)
    return A

Y = [np.load(f) for f in F_TRAJ]
P = np.array([weights(X,y) for y in Y])
P /= P.sum(axis=0)

dA = [diff_free_energy(X,y,
                       k=umbrella_strength,
                       x0=c) for c,y in zip(centers, Y)]
dA_avg = (dA*P).sum(axis=0)

from scipy.integrate import cumtrapz
A = cumtrapz(dA_avg,X,initial=0)
A -= A.min()

plt.plot(X,A)
plt.ylim(-0.1, 1.2)

#print dA
#exit()
#print P
#exit()
#for p in P: plt.plot(X,p)
#for p in dA: plt.plot(X,p)

plt.show()
exit()
#    sns.distplot(pos)

