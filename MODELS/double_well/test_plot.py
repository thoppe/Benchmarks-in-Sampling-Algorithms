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
X = np.linspace(-3,3,1000)
Y = [np.load(f)[50000:] for f in F_TRAJ]
centers =  [-2,-1.5,-1,-0.5,0.0,0.5,1.5,2.0]
umbrella_strength = 2.0

# Try WHAM!
import uwham
Nk = [len(y) for y in Y]
print Nk



print "HERE!"
exit()


def weights(X, Y):
    mu,sigma = Y.mean(), np.sqrt(Y.var())
    w = np.exp(-0.5*(((X-mu)/sigma)**2))
    return w*(sigma*np.sqrt(2*np.pi))

def diff_free_energy(X, Y, k, x0, kT=1.0):
    mu,variance = Y.mean(), Y.var()
    A = kT*((X-mu)/variance) - k*(X-x0)
    return A

Y = [np.load(f)[:] for f in F_TRAJ]
#for k,y in enumerate(Y):
#    print k
#    label = "center = %s"%centers[k]
#    sns.distplot(y,label=label,hist=False)
#plt.legend()
#plt.show()
#exit()
#for y in Y:
#    print len(y)
#exit()


P = np.array([weights(X,y) for y in Y])
#P /= P.sum(axis=0)

dA = [diff_free_energy(X,y,
                       k=umbrella_strength,
                       x0=c) for c,y in zip(centers, Y)]
dA_avg = (dA*P).sum(axis=0)

from scipy.integrate import cumtrapz, simps
#A = cumtrapz(dA_avg,X)
#A-= A.min()
#print A
#exit()
'''
from scipy.integrate import cumtrapz, simps
F_Y, F_X = [],[]
k = 20
for x,y in zip(zip(X,X[k:]),zip(dA_avg,dA_avg[k:])):
    #x_grid = np.array([x0,x1])
    #y_grid = np.array([a0,a1])
    F_X.append(np.array(x).mean())
    F_Y.append(simps(y,x))
'''
#plt.plot(X,dA_avg)

DP = [np.log(p)/(X[0]-X[1]) for p in P]

for y in DP:
    plt.plot(X,y)


#plt.ylim(-0.1, 1.2)
#plt.plot(A)
plt.show()

exit()

#A = cumtrapz(dA_avg,X,initial=0)
#A -= A.min()

#plt.plot(X,A)
##plt.ylim(-0.1, 1.2)
#plt.show()
exit()
#print dA
#exit()
#print P
#exit()
for p in P: plt.plot(X,p)
#for p in dA: plt.plot(X,p)
plt.show()
exit()

plt.show()
exit()
#    sns.distplot(pos)

