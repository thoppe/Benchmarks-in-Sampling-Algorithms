from scipy.stats.kde import gaussian_kde as KDE
import numpy as np


def average_activation_energy(S):
    ''' 
    Computes the activation energy (\delta H) averaged across both sides 
    of the double well with minima and maximum at [-1,0,1] respectively. 
    The enthalpy is estimated using kernel density estimation (KDE),
    slow but more accurate then a histogram. 
    '''

    X = S.traj_x
    H = KDE(X)    
    def U_estimated(x): return -np.log(H(x))*S["kT"]
    Em0, Eb, Em1 = map(U_estimated, [-1,0,1])
    activation_energy = np.array([Eb-Em0, Eb-Em1])
    return activation_energy.mean()
