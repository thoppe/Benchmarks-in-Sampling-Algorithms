import numpy as np
import logging
from helper_functions import iterate_trajectory

def activation_energy(**S):
    ''' 
    Computes the activation energy (\delta U) over both sides  
    of the double well with minima and maximum at [-1,0,1] respectively. 
    A histogram H, is built from the saved trajectory and the potential energy
    is estimated by U = -np.log(H)*kT. 

    Parameters
    ----------
    f_trajectory   : filename of input trajectory stored as (t,x)
    histogram_bins : number of bins in the histogram
    histogram_min  : left side of the histogram
    histogram_max  : right side of the histogram
    
    Returns
    -------
    \delta H : np.array([measurements, 2])
        Measurements of the left and right side of the energy barrier at
        intervals of S["metric_check"]
    '''

    logging.info("Computing the activation energy for %s"%S["f_trajectory"])

    H,bins = np.histogram([],bins=S["histogram_bins"],
                          range=(S["histogram_min"],S["histogram_max"]))
    bins = np.linspace(S["histogram_min"],
                       S["histogram_max"],
                       S["histogram_bins"]+1)
    H    = np.ones(bins.size)

    def close_index(bins, value):
        return (np.abs(bins-value)).argmin()
    
    def U_estimated(idx): return -np.log(H[idx])*S["kT"]

    target_E = [-1,0,1]
    U_index  = [close_index(bins, E) for E in target_E]

    all_E = []
    all_time_step = []
    time_step = 0

    for t,x in iterate_trajectory(S):
        time_step += t.size
        H += np.bincount(np.digitize(x, bins),minlength=H.size)
        Em0, Eb, Em1 = map(U_estimated, U_index)
        activation_energy = np.array([Eb-Em0, Eb-Em1])

        all_E.append( activation_energy )
        all_time_step.append( time_step )

    return np.array(all_E), np.array(all_time_step)


###########################################################################


def compute_activation_error(S, well_locations = [-1,0,1]):
    # Compute the exact value for error measurements
    Um0, Ub, Um1 = map(S["potential"], well_locations)
    exact_activation_energy = np.array([Ub-Um0, Ub-Um1])

    # Measure the barrier height from the trajectory
    estimated_activation_energy, time_steps = activation_energy(**S)

    # Compute the error, summed over both wells
    epsilon = np.abs(estimated_activation_energy-exact_activation_energy)
    epsilon = epsilon.sum(axis=1)

    return time_steps, epsilon

