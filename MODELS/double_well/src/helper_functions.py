import numpy as np
import json, os

import metrics_double_well as metrics

def load_parameters(f_json):
    '''
    Helper function that loads a JSON file with the system parameters in it.
    Checks if the metric function is a valid known function.
    Also creates the directory "results" if it doesn't exist.

    Parameters
    ----------
    f_json : string
        Filename of the input json file

    Raises
    -------
    ValueError
        For problems loading the json file.

    KeyError
        For problems loading the metric.

    Returns
    -------
    params : dict
        A dictionary with the paramaters loaded.
    '''

    with open(f_json) as FIN:
        try:
            params = json.loads(FIN.read())
        except Exception as ex:
            err_msg = "Problem with json file: {} {}".format(f_json, ex)
            raise ValueError(err_msg)

        # Turn the loaded string SIM_metric_func into a function
        try:
            func_name = "{}.{}".format("metrics", params["SIM_metric_func"])
            params["SIM_metric_func"] = eval(func_name)
        except Exception as ex:
            err_msg = "Problem with metric: {} {}"
            err_msg =  err_msg.format(func_name, ex)
            raise KeyError(err_msg)

    target_dir = "results"
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    return params

###########################################################################

def save_results(f_results, time_steps, epsilon):
    # Save the result data to f_results
    data = np.array([time_steps,epsilon]).T
    return np.savetxt(f_results,data)

def load_results(S):
    # Loads the results the file S["f_results"]
    # Returns epsilon_time_step, epsilon
    return np.loadtxt(S["f_results"],unpack=True)

def load_trajectory(S):
    # Loads the trajectory in the file S["f_trajectory"]
    # Returns T, X
    return np.loadtxt(S["f_trajectory"],unpack=True)

def iterate_trajectory(S):
    ''' 
    Iterates through the file S["f_trajectory"] in steps of S["metric_check"]
    '''
    T,X   = load_trajectory(S)
    d_idx = S["metric_check"]

    for idx in xrange(0, X.size, d_idx):
        x_chunk = X[idx:idx+d_idx]
        t_chunk = T[idx:idx+d_idx]
        if x_chunk.size == d_idx:
            yield t_chunk,x_chunk

###########################################################################

def compute_activation_error(S):
    # Compute the exact value for error measurements
    Um0, Ub, Um1 = map(S["potential"], [-1,0,1])
    exact_activation_energy = np.array([Ub-Um0, Ub-Um1])

    # Measure the barrier height from the trajectory
    H = metrics.activation_energy
    estimated_activation_energy, time_steps = H(**S)

    # Compute the error, summed over both wells
    epsilon = np.abs(estimated_activation_energy-exact_activation_energy)
    epsilon = epsilon.sum(axis=1)

    return epsilon, time_steps

