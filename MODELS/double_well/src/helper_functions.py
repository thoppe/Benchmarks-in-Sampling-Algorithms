import numpy as np
import metrics_double_well as metrics
import json, os, argparse, logging

###########################################################################

def startup_simulation(desc="",
                       logging_level=logging.INFO,
                       format_filenames=True):
    '''
    Helper function that runs command-line parsing, 
    and sets the logging level.

    Parameters
    ----------
    desc : string
        Program description

    Returns
    -------
    params : dict
        A dictionary with the paramaters loaded.
    '''

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('parameter_file_json')
    parser.add_argument('--replica_n', type=int, default=0,
                        help="Formats this number into the results file")
    cargs = vars(parser.parse_args())
    
    # Start the logger
    logging.root.setLevel(logging_level)

    # Load the simulation parameters
    params = load_parameters(cargs["parameter_file_json"])

    if format_filenames:
        params["f_trajectory"] = params["f_trajectory"].format(**cargs)
        params["f_results"]    = params["f_results"].format(**cargs)

    return params


def finalize_simulation(S, metric_function):
    '''
    Helper function finishes a simulation.

    Parameters
    ----------

    S              : sim_double_well
        Completed simulation

    metric_funcion : function
        Function to compute the metric over, must return
        (time_steps, epsilon)

    Returns
    -------
    params : dict
        A dictionary with the paramaters loaded.
    '''
    # Close any open files
    S.close()

    time_steps, epsilon = metric_function(S)

    # Save the results to file
    save_results(S["f_results"], time_steps, epsilon)

    # Plot the results if asked
    if "show_plot" in S and S["show_plot"]:
        from plots_double_well import plot_simulation
        plot_simulation(S)
   


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
