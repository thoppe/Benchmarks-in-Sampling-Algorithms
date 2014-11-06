import logging, os
import numpy as np
from numpy.random import normal
import metrics_double_well as metrics

'''
Overdamped Langevin dynamics in a double well potential
Euler-Maruyama SDE
Author: Travis Hoppe

The classes build off each other, to run a simulation define a metric function,
set some parameters and run the simulation by passing a dictionary:

args = {"kT":1.5, 
        "friction_coeff":.1,
        "simulation_time": 30,
        "SIM_metric_func":average_activation_energy}

S = sim_double_well(**args)
S.run()
'''


##############################################################################
# http://stackoverflow.com/a/3653049/249341
##############################################################################


import inspect
import functools

def autoargs(*include,**kwargs):   
    def _autoargs(func):
        attrs,varargs,varkw,defaults=inspect.getargspec(func)

        @functools.wraps(func)
        def wrapper(self,*args,**kwargs):

            # Do not set default values if already given a value
            # by a child-class; handle default values
            if defaults:
                for attr,val in zip(reversed(attrs),reversed(defaults)):
                    if attr not in self.__dict__:
                        setattr(self, attr, val)

            # handle varkw
            if kwargs:
                for key,val in kwargs.iteritems():
                    if key in attrs:
                        setattr(self,key,val)                       

            return func(self,*args,**kwargs)
        return wrapper
    return _autoargs


##############################################################################

class SDE_euler_maruyama(object):
    ''' 
    Integrates a stochastic differential equation of the form:
        dx(t) = f(x,t)*dt + g(x,t)*dW
    here W is a Wigner process. Euler-Maruyama advances the SDE by the scheme:
        x(t+dt) = x(t) + f(x,t)*dt + g(x,t)*normal()*np.sqrt(dt)

    Parameters
    ----------
    dt : Integration timestep  [default = 0.001]
    xi : Initial x position    [default = 0.0]
    ti : Initial time          [default = 0.0]

    SDE_f : f(x,t) : Must be set before step is called.
    SDE_g : g(x,t) : Must be set before step is called.

    Internal Parameters
    ----------
    n  : Number of integration calls [default = 0]

    Raises
    -------
    SyntaxError
        If functions f and g are not defined when running.
    '''
    @autoargs()
    def __init__(self, 
                 dt = 0.001,
                 xi = 0.0,
                 ti = 0.0,
                 n  = 0,
                 SDE_f = None,
                 SDE_g = None,
                 **simulation_args):
        super(SDE_euler_maruyama, self).__init__()

    def check(self):
        # Make sure functions have been defined
        if not (self.SDE_f or self.SDE_g):
            msg = "Must define functions SDE_f and SDE_g"
            raise SyntaxError(msg)

    def step(self):
        self.check()

        x,t,dt = self.xi, self.ti, self.dt

        # Euler-Maruyama scheme, first compute f(x,t)*dt
        dx  = self.SDE_f(x,t)*dt

        # then compute g(x,t)*normal()*np.sqrt(dt)
        dx += self.SDE_g(x,t)*normal()*np.sqrt(dt)

        # now adjust the position by this amount
        self.xi += dx
        self.ti += dt
        self.n  += 1

#params = {"dt":2}
#S = SDE_euler_maruyama(**params)
#print S.__dict__
#exit()

class overdamped_langevin(SDE_euler_maruyama):
    ''' 
    Defines a generic one-dimensional Langevin overdamped system
    integrated by Euler-Maruyama.

    This sets SDE_g to brownian motion defined by the friction_coeff.

    Parameters
    ----------
    all parameters defined by SDE_euler_maruyama AND

    friction_coeff : zeta      [default = 1.0]
    kT : System temperature    [default = 1.0]
    '''

    def brownian_motion(self, x, t,**kw):
        return np.sqrt(2*self.kT/self.friction_coeff)

    @autoargs()
    def __init__(self, 
                 kT = 1.0,
                 friction_coeff = 1.0,
                 **simulation_args):

        super(overdamped_langevin, self).__init__(**simulation_args)

        # The overdamped langevin experiences constant brownian motion
        self.SDE_g = self.brownian_motion

#params = {"dt":2}
#S = overdamped_langevin(**params)
#print S.__dict__
#print S.SDE_g(2,2)
#exit()

class double_well(overdamped_langevin):
    ''' 
    Defines a symmetric double well of height 1kT, under a 
    Langevin overdamped system integrated by Euler-Maruyama.

    This sets SDE_f to the gradient of the double well potential.

    Parameters
    ----------
    all parameters defined by SDE_euler_maruyama, overdamped_langevin

    bias_potential             [default: returns 0.0]
    bias_force                 [default: returns 0.0]
    '''

    def invariant_measure(self,x):
        return np.exp(-self.U(x)/self.kT)

    def U(self,x,t=None,**kw):
        E = (x**2-1.0)**2 
        E_bias = self.bias_potential(x,t)
        return E + E_bias
    
    def force(self,x,t=None,**kw):
        dU      = -4*x*(x**2-1)
        dU_bias = self.bias_force(x,t)
        return (dU+dU_bias)/self.friction_coeff
        
    def __zero_function(x,y,**kw): 
        return 0.0

    @autoargs()
    def __init__(self, 
                 bias_potential = __zero_function,
                 bias_force     = __zero_function,
                 **simulation_args):

        super(double_well, self).__init__(**simulation_args)

        self.potential = self.U
        self.SDE_f     = self.force


#params = {"dt":2}
#S = double_well(**params)
#print S.invariant_measure(2)
#print S.SDE_g(2,2)
#exit()

class sim_double_well(double_well):
    ''' 
    Runs the simulation for the double well, 
    saves the output to file.

    Parameters
    ----------
    all parameters defined by SDE_euler_maruyama, overdamped_langevin and

    simulation_time : total integrated time     [default = 300.0]
    
    f_trajectory    : output file to save trajectory data 
                      [default=None]
    simulation_time : total time to integrate if self.run is called

    Internal Parameters
    ----------
    simulation_step : Current count of simulation steps [default = 0]

    '''

    @autoargs()
    def __init__(self, 
                 simulation_time = 300.0,
                 simulation_step = 0,
                 f_trajectory    = "output.txt",
                 **simulation_args):

        super(sim_double_well, self).__init__(**simulation_args)

        # Create the file object, if not possible skip
        self.FOUT = None
        try:
            self.FOUT = open(self.f_trajectory,'w')
        except:
            pass

    # __enter__ and __exit__ allow the class to be used as a context manager
    def __enter__(self):  
        return self

    def __exit__ (self, exception_type, exception_value, traceback):  
        self.close()

    def record_traj(self):
        if self.FOUT:
            out_str = "{:.5f} {:.7f}\n"
            self.FOUT.write(out_str.format(self.ti,self.xi))

    def is_complete(self):
        return self.ti > self.simulation_time

    def run(self, fixed_time=None, record=True):

        run_counter = 0
        if record: self.record_traj()

        start_time = self.ti

        # If a fixed time is not set, run the length of the simulation
        if fixed_time == None:
            fixed_time = self.simulation_time

        while self.ti - start_time < fixed_time:
            
            self.step()

            if record: self.record_traj()
            self.simulation_step += 1

            run_counter += 1

            if run_counter and run_counter%1000==0:
                if record:
                    logging.info("Simulation time %f" % self.ti)
                if not record:
                    logging.info("Warmup time     %f" % self.ti)

    def close(self):
        if self.FOUT:
            self.FOUT.close()

#params = {"dt":2}
#S = sim_double_well(**params)
#print S.invariant_measure(2)
#print S.SDE_g(2,2)
#exit()
