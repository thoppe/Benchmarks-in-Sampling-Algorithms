# Double well

+ **Model**: [Overdamped Langevin dynamics](http://en.wikipedia.org/wiki/Langevin_dynamics) in a double well potential
+ **Integrator**: [Euler-Maruyama SDE](http://en.wikipedia.org/wiki/Euler-Maruyama)
+ **Author**: [Travis Hoppe](https://github.com/thoppe)

This is the initial setup of the toy problem of a one-dimensional symmetric double well potential, with a barrier height of 1kT. The potential energy is:

<p align="center" class="mdequation"><img src=".equations/8df0e476ded358dc7b2784f8d8ab5861a383631a73dd32f678a8dcd46c07a3db.png" alt="$U(x) = (x^2 - 1)^2$" /></p>

The motion is overdamped and stochastic, hence the instantaneous momentum is simply a combination of Brownian motion and the underlying potential.
The system evolves according to the stochastic differential equation:

<p align="center" class="mdequation"><img src=".equations/b2304f4808aa680f5d85d94b293e3cc91bea3288231acb61cbcf6823d983adf2.png" alt="$ \dot{x}(t) = - \nabla U(x)/\zeta^{-1} + \sqrt{2 kT \zeta^{-1}} W(t)  $" /></p>

where zeta is the frictional coefficient times the mass and W is a delta-correlated stationary Gaussian process with zero-mean (simulating random thermal motion). 
Given enough time, the trajectory of the particle samples the the invariant measure

<p align="center" class="mdequation"><img src=".equations/7083ff27c4483a6ecaa6326ce60b0df098270a31d358922510406337b8f46a22.png" alt="$\mu(x) = e ^{-U(x)/kT} \mathcal{Z}^{-1}$" /></p>

## Metric: Energy Barrier

The energy barrier is estimated on both sides and compared to the exact value.
The error term is the L1 average of these differences.

<p align="center" class="mdequation"><img src=".equations/7aa8e55b6db2940e1365f4811d7350f072c259abe1f681bf86a7adc2766ae28d.png" alt="$ \epsilon = \frac{1}{2} (|\Delta V(-1) - 1| + |\Delta V(1) - 1|)$" /></p>

Here V(x) is the estimated potential and \Delta V(x) measures the difference from x to the top of the well.
Given a histogram H(x), that records the number of visits to each state one can estimate the potential as

<p align="center" class="mdequation"><img src=".equations/85b3044977c3d55392496cd69131615af66aca3d41af7c3d9981ddd0c6447e8e.png" alt="$V(x) \approx -kT \log( H(x) ) $" /></p>

The histogram can be a simple bin, or alternatively one can use a kernel density estimation on the entire trajectory.

#### Results: Energy Barrier

The results of the sampling algorithms are shown below:

![](figures/convergence_EnergyBarrier.png)

**Sampling Algorithm**: None

The simulation shows a naive way of calculating the energy barrier, simply let the system evolve. 
Shown below is a sample trajectory, the estimated potential, the error and the observed versus expected visits to each position. 
While the estimated potential has a large absolute error, the estimated energy difference between the two wells converges quickly since the potential is simple.

![](figures/example_traj.png)

The simulation can be repeated by running:

    python alg_EnergyBarrier_equilibrium.py simulation_setups/example_EnergyBarrier.json

Where the configuration file [example_EnergyBarrier.json](simulation_setups/example_EnergyBarrier.json) is given by (TO DO: Comment on parameters):

```JSON
{
    "kT": 1.0,
    "friction_coeff" : 0.1,

    "dt" : 0.0001,

    "simulation_time": 2000,
    "metric_check"   : 10000,
 
    "SIM_metric_func": "average_activation_energy",
    "f_results": "results/example_EnergyBarrier_r{replica_n}.txt",

    "show_plot" : false
}
```

**Sampling Algorithm**: Replica Exchange

The simulation can be repeated by running:

    python alg_EnergyBarrier_replicaEx.py simulation_setups/replicaEx_EnergyBarrier.json

In addition to the parameters set by the simple sampling algorithm, the following options are accepted:

```JSON
{
    "kT_list" : [0.8,0.9,1,1.1,1.2],
    "exchange_steps" : 30000,
}
```





