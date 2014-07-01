# Benchmarks in Sampling Algorithms
Originally proposed at the 2014 Telluride Science Research Center: Advances in Sampling algoirthms.

The goal of this project is to develop a suite of cannonical test problems where various sampling algorithms can be compaired quantatively.
Since sampling alorithms are not always directly comparable, we seperate the problem sets into several groups.

### Simple acceptence ratios on Energy: 

+ Double well in 1D, 2D, 3D
+ Ising/Potts models on fixed lattice sizes, either 1D or 2D periodic boundary conditions
+ LJ gas, with a fixed number of atoms (48 I think) there are two almost equal energy minima at completely different configurations
+ On/Off lattice polymers with simple potentials e.g. Go models or MJ interactions

It is important to determine a metric for the right answer. These simple system are proposed because they generally have fully enumerable states or an analytical solution. That said, a few reasonable benchmarks are:

+ Quality of the density of states
+ Time to first find first global energy minima
+ Tunneling time from two extremal states
+ Measuring of the free energy barrier (when relevant)