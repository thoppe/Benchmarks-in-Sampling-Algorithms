# Benchmarks in Sampling Algorithms
Originally proposed at the 2014 Telluride Science Research Center: [_Advances in Enhanced Sampling Algorithms_](https://www.telluridescience.org/meetings/workshop-details?wid=422).

The goal of this project is to develop a suite of canonical test problems where various sampling algorithms can be compared quantitatively.
Since sampling algorithms are not always directly comparable, we separate the problem sets into several groups.

## Proposed model systems

### Simple potentials

+ Double well in 1D, 2D, 3D
+ On/Off lattice polymers with simple potentials e.g. Go models or MJ interactions
+ Muller-Brown like: polynomial or exponential version

### Biological systems, Molecular Dynamics

+ Alanine dipeptide with and without solvent
+ Met-enkephalin 
+ Trp (1L2Y) cage in water
+ Deca alanine 
+ RNA tetramer
+ Cyclodextrin (ligand binding), host-guest systems
+ Lysozyme
+ BPTI (bovine pancreatic trypsin inhibitor)

### Glassy systems

+ Ising/Potts models on fixed lattice sizes, either 1D or 2D periodic boundary conditions
+ LJ gas, with a fixed number of atoms, particular choices of N have almost equal energy minima at completely different configurations 

### Liquid simulations

+ Square-well, LJ phase separations
+ bulk water (TIP3)

-------------------------------------------------------------

It is important to determine a metric for the right answer, different algorithms tend to concentrate on solving a particular problems.
The simple system are proposed because they generally have fully enumerable states or an analytical solution.
The larger systems, especially the biological ones, may need to be compared to each other or an agreed gold-standard (ANTON).
The "unit" of measure is the number of potentials calls, since it is assumed that this is an order of magnitude more computationally expensive than the acceptance function. 

Suggested metrics:

+ Quality of the density of states
+ Time to first find first global energy minima
+ Tunneling time from two extreme states
+ Measuring of the free energy barrier (when relevant)
+ Transition state/boundary
+ Phase separation densities (when relevant)

---------------

Suggested Organization:
 
Each model system should have a README file detailing the basic implementation and bibliography describing the problem.
Each model should have an example sampling algorithm showing how to interface with the model and an example benchmark result.


