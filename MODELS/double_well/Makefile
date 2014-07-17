# Run the simulations

sim_equilibrium:
	python alg_EnergyBarrier_equilibrium.py simulation_setups/example_EnergyBarrier.json

parallel_equilibrium:
	parallel python alg_EnergyBarrier_equilibrium.py simulation_setups/example_EnergyBarrier.json --replica_n ::: 0 1 2 3

sim_replicaEx:
	python alg_EnergyBarrier_replicaEx.py simulation_setups/replicaEx_EnergyBarrier.json

plot_results:
	python src/plots_results.py
