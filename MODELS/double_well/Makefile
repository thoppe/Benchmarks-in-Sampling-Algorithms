# Run the simulations

sim_umbrella:
	python alg_EnergyBarrier_umbrella.py simulation_setups/umbrella_EnergyBarrier.json

sim_equilibrium:
	python alg_EnergyBarrier_equilibrium.py simulation_setups/equilibrium_EnergyBarrier.json

parallel_equilibrium:
	parallel python alg_EnergyBarrier_equilibrium.py simulation_setups/equilibrium_EnergyBarrier.json --replica_n ::: 0 1 2 3 5

sim_replicaEx:
	python alg_EnergyBarrier_replicaEx.py simulation_setups/replicaEx_EnergyBarrier.json

sim_example:
	python alg_EnergyBarrier_equilibrium.py simulation_setups/example_EnergyBarrier.json

plot_results:
	python src/plots_results.py

clean:
	rm -rvf results/* trajectory/*.txt trajectory/*.npy

