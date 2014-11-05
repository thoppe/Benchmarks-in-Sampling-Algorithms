# Run the simulations

sim_example:
	python method_equilibrium.py simulation_setups/equilibrium_example.json

sim_equilibrium:
	python method_equilibrium.py simulation_setups/equilibrium_sim.json

parallel_equilibrium:
	parallel python method_equilibrium.py simulation_setups/equilibrium_sim.json --replica_n ::: 0 1 2 3 5

sim_replicaEx:
	python method_replicaEx.py simulation_setups/replicaEx_sim.json

sim_umbrella:
	python method_umbrella.py simulation_setups/umbrella_sim.json

############################################################################

plot_results:
	python src/plots_results.py

############################################################################

clean:
	rm -rvf results/* trajectory/*.txt trajectory/*.npy

