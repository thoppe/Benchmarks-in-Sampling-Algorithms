import pylab as plt
import seaborn as sns
import numpy as np
import glob

results_list = ["results/example_EnergyBarrier*",
                "results/replicaEx_EnergyBarrier**",]
labels_list  = ["Equilibrium sampling", "Replica Exchange"]

color_list   = sns.color_palette("muted", len(results_list))

for c,label,phrase in zip(color_list, 
                          labels_list,
                          results_list):

    file_list = glob.glob(phrase)

    raw_data = np.array([np.loadtxt(f_result) for f_result in file_list])
    err  = raw_data[:,:, 1]
    time = raw_data[0,:, 0]

    # Manual correction
    if "replica" in phrase: 
        time *= raw_data.shape[0]

    sns.tsplot(err,time, color=c, condition=label)

plt.ylim(0,0.5)
plt.legend(loc="best")
plt.show()

               



