# Alanine dipeptide
+ **Model**: 
Direct MD of [Alaine Dipeptide](http://en.wikipedia.org/wiki/Dipeptide): 
Gromacs + Charmm27
+ **Author**: Bradley M Dickson
+ **NOTE**: Preliminary steps for alaine dipeptide.

This is a 500ns trajectory of alanine dipeptide in tip3p water.
The simulation was run under GROMACS 4.5.5.
(500ns.zip)[500ns.zip] contains a file `traj.pdb` which contains the coordinates from the simulation. 
You can make any analysis of this data that you'd like, and we may follow this up by posting a phi-psi trace of the trajectory and the corresponding free energy estimate.
   
We also provide `isob.gro`, `isob.cpt`, `topol.top` and `md.mdp` which can be used to generate inputs for a new GROMACS run. 
If you have installed GROMACS, use:

    PATHTOYOUR/GROMACS/bin/grompp -f md.mdp -c isob.gro -t isob.cpt -p topol.top -o Run1.tpr

Now, skipping all the details you'll have to sort for your local machines, you can launch a new simulation with: 

    PATHTOYOUR/GROMACS/bin/mdrun -deffnm Run1

Nothing special was done for "md.mdp" or the NPT equilibration. 
The "isob" files are from the last frame of 3ns NPT. 

**Please read (topol.top)[topol.top] for force field details.**
