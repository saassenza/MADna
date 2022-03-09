# MADna-LAMMPS
LAMMPS implementation of MADna, a coarse-grained model for sequence-dependent elasticity and conformation of DNA

There are two folders, for which more details are provided below:
- **Initialization**: This folder contains all the scripts needed to create the topology and initial configuration of a double-stranded DNA molecule from its sequence
- **Analysis**: This folder contains analysis scripts to determine geometrical parameters of DNA (grooves, hrise, htwist, diameter, crookedness, extension) and the tangent-tangent correlation function needed to compute the persistence length

## Initialization
This folder contains scripts for generation of MADna topology and initial configuration of a double-stranded DNA molecule in LAMMPS format. To execute the scripts, **python3 is needed with libraries sys, os, numpy, copy**.  
From the user's perspective, the 
