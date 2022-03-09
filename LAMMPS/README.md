# MADna-LAMMPS
LAMMPS implementation of MADna, a coarse-grained model for sequence-dependent elasticity and conformation of DNA. In order to run MADna within LAMMPS, the latter has to be built with the packages EXTRA-MOLECULE and EXTRA-PAIR. Morover, for application of torque according to the scripts provided in the folder "additional" (see below), also the package EXTRA-FIX is needed. 

There are two folders, for which more details are provided below:
- **Initialization**: This folder contains all the scripts needed to create the topology and initial configuration of a double-stranded DNA molecule from its sequence
- **Analysis**: This folder contains analysis scripts to determine geometrical parameters of DNA (grooves, hrise, htwist, diameter, crookedness, extension) and the tangent-tangent correlation function needed to compute the persistence length

## Initialization
This folder contains scripts for generation of MADna topology and initial configuration of a double-stranded DNA molecule in LAMMPS format. To execute the scripts, **python3 is needed with libraries sys, os, numpy, copy**.  
From the user's perspective, the only relevant script is **Initialization.py**, which imports the other ones as libraries (further useful scripts are found in the subfolder "additional", as commented below). Usage of the script is pretty straightforward:
```
Initialization/Initialization.py sequence folder ionic_strength temperature
```
Hence, the script needs four inputs:
- sequence: the sequence of the leading strand in the 5'-3' direction
- folder: the folder where the various files will be stored (if the folder does not exist, it is created automatically)
- ionic strength: the ionic strength of the solution in mM
- temperature: the temperature of the system in K

For instance, the following command
```
Initialization/Initialization.py CAAGATGC mySim/initialization 150 300
```
creates the files needed to simulate a DNA molecule with sequence 5'-CAAGATGC-3' (corresponding to 5'-GCATCTTG-3' on the other strand) embedded in a solution containing 150 mM of salt and at temperature 300 K, and stores the file in the folder mySim/initialization

The script generates four files:
- **sequence.dat** contains the chosen sequence as a simple text
- **stdump.lammpstrj** contains the initial coordinates of the beads within the molecule, written in LAMMPS format
- **chain.dat** contains the topology of the molecule in LAMMPS format (units correspond to the "real" format in LAMMPS)
- **lammps.in** is a minimal script for simulation of the system in LAMMPS, where the dynamics is run for 10 ns with a dump every 10 ps. Note that a group of commented lines indicates the point of the script where additional features can be added (see below)
