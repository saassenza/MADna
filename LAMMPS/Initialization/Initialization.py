#!/usr/bin/env python3

import sys
import os
import numpy as np

import CreateCoordinates
import CreateBonds
import CreateAngles
import CreateDihedrals

if (len(sys.argv) != 4+1):
    print("Usage:") 
    print("# 1 -> input sequence")
    print("# 2 -> working folder")
    print("# 3 -> ionic strength (mM)")
    print("# 4 -> temperature (K)")
    sys.exit(0)

seqstring = sys.argv[1]
foldbase = sys.argv[2]
ionic_strength = float(sys.argv[3])
temperature = float(sys.argv[4])

kDebye = np.sqrt(ionic_strength/temperature)/5.57195
cutDebye = 4.0/kDebye

mass = {
	1: 76.05, # Sugar
	2: 94.97, # Phosphate
    3: 130.09, # ADE
    4: 106.06, # CYT
    5: 146.09, # GUA
    6: 120.07 # THY
}

charge = {
	1: 0, # Sugar
	2: -0.6, # Phosphate
    3: 0, # ADE
    4: 0, # CYT
    5: 0, # GUA
    6: 0 # THY
}

name = {
	1: "Sugar",
	2: "Phosphate",
    3: "ADE",
    4: "CYT", 
    5: "GUA",
    6: "THY"
}

cutLJ = {
	1: 6.40, # Sugar
	2: 4.50, # Phosphate
    3: 5.40, # ADE
    4: 6.40, # CYT
    5: 4.90, # GUA
    6: 7.10 # THY
}

### define filenames
if not os.path.exists(foldbase):
	os.makedirs(foldbase)
foldbase = os.path.realpath(foldbase)
datafile = foldbase+"/stdump.lammpstrj"
chainfile = foldbase+"/chain.dat"
sequencefile = foldbase+"/sequence.dat"
lammpsfile = foldbase+"/lammps.in"
dumpfile = foldbase+"/dump.lammpstrj"


### create coordinates and topology
coordlist, typelist = CreateCoordinates.create_coordinates(seqstring, datafile)
bondtype, i1bond, i2bond = CreateBonds.create_bonds(seqstring)
angletype, i1angle, i2angle, i3angle = CreateAngles.create_angles(seqstring)
dihedraltype, i1dihedral, i2dihedral, i3dihedral, i4dihedral = CreateDihedrals.create_dihedrals(seqstring)

### parameters
lseq = len(seqstring)
natoms = len(coordlist)-1
nbonds = len(bondtype)
nangles = len(angletype)
ndihedrals = len(dihedraltype)

### print chain file
file_writer=open(chainfile,'w')

print("MADna data file for LAMMPS", file=file_writer)
print(str(natoms)+" atoms", file=file_writer)
print(str(nbonds)+" bonds", file=file_writer)
print(str(nangles)+" angles", file=file_writer)
print(str(ndihedrals)+" dihedrals", file=file_writer)

print("6 atom types", file=file_writer)
print("54 bond types", file=file_writer)
print("52 angle types", file=file_writer)
print("96 dihedral types", file=file_writer)

print("-500000 500000 xlo xhi", file=file_writer)
print("-500000 500000 ylo yhi", file=file_writer)
print("-500000 500000 zlo zhi", file=file_writer)

print("", file=file_writer)
print("Masses", file=file_writer)
print("", file=file_writer)
for i in range(1, 7):
    print('{:d} {:.2f} # {:s}'.format(i, mass[i], name[i]), file=file_writer)

print("", file=file_writer)
print("Atoms", file=file_writer)
print("", file=file_writer)
for i in range(1,natoms+1):
    print('{:d} 1 {:d} {:.1f} {:.4f} {:.4f} {:.4f}'.format(i, typelist[i], charge[typelist[i]], coordlist[i,0], coordlist[i,1], coordlist[i,2]), file=file_writer)

print("", file=file_writer)
print("Pair Coeffs", file=file_writer)
print("", file=file_writer)
for i in range(1, 7):
    print('{:d} 1.0 {:.4f} {:.4f} {:.4f}'.format(i, cutLJ[i]/np.power(2, 1./6), cutLJ[i], cutDebye), file=file_writer)


print("", file=file_writer)
print("Bonds", file=file_writer)
print("", file=file_writer)
for i in range(0, nbonds):
	i1 = i1bond[i]
	i2 = i2bond[i]
	tipo = bondtype[i]
	print('{:d} {:d} {:d} {:d}'.format(i+1, tipo, i1, i2), file=file_writer)

print("", file=file_writer)
print("Bond Coeffs", file=file_writer)
print("", file=file_writer)
print("1 54.4248 3.73768  #  Bonds/SP/AA", file=file_writer)
print("2 58.0465 3.73985  #  Bonds/SP/AC", file=file_writer)
print("3 65.2756 3.75037  #  Bonds/SP/AG", file=file_writer)
print("4 67.1491 3.74705  #  Bonds/SP/AT", file=file_writer)
print("5 58.8978 3.75846  #  Bonds/SP/CA", file=file_writer)
print("6 79.3592 3.76515  #  Bonds/SP/CC", file=file_writer)
print("7 60.9381 3.76047  #  Bonds/SP/CG", file=file_writer)
print("8 75.034 3.75193  #  Bonds/SP/CT", file=file_writer)
print("9 49.4744 3.73806  #  Bonds/SP/GA", file=file_writer)
print("10 41.3303 3.70721  #  Bonds/SP/GC", file=file_writer)
print("11 75.1213 3.76077  #  Bonds/SP/GG", file=file_writer)
print("12 54.8989 3.74121  #  Bonds/SP/GT", file=file_writer)
print("13 69.769 3.75839  #  Bonds/SP/TA", file=file_writer)
print("14 73.8681 3.76381  #  Bonds/SP/TC", file=file_writer)
print("15 57.77 3.7583  #  Bonds/SP/TG", file=file_writer)
print("16 79.9486 3.75531  #  Bonds/SP/TT", file=file_writer)
print("17 14.5891 4.08358  #  Bonds/PS/AA", file=file_writer)
print("18 15.0319 4.07002  #  Bonds/PS/AC", file=file_writer)
print("19 11.511 4.09466  #  Bonds/PS/AG", file=file_writer)
print("20 15.8717 4.11528  #  Bonds/PS/AT", file=file_writer)
print("21 14.8208 4.12779  #  Bonds/PS/CA", file=file_writer)
print("22 19.0111 4.1697  #  Bonds/PS/CC", file=file_writer)
print("23 10.4262 4.08759  #  Bonds/PS/CG", file=file_writer)
print("24 18.6914 4.11997  #  Bonds/PS/CT", file=file_writer)
print("25 15.1302 4.10882  #  Bonds/PS/GA", file=file_writer)
print("26 19.2105 4.03813  #  Bonds/PS/GC", file=file_writer)
print("27 16.7518 4.17386  #  Bonds/PS/GG", file=file_writer)
print("28 17.1072 4.11318  #  Bonds/PS/GT", file=file_writer)
print("29 19.2661 4.14302  #  Bonds/PS/TA", file=file_writer)
print("30 19.8115 4.16603  #  Bonds/PS/TC", file=file_writer)
print("31 16.3902 4.13393  #  Bonds/PS/TG", file=file_writer)
print("32 23.6729 4.13829  #  Bonds/PS/TT", file=file_writer)
print("33 39.2297 4.8948  #  Bonds/SB/A", file=file_writer)
print("34 44.5146 4.39384  #  Bonds/SB/C", file=file_writer)
print("35 49.1322 5.01262  #  Bonds/SB/G", file=file_writer)
print("36 45.2514 4.4579  #  Bonds/SB/T", file=file_writer)
print("37 21.484 6.08987  #  Bonds/BB-inter/AT", file=file_writer)
print("38 35.9544 5.70003  #  Bonds/BB-inter/CG", file=file_writer)
print("39 5.86165 3.87173 # Bonds/BB-intra/AA", file=file_writer)
print("40 8.01652 3.7331 # Bonds/BB-intra/AC", file=file_writer)
print("41 3.47674 4.12332 # Bonds/BB-intra/AG", file=file_writer)
print("42 18.5692 3.69623 # Bonds/BB-intra/AT", file=file_writer)
print("43 5.83909 4.23581 # Bonds/BB-intra/CA", file=file_writer)
print("44 2.27823 4.18529 # Bonds/BB-intra/CC", file=file_writer)
print("45 4.46467 4.25942 # Bonds/BB-intra/CG", file=file_writer)
print("46 4.00431 3.94602 # Bonds/BB-intra/CT", file=file_writer)
print("47 6.07068 3.81948 # Bonds/BB-intra/GA", file=file_writer)
print("48 9.52197 3.70822 # Bonds/BB-intra/GC", file=file_writer)
print("49 3.07702 4.14853 # Bonds/BB-intra/GG", file=file_writer)
print("50 13.188 3.69232 # Bonds/BB-intra/GT", file=file_writer)
print("51 8.32623 4.35837 # Bonds/BB-intra/TA", file=file_writer)
print("52 2.95596 4.18839 # Bonds/BB-intra/TC", file=file_writer)
print("53 3.96373 4.48331 # Bonds/BB-intra/TG", file=file_writer)
print("54 4.65664 4.01642 # Bonds/BB-intra/TT", file=file_writer)

print("", file=file_writer)
print("Angles", file=file_writer)
print("", file=file_writer)
for i in range(0, nangles):
    i1 = i1angle[i]
    i2 = i2angle[i]
    i3 = i3angle[i]
    tipo = angletype[i]
    print('{:d} {:d} {:d} {:d} {:d}'.format(i+1, tipo, i1, i2, i3), file=file_writer)

print("", file=file_writer)
print("Angle Coeffs", file=file_writer)
print("", file=file_writer)
print("1 26.1487 94.1182 # Angles/SPS/AA", file=file_writer)
print("2 27.7272 91.8136 # Angles/SPS/AC", file=file_writer)
print("3 23.8434 94.1068 # Angles/SPS/AG", file=file_writer)
print("4 33.0389 92.8322 # Angles/SPS/AT", file=file_writer)
print("5 22.6895 96.2427 # Angles/SPS/CA", file=file_writer)
print("6 32.4694 94.5246 # Angles/SPS/CC", file=file_writer)
print("7 21.8019 95.4622 # Angles/SPS/CG", file=file_writer)
print("8 38.1845 92.8668 # Angles/SPS/CT", file=file_writer)
print("9 22.9078 95.1377 # Angles/SPS/GA", file=file_writer)
print("10 25.4598 92.6603 # Angles/SPS/GC", file=file_writer)
print("11 34.2773 94.2534 # Angles/SPS/GG", file=file_writer)
print("12 28.775 94.4556 # Angles/SPS/GT", file=file_writer)
print("13 26.4028 94.3862 # Angles/SPS/TA", file=file_writer)
print("14 34.1809 93.8147 # Angles/SPS/TC", file=file_writer)
print("15 29.8587 93.1195 # Angles/SPS/TG", file=file_writer)
print("16 46.623 92.6595 # Angles/SPS/TT", file=file_writer)
print("17 51.5777 154.308 # Angles/SBB/AT", file=file_writer)
print("18 53.8562 138.351 # Angles/SBB/CG", file=file_writer)
print("19 52.1523 158.754 # Angles/SBB/GC", file=file_writer)
print("20 49.6802 132.799 # Angles/SBB/TA", file=file_writer)
print("21 32.4685 115.327 # Angles/3PSB5/AA", file=file_writer)
print("22 40.5072 116.095 # Angles/3PSB5/AC", file=file_writer)
print("23 35.9085 115.275 # Angles/3PSB5/AG", file=file_writer)
print("24 43.5787 114.732 # Angles/3PSB5/AT", file=file_writer)
print("25 38.059 119.218 # Angles/3PSB5/CA", file=file_writer)
print("26 51.3764 117.559 # Angles/3PSB5/CC", file=file_writer)
print("27 37.7734 120.498 # Angles/3PSB5/CG", file=file_writer)
print("28 51.275 117.209 # Angles/3PSB5/CT", file=file_writer)
print("29 29.2664 111.774 # Angles/3PSB5/GA", file=file_writer)
print("30 33.1859 111.938 # Angles/3PSB5/GC", file=file_writer)
print("31 52.7018 109.64 # Angles/3PSB5/GG", file=file_writer)
print("32 37.1753 111.275 # Angles/3PSB5/GT", file=file_writer)
print("33 40.8179 119.619 # Angles/3PSB5/TA", file=file_writer)
print("34 48.5882 121.009 # Angles/3PSB5/TC", file=file_writer)
print("35 42.7175 120.089 # Angles/3PSB5/TG", file=file_writer)
print("36 57.9152 119.888 # Angles/3PSB5/TT", file=file_writer)
print("37 6.21218 113.213 # Angles/5PSB3/AA", file=file_writer)
print("38 20.7368 109.287 # Angles/5PSB3/AC", file=file_writer)
print("39 2.43994 119.238 # Angles/5PSB3/AG", file=file_writer)
print("40 20.242 105.789 # Angles/5PSB3/AT", file=file_writer)
print("41 13.4628 110.492 # Angles/5PSB3/CA", file=file_writer)
print("42 19.1926 111.333 # Angles/5PSB3/CC", file=file_writer)
print("43 8.2146 113.341 # Angles/5PSB3/CG", file=file_writer)
print("44 4.25919 106.691 # Angles/5PSB3/CT", file=file_writer)
print("45 1.56112 108.341 # Angles/5PSB3/GA", file=file_writer)
print("46 10.4942 106.681 # Angles/5PSB3/GC", file=file_writer)
print("47 11.4463 121.937 # Angles/5PSB3/GG", file=file_writer)
print("48 26.4827 104.883 # Angles/5PSB3/GT", file=file_writer)
print("49 8.78587 114.263 # Angles/5PSB3/TA", file=file_writer)
print("50 3.25019 107.486 # Angles/5PSB3/TC", file=file_writer)
print("51 15.0148 118.632 # Angles/5PSB3/TG", file=file_writer)
print("52 13.1548 104.955 # Angles/5PSB3/TT", file=file_writer)

print("", file=file_writer)
print("Dihedrals", file=file_writer)
print("", file=file_writer)
for i in range(0, ndihedrals):
    i1 = i1dihedral[i]
    i2 = i2dihedral[i]
    i3 = i3dihedral[i]
    i4 = i4dihedral[i]
    tipo = dihedraltype[i]
    print('{:d} {:d} {:d} {:d} {:d} {:d}'.format(i+1, tipo, i1, i2, i3, i4), file=file_writer)

print("", file=file_writer)
print("Dihedral Coeffs", file=file_writer)
print("", file=file_writer)
print("1 1 59.5504 1 360.429 # Dihedrals/SPSP/AA", file=file_writer)
print("2 1 23.66 1 358.496 # Dihedrals/SPSP/AC", file=file_writer)
print("3 1 47.3402 1 357.837 # Dihedrals/SPSP/AG", file=file_writer)
print("4 1 19.1027 1 359.127 # Dihedrals/SPSP/AT", file=file_writer)
print("5 1 25.628 1 357.359 # Dihedrals/SPSP/CA", file=file_writer)
print("6 1 44.8217 1 356.309 # Dihedrals/SPSP/CC", file=file_writer)
print("7 1 54.62 1 359.207 # Dihedrals/SPSP/CG", file=file_writer)
print("8 1 57.4255 1 359.497 # Dihedrals/SPSP/CT", file=file_writer)
print("9 1 50.0675 1 359.525 # Dihedrals/SPSP/GA", file=file_writer)
print("10 1 52.4038 1 360.979 # Dihedrals/SPSP/GC", file=file_writer)
print("11 1 42.7836 1 356.833 # Dihedrals/SPSP/GG", file=file_writer)
print("12 1 25.742 1 360.701 # Dihedrals/SPSP/GT", file=file_writer)
print("13 1 19.6327 1 357.142 # Dihedrals/SPSP/TA", file=file_writer)
print("14 1 60.4135 1 360.068 # Dihedrals/SPSP/TC", file=file_writer)
print("15 1 32.6503 1 357.335 # Dihedrals/SPSP/TG", file=file_writer)
print("16 1 77.5541 1 362.02 # Dihedrals/SPSP/TT", file=file_writer)
print("17 1 4.7246 1 382.36 # Dihedrals/PSPS/AA", file=file_writer)
print("18 1 4.82176 1 387.392 # Dihedrals/PSPS/AC", file=file_writer)
print("19 1 5.55737 1 386.729 # Dihedrals/PSPS/AG", file=file_writer)
print("20 1 7.50281 1 384.423 # Dihedrals/PSPS/AT", file=file_writer)
print("21 1 5.31944 1 380.003 # Dihedrals/PSPS/CA", file=file_writer)
print("22 1 11.5145 1 383.608 # Dihedrals/PSPS/CC", file=file_writer)
print("23 1 4.11062 1 386.547 # Dihedrals/PSPS/CG", file=file_writer)
print("24 1 10.4625 1 384.682 # Dihedrals/PSPS/CT", file=file_writer)
print("25 1 4.59892 1 379.194 # Dihedrals/PSPS/GA", file=file_writer)
print("26 1 2.74835 1 377.707 # Dihedrals/PSPS/GC", file=file_writer)
print("27 1 10.89 1 383.545 # Dihedrals/PSPS/GG", file=file_writer)
print("28 1 5.60868 1 379.033 # Dihedrals/PSPS/GT", file=file_writer)
print("29 1 8.34713 1 382.107 # Dihedrals/PSPS/TA", file=file_writer)
print("30 1 10.7596 1 379.317 # Dihedrals/PSPS/TC", file=file_writer)
print("31 1 7.74714 1 386.504 # Dihedrals/PSPS/TG", file=file_writer)
print("32 1 15.6575 1 380.518 # Dihedrals/PSPS/TT", file=file_writer)
print("33 1 13.8177 1 217.662 # Dihedrals/SPSB53/AA", file=file_writer)
print("34 1 19.0939 1 218.714 # Dihedrals/SPSB53/AC", file=file_writer)
print("35 1 10.1199 1 214.211 # Dihedrals/SPSB53/AG", file=file_writer)
print("36 1 27.1357 1 222.578 # Dihedrals/SPSB53/AT", file=file_writer)
print("37 1 10.957 1 215.366 # Dihedrals/SPSB53/CA", file=file_writer)
print("38 1 21.61 1 219.639 # Dihedrals/SPSB53/CC", file=file_writer)
print("39 1 7.56341 1 210.442 # Dihedrals/SPSB53/CG", file=file_writer)
print("40 1 27.5272 1 223.282 # Dihedrals/SPSB53/CT", file=file_writer)
print("41 1 16.5934 1 220.19 # Dihedrals/SPSB53/GA", file=file_writer)
print("42 1 21.4059 1 222.059 # Dihedrals/SPSB53/GC", file=file_writer)
print("43 1 13.6335 1 218.983 # Dihedrals/SPSB53/GG", file=file_writer)
print("44 1 28.0474 1 224.802 # Dihedrals/SPSB53/GT", file=file_writer)
print("45 1 14.2121 1 216.144 # Dihedrals/SPSB53/TA", file=file_writer)
print("46 1 23.5475 1 221.985 # Dihedrals/SPSB53/TC", file=file_writer)
print("47 1 9.02256 1 213.077 # Dihedrals/SPSB53/TG", file=file_writer)
print("48 1 34.024 1 224.39 # Dihedrals/SPSB53/TT", file=file_writer)
print("49 1 10.3949 1 163.974 # Dihedrals/SPSB35/AA", file=file_writer)
print("50 1 10.9112 1 166.223 # Dihedrals/SPSB35/AC", file=file_writer)
print("51 1 8.93809 1 165.093 # Dihedrals/SPSB35/AG", file=file_writer)
print("52 1 14.8604 1 164.452 # Dihedrals/SPSB35/AT", file=file_writer)
print("53 1 6.91255 1 155.621 # Dihedrals/SPSB35/CA", file=file_writer)
print("54 1 13.9231 1 157.368 # Dihedrals/SPSB35/CC", file=file_writer)
print("55 1 5.55897 1 159.631 # Dihedrals/SPSB35/CG", file=file_writer)
print("56 1 13.4308 1 158.951 # Dihedrals/SPSB35/CT", file=file_writer)
print("57 1 11.5279 1 165.397 # Dihedrals/SPSB35/GA", file=file_writer)
print("58 1 11.0388 1 165.796 # Dihedrals/SPSB35/GC", file=file_writer)
print("59 1 15.0442 1 166.73 # Dihedrals/SPSB35/GG", file=file_writer)
print("60 1 14.5422 1 165.124 # Dihedrals/SPSB35/GT", file=file_writer)
print("61 1 8.69143 1 152.127 # Dihedrals/SPSB35/TA", file=file_writer)
print("62 1 12.006 1 149.365 # Dihedrals/SPSB35/TC", file=file_writer)
print("63 1 6.70681 1 154.621 # Dihedrals/SPSB35/TG", file=file_writer)
print("64 1 15.7239 1 151.877 # Dihedrals/SPSB35/TT", file=file_writer)
print("65 1 10.4952 1 380.936 # Dihedrals/PSBB53/AA", file=file_writer)
print("66 1 13.6367 1 379.792 # Dihedrals/PSBB53/AC", file=file_writer)
print("67 1 6.69763 1 382.429 # Dihedrals/PSBB53/AG", file=file_writer)
print("68 1 23.8987 1 381.437 # Dihedrals/PSBB53/AT", file=file_writer)
print("69 1 7.62427 1 383.127 # Dihedrals/PSBB53/CA", file=file_writer)
print("70 1 17.9107 1 368.613 # Dihedrals/PSBB53/CC", file=file_writer)
print("71 1 6.38959 1 398.225 # Dihedrals/PSBB53/CG", file=file_writer)
print("72 1 18.0406 1 381.424 # Dihedrals/PSBB53/CT", file=file_writer)
print("73 1 10.6723 1 381.206 # Dihedrals/PSBB53/GA", file=file_writer)
print("74 1 15.3054 1 383.863 # Dihedrals/PSBB53/GC", file=file_writer)
print("75 1 7.2435 1 370.348 # Dihedrals/PSBB53/GG", file=file_writer)
print("76 1 24.1227 1 382.28 # Dihedrals/PSBB53/GT", file=file_writer)
print("77 1 8.40469 1 387.757 # Dihedrals/PSBB53/TA", file=file_writer)
print("78 1 12.935 1 379.824 # Dihedrals/PSBB53/TC", file=file_writer)
print("79 1 6.90351 1 395.889 # Dihedrals/PSBB53/TG", file=file_writer)
print("80 1 21.2987 1 381.428 # Dihedrals/PSBB53/TT", file=file_writer)
print("81 1 10.116 1 235.218 # Dihedrals/PSBB35/AA", file=file_writer)
print("82 1 7.4467 1 239.788 # Dihedrals/PSBB35/AC", file=file_writer)
print("83 1 8.55808 1 239.834 # Dihedrals/PSBB35/AG", file=file_writer)
print("84 1 8.17615 1 243.855 # Dihedrals/PSBB35/AT", file=file_writer)
print("85 1 8.98152 1 242.991 # Dihedrals/PSBB35/CA", file=file_writer)
print("86 1 13.9191 1 233.617 # Dihedrals/PSBB35/CC", file=file_writer)
print("87 1 8.85158 1 246.915 # Dihedrals/PSBB35/CG", file=file_writer)
print("88 1 10.0393 1 240.077 # Dihedrals/PSBB35/CT", file=file_writer)
print("89 1 7.92129 1 235.213 # Dihedrals/PSBB35/GA", file=file_writer)
print("90 1 6.7124 1 245.904 # Dihedrals/PSBB35/GC", file=file_writer)
print("91 1 7.91363 1 229.423 # Dihedrals/PSBB35/GG", file=file_writer)
print("92 1 7.55626 1 246.524 # Dihedrals/PSBB35/GT", file=file_writer)
print("93 1 12.1136 1 247.238 # Dihedrals/PSBB35/TA", file=file_writer)
print("94 1 12.9901 1 247.848 # Dihedrals/PSBB35/TC", file=file_writer)
print("95 1 10.7031 1 249.333 # Dihedrals/PSBB35/TG", file=file_writer)
print("96 1 14.6266 1 244.19 # Dihedrals/PSBB35/TT", file=file_writer)

file_writer.close()

####### write sequence
file_writer=open(sequencefile,'w')
print(seqstring, file=file_writer)
file_writer.close()


###### write sample input file for lammps
file_writer=open(lammpsfile,'w')
print("units real", file=file_writer)
print("boundary p p p", file=file_writer)
print("dimension 3", file=file_writer)
print("atom_style full", file=file_writer)
print("atom_modify sort 0 0", file=file_writer)
print("bond_style harmonic", file=file_writer)
print("angle_style harmonic", file=file_writer)
print("dihedral_style fourier", file=file_writer)
print('pair_style lj/cut/coul/debye {:.4f} 20 {:.4f}'.format(kDebye, cutDebye), file=file_writer)
print("pair_modify shift yes mix arithmetic ", file=file_writer)
print("special_bonds lj 0.0 0.0 0.0 coul 1.0 1.0 1.0", file=file_writer)
print("dielectric 78", file=file_writer)
print("read_data "+chainfile, file=file_writer)
print("", file=file_writer)
print("neighbor 2 nsq", file=file_writer)
print("neigh_modify delay 0 every 1 check yes page 500000 one 50000", file=file_writer)
print("", file=file_writer)
print("########## Additional features (e.g. pulling force) should be added below this line", file=file_writer)
print("########## ", file=file_writer)
print("########## Additional features (e.g. pulling force) should be added above this line", file=file_writer)
print("", file=file_writer)
print("fix rkick all langevin {:.4f} {:.4f} 20000 {:d}".format(temperature, temperature, np.random.randint(1, high=32768)), file=file_writer)
print("fix evol all nve", file=file_writer)
print("", file=file_writer)
print("thermo 500", file=file_writer)
print("variable tns equal time/1e6", file=file_writer)
print("variable cpuh equal cpuremain/3600", file=file_writer)
print("thermo_style custom v_tns temp evdwl ecoul ebond eangle edihed pe v_cpuh", file=file_writer)
print("velocity all create {:.4f} {:d}".format(temperature, np.random.randint(1, high=32768)), file=file_writer)
print("", file=file_writer)
print("timestep 20", file=file_writer)
print("", file=file_writer)
print("dump coord all custom 500 {:s} id type xu yu zu".format(dumpfile), file=file_writer)
print("dump_modify coord sort id", file=file_writer)
print("run 500000", file=file_writer)
file_writer.close()
