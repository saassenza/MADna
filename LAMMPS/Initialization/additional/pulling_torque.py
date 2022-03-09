#!/usr/bin/env python3

import sys

if (len(sys.argv) != 3+1):
    print("Usage:")
    print("# 1 -> input sequence")
    print("# 2 -> pulling force (pN)")
    print("# 3 -> applied torque (pN nm)")
    sys.exit(0)

seqstring = sys.argv[1]
pulling_force_pN = float(sys.argv[2])
torque_pNnm = float(sys.argv[3])

pulling_force_lammps = pulling_force_pN*0.5/69.5118
torque_lammps = torque_pNnm*0.14386

lseq = len(seqstring)
natoms = 6*lseq - 2
idA1 = 4
idA2 = natoms - idA1
idB1 = int(natoms/2) - 4
idB2 = natoms - idB1

print ("compute xu all property/atom xu")
print ("compute yu all property/atom yu")
print ("variable dx equal -0.5*(c_xu[{:d}]+c_xu[{:d}])".format(idB1, idB2))
print ("variable dy equal -0.5*(c_yu[{:d}]+c_yu[{:d}])".format(idB1, idB2))
print ("thermo_style custom v_dx v_dy")
print ("run 0")
print ("displace_atoms all move v_dx v_dy 0 units box")
print ("thermo_style one")
print ("group blockA1 id <= {:d}".format(idA1+1))
print ("group blockA2 id >= {:d}".format(idA2))
print ("group blockA union blockA1 blockA2")
print ("group fB id {:d} {:d}".format(idB1, idB2))
print ("group torque id <> {:d} {:d}".format(idB1, idB2))
print ("fix tetherA blockA spring/self 100 xyz")
print ("variable fxB equal -100*(c_xu[{:d}]+c_xu[{:d}])".format(idB1, idB2))
print ("variable fyB equal -100*(c_yu[{:d}]+c_yu[{:d}])".format(idB1, idB2))
print ("fix fB fB addforce v_fxB v_fyB {:.4f}".format(pulling_force_lammps))
print ("fix torque torque addtorque 0 0 {:.4f}".format(torque_lammps))
