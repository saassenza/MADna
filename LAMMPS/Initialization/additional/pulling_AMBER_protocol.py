#!/usr/bin/env python3

import sys

if (len(sys.argv) != 2+1):
    print("Usage:")
    print("# 1 -> input sequence")
    print("# 2 -> pulling force (pN)")
    sys.exit(0)

seqstring = sys.argv[1]
pulling_force_pN = float(sys.argv[2])

pulling_force_lammps = pulling_force_pN*0.5/69.5118
lseq = len(seqstring)
natoms = 6*lseq - 2
idA1 = 4
idA2 = natoms - idA1
idB1 = int(natoms/2) - 4
idB2 = natoms - idB1

print ("compute xu all property/atom xu")
print ("compute yu all property/atom yu")
print ("compute zu all property/atom zu")
print ("variable xcmA equal 0.5*(c_xu[{:d}]+c_xu[{:d}])".format(idA1, idA2))
print ("variable ycmA equal 0.5*(c_yu[{:d}]+c_yu[{:d}])".format(idA1, idA2))
print ("variable zcmA equal 0.5*(c_zu[{:d}]+c_zu[{:d}])".format(idA1, idA2))
print ("variable xcmB equal 0.5*(c_xu[{:d}]+c_xu[{:d}])".format(idB1, idB2))
print ("variable ycmB equal 0.5*(c_yu[{:d}]+c_yu[{:d}])".format(idB1, idB2))
print ("variable zcmB equal 0.5*(c_zu[{:d}]+c_zu[{:d}])".format(idB1, idB2))
print ("variable dr equal sqrt((v_xcmB-v_xcmA)*(v_xcmB-v_xcmA)+(v_ycmB-v_ycmA)*(v_ycmB-v_ycmA)+(v_zcmB-v_zcmA)*(v_zcmB-v_zcmA))")
print ("variable fxA equal -(v_xcmB-v_xcmA)*{:.4f}/v_dr".format(pulling_force_lammps))
print ("variable fyA equal -(v_ycmB-v_ycmA)*{:.4f}/v_dr".format(pulling_force_lammps))
print ("variable fzA equal -(v_zcmB-v_zcmA)*{:.4f}/v_dr".format(pulling_force_lammps))
print ("variable fxB equal -v_fxA")
print ("variable fyB equal -v_fyA")
print ("variable fzB equal -v_fzA")
print ("group fA id {:d} {:d}".format(idA1, idA2))
print ("group fB id {:d} {:d}".format(idA1, idA2))
print ("fix fA fA addforce v_fxA v_fyA v_fzA")
print ("fix fB fB addforce v_fxB v_fyB v_fzB")
