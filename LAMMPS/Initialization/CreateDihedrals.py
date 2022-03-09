#!/usr/bin/env python3

WC = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A'
}

typeDihedral = {
    'SPSP/AA': 1,
    'SPSP/AC': 2,
    'SPSP/AG': 3,
    'SPSP/AT': 4,
    'SPSP/CA': 5,
    'SPSP/CC': 6,
    'SPSP/CG': 7,
    'SPSP/CT': 8,
    'SPSP/GA': 9,
    'SPSP/GC': 10,
    'SPSP/GG': 11,
    'SPSP/GT': 12,
    'SPSP/TA': 13,
    'SPSP/TC': 14,
    'SPSP/TG': 15,
    'SPSP/TT': 16,
    'PSPS/AA': 17,
    'PSPS/AC': 18,
    'PSPS/AG': 19,
    'PSPS/AT': 20,
    'PSPS/CA': 21,
    'PSPS/CC': 22,
    'PSPS/CG': 23,
    'PSPS/CT': 24,
    'PSPS/GA': 25,
    'PSPS/GC': 26,
    'PSPS/GG': 27,
    'PSPS/GT': 28,
    'PSPS/TA': 29,
    'PSPS/TC': 30,
    'PSPS/TG': 31,
    'PSPS/TT': 32,
    'SPSB53/AA': 33,
    'SPSB53/AC': 34,
    'SPSB53/AG': 35,
    'SPSB53/AT': 36,
    'SPSB53/CA': 37,
    'SPSB53/CC': 38,
    'SPSB53/CG': 39,
    'SPSB53/CT': 40,
    'SPSB53/GA': 41,
    'SPSB53/GC': 42,
    'SPSB53/GG': 43,
    'SPSB53/GT': 44,
    'SPSB53/TA': 45,
    'SPSB53/TC': 46,
    'SPSB53/TG': 47,
    'SPSB53/TT': 48,
    'SPSB35/AA': 49,
    'SPSB35/AC': 50,
    'SPSB35/AG': 51,
    'SPSB35/AT': 52,
    'SPSB35/CA': 53,
    'SPSB35/CC': 54,
    'SPSB35/CG': 55,
    'SPSB35/CT': 56,
    'SPSB35/GA': 57,
    'SPSB35/GC': 58,
    'SPSB35/GG': 59,
    'SPSB35/GT': 60,
    'SPSB35/TA': 61,
    'SPSB35/TC': 62,
    'SPSB35/TG': 63,
    'SPSB35/TT': 64,
    'PSBB53/AA': 65,
    'PSBB53/AC': 66,
    'PSBB53/AG': 67,
    'PSBB53/AT': 68,
    'PSBB53/CA': 69,
    'PSBB53/CC': 70,
    'PSBB53/CG': 71,
    'PSBB53/CT': 72,
    'PSBB53/GA': 73,
    'PSBB53/GC': 74,
    'PSBB53/GG': 75,
    'PSBB53/GT': 76,
    'PSBB53/TA': 77,
    'PSBB53/TC': 78,
    'PSBB53/TG': 79,
    'PSBB53/TT': 80,
    'PSBB35/AA': 81,
    'PSBB35/AC': 82,
    'PSBB35/AG': 83,
    'PSBB35/AT': 84,
    'PSBB35/CA': 85,
    'PSBB35/CC': 86,
    'PSBB35/CG': 87,
    'PSBB35/CT': 88,
    'PSBB35/GA': 89,
    'PSBB35/GC': 90,
    'PSBB35/GG': 91,
    'PSBB35/GT': 92,
    'PSBB35/TA': 93,
    'PSBB35/TC': 94,
    'PSBB35/TG': 95,
    'PSBB35/TT': 96
}

def create_dihedrals(seqstring):
   lseq = len(seqstring)
   dihedraltype = []
   i1list = [] 
   i2list = [] 
   i3list = [] 
   i4list = [] 
   
   firstS1 = 1
   lastS1 = firstS1 + (lseq-1)*3
   firstP1 = 3
   lastP1 = firstP1 + (lseq-2)*3
   firstB1 = 2
   lastB1 = firstB1 + (lseq-1)*3
   
   firstS2 = lastB1 + 1
   lastS2 = firstS2 + (lseq-1)*3
   firstP2 = lastB1 + 3
   lastP2 = firstP2 + (lseq-2)*3
   firstB2 = lastB1 + 2
   lastB2 = firstB2 + (lseq-1)*3
   
   idihedral = 1
   
   #### SPSP
   for i in range(1,lseq-1):
       i1 = 3*i - 2
       i2 = i1 + 2
       i3 = i1 + 3
       i4 = i1 + 5
       tipo = typeDihedral['SPSP/'+seqstring[i-1]+seqstring[i]]
       i1list.append(i1)
       i2list.append(i2)
       i3list.append(i3)
       i4list.append(i4)
       dihedraltype.append(tipo)
   for i in range(1,lseq-1):
       i1 = lastB1 + 3*i - 2
       i2 = i1 + 2
       i3 = i1 + 3
       i4 = i1 + 5
       tipo = typeDihedral['SPSP/'+WC[seqstring[lseq-1-(i-1)]]+WC[seqstring[lseq-1-i]]]
       i1list.append(i1)
       i2list.append(i2)
       i3list.append(i3)
       i4list.append(i4)
       dihedraltype.append(tipo)
   
   #### PSPS
   for i in range(2,lseq):
       i1 = 3*i - 3
       i2 = i1 + 1
       i3 = i1 + 3
       i4 = i1 + 4
       tipo = typeDihedral['PSPS/'+seqstring[i-1]+seqstring[i]]
       i1list.append(i1)
       i2list.append(i2)
       i3list.append(i3)
       i4list.append(i4)
       dihedraltype.append(tipo)
   for i in range(2,lseq):
       i1 = lastB1 + 3*i - 3
       i2 = i1 + 1
       i3 = i1 + 3
       i4 = i1 + 4
       tipo = typeDihedral['PSPS/'+WC[seqstring[lseq-1-(i-1)]]+WC[seqstring[lseq-1-i]]]
       i1list.append(i1)
       i2list.append(i2)
       i3list.append(i3)
       i4list.append(i4)
       dihedraltype.append(tipo)
   
   #### SPSB53
   for i in range(1,lseq):
       i1 = 3*i - 2
       i2 = i1 + 2
       i3 = i1 + 3
       i4 = i1 + 4
       tipo = typeDihedral['SPSB53/'+seqstring[i-1]+seqstring[i]]
       i1list.append(i1)
       i2list.append(i2)
       i3list.append(i3)
       i4list.append(i4)
       dihedraltype.append(tipo)
   for i in range(1,lseq):
       i1 = lastB1 + 3*i - 2
       i2 = i1 + 2
       i3 = i1 + 3
       i4 = i1 + 4
       tipo = typeDihedral['SPSB53/'+WC[seqstring[lseq-1-(i-1)]]+WC[seqstring[lseq-1-i]]]
       i1list.append(i1)
       i2list.append(i2)
       i3list.append(i3)
       i4list.append(i4)
       dihedraltype.append(tipo)
   
   #### SPSB35
   for i in range(1,lseq):
       i1 = 3*i - 1
       i2 = i1 - 1
       i3 = i1 + 1
       i4 = i1 + 2
       tipo = typeDihedral['SPSB35/'+seqstring[i-1]+seqstring[i]]
       i1list.append(i1)
       i2list.append(i2)
       i3list.append(i3)
       i4list.append(i4)
       dihedraltype.append(tipo)
   for i in range(1,lseq):
       i1 = lastB1 + 3*i - 1
       i2 = i1 - 1
       i3 = i1 + 1
       i4 = i1 + 2
       tipo = typeDihedral['SPSB35/'+WC[seqstring[lseq-1-(i-1)]]+WC[seqstring[lseq-1-i]]]
       i1list.append(i1)
       i2list.append(i2)
       i3list.append(i3)
       i4list.append(i4)
       dihedraltype.append(tipo)
   
   #### PSBB53
   for i in range(1,lseq):
       i1 = 3*i
       i2 = i1 + 1
       i3 = i1 + 2
       i4 = lastB2 + 2 - i3
       tipo = typeDihedral['PSBB53/'+seqstring[i-1]+seqstring[i]]
       i1list.append(i1)
       i2list.append(i2)
       i3list.append(i3)
       i4list.append(i4)
       dihedraltype.append(tipo)
   for i in range(1,lseq):
       i1 = lastB1 + 3*i
       i2 = i1 + 1
       i3 = i1 + 2
       i4 = lastB2 + 2 - i3
       tipo = typeDihedral['PSBB53/'+WC[seqstring[lseq-1-(i-1)]]+WC[seqstring[lseq-1-i]]]
       i1list.append(i1)
       i2list.append(i2)
       i3list.append(i3)
       i4list.append(i4)
       dihedraltype.append(tipo)
   
   #### PSBB35
   for i in range(1,lseq):
       i1 = 3*i
       i2 = i1 - 2
       i3 = i1 - 1
       i4 = lastB2 + 2 - i3
       tipo = typeDihedral['PSBB35/'+seqstring[i-1]+seqstring[i]]
       i1list.append(i1)
       i2list.append(i2)
       i3list.append(i3)
       i4list.append(i4)
       dihedraltype.append(tipo)
   for i in range(1,lseq):
       i1 = lastB1 + 3*i
       i2 = i1 - 2
       i3 = i1 - 1
       i4 = lastB2 + 2 - i3
       tipo = typeDihedral['PSBB35/'+WC[seqstring[lseq-1-(i-1)]]+WC[seqstring[lseq-1-i]]]
       i1list.append(i1)
       i2list.append(i2)
       i3list.append(i3)
       i4list.append(i4)
       dihedraltype.append(tipo)
   
   return dihedraltype, i1list, i2list, i3list, i4list
