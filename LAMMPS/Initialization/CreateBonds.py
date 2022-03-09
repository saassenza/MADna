#!/usr/bin/env python3

WC = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A'
}

typeBond = {
    'SP/AA': 1,
    'SP/AC': 2,
    'SP/AG': 3,
    'SP/AT': 4,
    'SP/CA': 5,
    'SP/CC': 6,
    'SP/CG': 7,
    'SP/CT': 8,
    'SP/GA': 9,
    'SP/GC': 10,
    'SP/GG': 11,
    'SP/GT': 12,
    'SP/TA': 13,
    'SP/TC': 14,
    'SP/TG': 15,
    'SP/TT': 16,
    'PS/AA': 17,
    'PS/AC': 18,
    'PS/AG': 19,
    'PS/AT': 20,
    'PS/CA': 21,
    'PS/CC': 22,
    'PS/CG': 23,
    'PS/CT': 24,
    'PS/GA': 25,
    'PS/GC': 26,
    'PS/GG': 27,
    'PS/GT': 28,
    'PS/TA': 29,
    'PS/TC': 30,
    'PS/TG': 31,
    'PS/TT': 32,
    'SB/A': 33,
    'SB/C': 34,
    'SB/G': 35,
    'SB/T': 36,
    'BB-inter/AT': 37,
    'BB-inter/CG': 38,
    'BB-inter/GC': 38,
    'BB-inter/TA': 37,
    'BB-intra/AA': 39,
    'BB-intra/AC': 40,
    'BB-intra/AG': 41,
    'BB-intra/AT': 42,
    'BB-intra/CA': 43,
    'BB-intra/CC': 44,
    'BB-intra/CG': 45,
    'BB-intra/CT': 46,
    'BB-intra/GA': 47,
    'BB-intra/GC': 48,
    'BB-intra/GG': 49,
    'BB-intra/GT': 50,
    'BB-intra/TA': 51,
    'BB-intra/TC': 52,
    'BB-intra/TG': 53,
    'BB-intra/TT': 54
}

def create_bonds(seqstring):
    lseq = len(seqstring)
    bondtype = []
    i1list = []
    i2list = []

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

    #### SP
    for i in range(1,lseq):
        i1 = 3*i - 2
        i2 = i1 + 2
        tipo = typeBond['SP/'+seqstring[i-1]+seqstring[i]]
        i1list.append(i1)
        i2list.append(i2)
        bondtype.append(tipo)
    for i in range(1,lseq):
        i1 = lastB1 + 3*i - 2
        i2 = i1 + 2
        tipo = typeBond['SP/'+WC[seqstring[lseq-1-(i-1)]]+WC[seqstring[lseq-1-i]]]
        i1list.append(i1)
        i2list.append(i2)
        bondtype.append(tipo)
    
    #### PS
    for i in range(2,lseq+1):
        i1 = 3*i - 3
        i2 = i1 + 1
        tipo = typeBond['PS/'+seqstring[i-2]+seqstring[i-1]]
        i1list.append(i1)
        i2list.append(i2)
        bondtype.append(tipo)
    for i in range(2,lseq+1):
        i1 = lastB1 + 3*i - 3
        i2 = i1 + 1
        tipo = typeBond['PS/'+WC[seqstring[lseq-1-(i-2)]]+WC[seqstring[lseq-1-(i-1)]]]
        i1list.append(i1)
        i2list.append(i2)
        bondtype.append(tipo)
    
    #### SB
    for i in range(1,lseq+1):
        i1 = 3*i - 2
        i2 = i1 + 1
        tipo = typeBond['SB/'+seqstring[i-1]]
        i1list.append(i1)
        i2list.append(i2)
        bondtype.append(tipo)
    for i in range(1,lseq+1):
        i1 = lastB1 + 3*i - 2
        i2 = i1 + 1
        tipo = typeBond['SB/'+WC[seqstring[lseq-1-(i-1)]]]
        i1list.append(i1)
        i2list.append(i2)
        bondtype.append(tipo)
    
    #### BB-inter
    for i in range(1,lseq+1):
        i1 = 3*i - 1
        i2 = lastB2 + 2 - i1
        tipo = typeBond['BB-inter/'+seqstring[i-1]+WC[seqstring[i-1]]]
        i1list.append(i1)
        i2list.append(i2)
        bondtype.append(tipo)
    
    #### BB-intra
    for i in range(1,lseq):
        i1 = 3*i - 1
        i2 = i1 + 3
        tipo = typeBond['BB-intra/'+seqstring[i-1]+seqstring[i]]
        i1list.append(i1)
        i2list.append(i2)
        bondtype.append(tipo)
    for i in range(1,lseq):
        i1 = lastB1 + 3*i - 1
        i2 = i1 + 3
        tipo = typeBond['BB-intra/'+WC[seqstring[lseq-1-(i-1)]]+WC[seqstring[lseq-1-i]]]
        i1list.append(i1)
        i2list.append(i2)
        bondtype.append(tipo)
    
    return bondtype, i1list, i2list
