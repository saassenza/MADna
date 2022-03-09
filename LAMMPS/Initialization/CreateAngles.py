#!/usr/bin/env python3

WC = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A'
}

typeAngle = {
    'SPS/AA': 1,
    'SPS/AC': 2,
    'SPS/AG': 3,
    'SPS/AT': 4,
    'SPS/CA': 5,
    'SPS/CC': 6,
    'SPS/CG': 7,
    'SPS/CT': 8,
    'SPS/GA': 9,
    'SPS/GC': 10,
    'SPS/GG': 11,
    'SPS/GT': 12,
    'SPS/TA': 13,
    'SPS/TC': 14,
    'SPS/TG': 15,
    'SPS/TT': 16,
    'SBB/AT': 17,
    'SBB/CG': 18,
    'SBB/GC': 19,
    'SBB/TA': 20,
    '3PSB5/AA': 21,
    '3PSB5/AC': 22,
    '3PSB5/AG': 23,
    '3PSB5/AT': 24,
    '3PSB5/CA': 25,
    '3PSB5/CC': 26,
    '3PSB5/CG': 27,
    '3PSB5/CT': 28,
    '3PSB5/GA': 29,
    '3PSB5/GC': 30,
    '3PSB5/GG': 31,
    '3PSB5/GT': 32,
    '3PSB5/TA': 33,
    '3PSB5/TC': 34,
    '3PSB5/TG': 35,
    '3PSB5/TT': 36,
    '5PSB3/AA': 37,
    '5PSB3/AC': 38,
    '5PSB3/AG': 39,
    '5PSB3/AT': 40,
    '5PSB3/CA': 41,
    '5PSB3/CC': 42,
    '5PSB3/CG': 43,
    '5PSB3/CT': 44,
    '5PSB3/GA': 45,
    '5PSB3/GC': 46,
    '5PSB3/GG': 47,
    '5PSB3/GT': 48,
    '5PSB3/TA': 49,
    '5PSB3/TC': 50,
    '5PSB3/TG': 51,
    '5PSB3/TT': 52
}

def create_angles(seqstring):
   lseq = len(seqstring)
   angletype = []
   i1list = [] 
   i2list = [] 
   i3list = [] 
   
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
   
   iangle = 1
   
   #### SPS
   for i in range(1,lseq):
       i1 = 3*i - 2
       i2 = i1 + 2
       i3 = i1 + 3
       tipo = typeAngle['SPS/'+seqstring[i-1]+seqstring[i]]
       i1list.append(i1)
       i2list.append(i2)
       i3list.append(i3)
       angletype.append(tipo)
   for i in range(1,lseq):
       i1 = lastB1 + 3*i - 2
       i2 = i1 + 2
       i3 = i1 + 3
       tipo = typeAngle['SPS/'+WC[seqstring[lseq-1-(i-1)]]+WC[seqstring[lseq-1-i]]]
       i1list.append(i1)
       i2list.append(i2)
       i3list.append(i3)
       angletype.append(tipo)
   
   #### SBB
   for i in range(1,lseq+1):
       i1 = 3*i - 2
       i2 = i1 + 1
       i3 = lastB2 + 2 - i2
       tipo = typeAngle['SBB/'+seqstring[i-1]+WC[seqstring[i-1]]]
       i1list.append(i1)
       i2list.append(i2)
       i3list.append(i3)
       angletype.append(tipo)
   for i in range(1,lseq+1):
       i1 = lastB1 + 3*i - 2
       i2 = i1 + 1
       i3 = lastB2 + 2 - i2
       tipo = typeAngle['SBB/'+WC[seqstring[lseq-1-(i-1)]]+seqstring[lseq-1-(i-1)]]
       i1list.append(i1)
       i2list.append(i2)
       i3list.append(i3)
       angletype.append(tipo)
   
   #### 3PSB5
   for i in range(1,lseq):
       i1 = 3*i
       i2 = i1 - 2
       i3 = i1 - 1
       tipo = typeAngle['3PSB5/'+seqstring[i-1]+seqstring[i]]
       i1list.append(i1)
       i2list.append(i2)
       i3list.append(i3)
       angletype.append(tipo)
   for i in range(1,lseq):
       i1 = lastB1 + 3*i
       i2 = i1 - 2
       i3 = i1 - 1
       tipo = typeAngle['3PSB5/'+WC[seqstring[lseq-1-(i-1)]]+WC[seqstring[lseq-1-i]]]
       i1list.append(i1)
       i2list.append(i2)
       i3list.append(i3)
       angletype.append(tipo)
   
   #### 5PSB3
   for i in range(1,lseq):
       i1 = 3*i
       i2= i1 + 1
       i3= i1 + 2
       tipo = typeAngle['5PSB3/'+seqstring[i-1]+seqstring[i]]
       i1list.append(i1)
       i2list.append(i2)
       i3list.append(i3)
       angletype.append(tipo)
   for i in range(1,lseq):
       i1 = lastB1 + 3*i
       i2 = i1 + 1
       i3 = i1 + 2
       tipo = typeAngle['5PSB3/'+WC[seqstring[lseq-1-(i-1)]]+WC[seqstring[lseq-1-i]]]
       i1list.append(i1)
       i2list.append(i2)
       i3list.append(i3)
       angletype.append(tipo)

   return angletype, i1list, i2list, i3list
