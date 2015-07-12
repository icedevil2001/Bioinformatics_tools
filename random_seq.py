__author__ = 'icedevil2001'

import random

def longestprefix(s1, s2):
    i = 0
    while i < len(s1) and i < len(s2) and s1[i] == s2[i]:
        i +=1
    return s1[:i]

def match(s1,s2):
    seq == seq1
    return seq


seq = ''
for _ in range(5000):
    seq += random.choice('A')
print 'seq=', seq

seq1 =''
for _ in range(5000):
    seq1 += random.choice('ATGC')
print 'seq1=', seq1

print 'prefix=',longestprefix(seq, seq1)