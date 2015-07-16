__author__ = 'icedevil2001'

#from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC


def naive(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

def reverseComplement(s):
    '''

    :param s: string of sequence ie template sequence
    :return: creates reverse complement strand

    '''
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

def readGenome(filename):
    '''

    :param filename: fasta file
    :return: reads fasta files

    '''
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

def readFastq(filename):
    '''

    :param filename: fastq file
    :return: spearates out the sequence and quals

    '''
    seqname =[]
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline().rstrip()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            #faname.append(seqname)
            sequences.append(seq)
            qualities.append(qual)


    #return sequences, qualities
    return qualities

def phred33ToQ(qual):
    return ord(qual) - 33


fastq='ERR037900_1.first1000.fastq'
print readFastq(fastq)


def createHist(qualityStrings):
    # Create a histogram of quality scores
    hist = [0]*50
    for read in qualityStrings:
        for phred in read:
            q = phred33ToQ(phred)
            hist[q] += 1
    return hist





'''
qual=readFastq(fastq)
#qual = ['#####IIIIIIII####IIIIIIIIII', 'IIIIIIIIIIIIIII']

for quals in qual:


    if createHist(quals)[:3] >= 10  :
        print quals
        print (createHist(quals))
        scores=(createHist(quals))
        print 'pos=',range(len(scores))
        #import matplotlib.pyplot as plt


        #plt.plot(range(len(scores)), scores)
        #plt.show()
'''






def findGCByPos(reads):
    ''' Find the GC ratio at each position in the read '''

    # Keep track of the number of G/C bases and the total number of bases at each position
    gc = [0] * 100
    totals = [0] * 100

    for read in reads:
        for i in range(len(read)):
            if read[i] == 'C' or read[i] == 'G':
                gc[i] += 1
            totals[i] += 1

    # Divide G/C counts by total counts to get the average at each position
    for i in range(len(gc)):
        if totals[i] > 0:
            gc[i] /= float(totals[i])

    return gc



def longestCommonPrefix(s1, s2):
    '''

    :param s1: sequence 1
    :param s2: sequence 2
    :return: find the longest Prefix

    '''
    i = 0
    while i < len(s1) and i < len(s2) and s1[i] == s2[i]:
        i += 1
    return s1[:i]
longestCommonPrefix('ACCATTG', 'ACCAAGTC')


def match(s1, s2):
    '''

    :param s1: sequence 1
    :param s2: sequence 2
    :return: find the match
    '''
    if not len(s1) == len(s2):
        return False
    for i in range(0, len(s1)):
        if not s1 == s2:
            return False
    return True

#match('ACCATTG', 'ACCATTG')

def countbases(seq):
    '''

    :param seq:
    :return: Count the number of occurences of each base
    '''
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    for base in seq:
        counts[base] +=1
    return (counts)

def occurencesbases(seq):
    '''

    :param seq:
    :return: Count the number of occurences of each base

    '''
    import collections
    counts=collections.Counter(seq)
    return counts


def naive_2mm(p,t):
    '''
    doesnt work at the moment
    :p: patten or quary sequence
    :t: target or genome sequence
    '''
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            print
            if t[i+j] > p[j]: # compare characters
                #match=False

            #elif t[i+j+1] !=p[j]:

                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences



#print naive_2mm('ACTTA', 'ACTTAACTTAGATAAAGT')

def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's suffx in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

import itertools



def scs(ss):
    """ Returns shortest common superstring of given strings,
        assuming no string is a strict substring of another """
    shortest_sup = None
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            #sup += ssperm[i+1][-(len(ssperm[i+1])-olen):]
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup  # found shorter superstring
    return shortest_sup  # return shortest



def pick_maximal_overlap(reads, k):
    """ Return a pair of reads from the list with a
        maximal suffix/prefix overlap >= k.  Returns
        overlap length 0 if there are no such overlaps."""
    reada, readb = None, None
    best_olen = 0
    for a, b in itertools.permutations(reads, 2):
        olen = overlap(a, b, min_length=k)
        if olen > best_olen:
            reada, readb = a, b
            best_olen = olen
    return reada, readb, best_olen



def greedy_scs(reads, k):
    """ Greedy shortest-common-superstring merge.
        Repeat until no edges (overlaps of length >= k)
        remain. """
    read_a, read_b, olen = pick_maximal_overlap(reads, k)
    while olen > 0:
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a + read_b[olen:])
        read_a, read_b, olen = pick_maximal_overlap(reads, k)
    return ''.join(reads)



print 'greedy', overlap('AGGGGGCATTCGATCGATGTACGTCGTAGCTAGTAGCCTG', 'CTAGTAGCCTGGTCAGTCGTAGTGTAGCTGCTAAGTCAT', 4)




#lambda_genome = readGenome('lambda_virus.fa')

#seq = 'AGTCAGGAGTCTAGTCGAGGTGTCGAGGTGGTCAGATCGATGACTAC'
#kmer = 'TTAA'



#temp=naive(kmer,lambda_genome)
#comp=naive(reverseComplement(kmer),lambda_genome)
#print 'kmer=', kmer, 'and RC=',reverseComplement(kmer)
#print 'temp=', len(temp)
#print 'RC=', len(comp)
'''


'''
#comp= (naive(reverseComplement(kmer),lambda_genome))
#print seq.reverse_complement
#temp= (naive(kmer,lambda_genome))
#print 'kmer=',kmer, 'RC_kmer=', reverseComplement(kmer)
#print 'temp=', temp[0]
#print 'RC=', comp[0]
#print 'total=',len(temp)+len(comp)
