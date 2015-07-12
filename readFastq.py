__author__ = 'icedevil2001'

#wget --no-check https://d28rh4a8wq0iu5.cloudfront.net/ads1/data/SRR835775_1.first1000.fastq

def readfastq(filename):
    seq_name = []
    sequence = []
    qualites = []

    with open(filename) as fh:
        while True:
            seqname= fh.readline().rstrip() # skip line
            seq = fh.readline().rstrip()
            fh.readline()
            qual = fh.readline().rstrip()
            if len(seq) == 0:
                break

            seq_name.append(seqname)
            sequence.append(seq)
            qualites.append(qual)

    return seq_name, sequence, qualites

seqnames, seqs, quals = readfastq('SRR835775_1.first1000.fastq')


def phred33ToQ(qual):
    return ord(qual) - 33




def createHist(qualityStrings):
    # Create a histogram of quality scores
    hist = [0]*50
    for read in qualityStrings:
        for phred in read:
            q = phred33ToQ(phred)
            hist[q] += 1
    return hist


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


import collections
count = collections.Counter()
for seq in seqs:
    count.update(seq)
count


h = createHist(quals)
print(h)


import matplotlib.pyplot as plt
plt.plot(range(len(h)), h)
plt.show()


gc = findGCByPos(seqs)
plt.plot(range(len(gc)), gc)
plt.show()




