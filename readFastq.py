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

