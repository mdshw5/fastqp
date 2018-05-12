import random
from math import log
random.seed(42)


def qual_to_prob(Q):
    return 10**(Q / -10)


def prob_to_qual(p):
    return -10 * log(p, 10)


nucs = ('A', 'C', 'T', 'G')
with open('example.fastq', 'w') as fastq:
    for i in range(1000):
        fastq.write("@%s\n" % i)
        seq = [
            random.choice(nucs) if random.random() < 0.99 else 'N'
            for _ in range(100)
        ]
        fastq.write("".join(seq) + "\n")
        fastq.write("+\n")
        qual = [
            chr(int(prob_to_qual(1.0 - random.random()) + 33))
            if n in nucs else '!' for n in seq
        ]
        fastq.write("".join(qual) + "\n")
