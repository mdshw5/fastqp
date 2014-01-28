""" 
Classes and functions for quality assessment of FASTQ and SAM format NGS reads
"""

from __future__ import division
import sys
import os
import gzip
import itertools
import functools
import math
import matplotlib as mpl
if sys.platform is not 'darwin':
    mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.mlab as mlab
import numpy as np
import array
import random
import binascii
from collections import defaultdict, Counter
        
class Fastq(object):
    """
    A class to hold features from fastq reads.
    """
    def __init__(self, name='', seq=tuple(), strand='+', qual=tuple(), meth=None):
        self.name = name
        self.seq = seq
        self.strand = strand
        self.qual = qual
        self.meth = meth
        self.i = int()
        assert isinstance(name, str)
        assert isinstance(seq, tuple)
        assert isinstance(qual, tuple)

    def __repr__(self):
        if self.meth:
            return '\n'.join([self.name, ''.join(self.seq), self.strand, ''.join(self.qual), ''.join(self.meth)])
        else:
            return '\n'.join([self.name, ''.join(self.seq), self.strand, ''.join(self.qual)])

    def __len__(self):
        return len(self.seq)
        
class Sam(object):
    """ Store fields in each line of a SAM file, provided as a tuple. """
    def __init__(self, fields):
        self.qname = str(fields[0])
        self.flag = int(fields[1])
        self.rname = fields[2]
        self.pos = int(fields[3])
        self.mapq = int(fields[4])
        self.cigar = str(fields[5])
        self.rnext = str(fields[6])
        self.pnext = int(fields[7])
        self.tlen = int(fields[8])
        self.seq = tuple(fields[9])
        self.qual = tuple(fields[10])
        self.tags = fields[11:]
        self.insertion = []
        self.mapped = False if self.rname == '*' else True
        self.secondary = True if self.get_tag('XS') else False
        self.reverse = True if self.flag == 16 else False
        self.i = int()

    def __gt__(self, other):
        if self.rname != other.rname:
            return self.rname > other.rname
        elif (self.rname == other.rname) and (self.pos != other.pos):
            return self.pos > other.pos
        else:
            return str(self) > str(other)

    def __lt__(self, other):
        if self.rname != other.rname:
            return self.rname < other.rname
        elif (self.rname == other.rname) and (self.pos != other.pos):
            return self.pos < other.pos
        else:
            return str(self) < str(other)

    def __eq__(self, other):
        if (self.rname == other.rname) and (self.pos == other.pos) and (str(self) != str(other)):
            return str(self) == str(other)
        else:
            return self.pos == other.pos

    def __str__(self):
        return '\t'.join((self.qname, str(self.flag), self.rname, str(self.pos),
                          str(self.mapq), str(self.cigar), self.rnext, str(self.pnext),
                          str(self.tlen), ''.join(self.seq), ''.join(self.qual)) + self.tags)
    def __repr__(self):
        return self.__str__() + '\n'

    def __len__(self):
        return len(self.seq)
        
    def get_tag(self, tag):
        """ Return the data encoded in SAM tag as a string """
        assert tag == str(tag)
        tags = [t[0:2] for t in self.tags]
        try:
            return str(self.tags[tags.index(tag)])
        except ValueError:
            return None
            
class Reader:
    """ 
    Read either a fastq or sam file and return an iterator.
    """
    def __init__(self, filename, format='fastq'):
        name, ext = os.path.splitext(filename)
        assert format in ['fastq', 'sam']
        if format == 'fastq':
            self.next = self.fq_next
            self.file = open(filename, 'r')
        elif format == 'sam':
            self.header = []
            with open(filename, 'r') as samfile:
                for line in samfile:
                    if line[0] == '@':
                        self.header.append(line)
                    else:
                        break
            self.file = open(filename, 'r')
            self.header = tuple(next(self.file) for i in range(len(self.header))) ## consume header lines    
            self.next = self.sam_next
            
            
    def __next__(self):
        return self.next()
        
    def __iter__(self):
        return self
            
    def fq_next(self):
        try:
            name = next(self.file).strip()
            seq = next(self.file).strip()
            strand = next(self.file).strip()
            qual = next(self.file).strip()
            return Fastq(name=name, seq=tuple(seq), strand=strand, qual=tuple(qual))
        except StopIteration:
            raise StopIteration

    def sam_next(self):
        try:
            line = next(self.file).strip()
            return Sam(tuple(line.rstrip('\n\r').split('\t')))
        except StopIteration:
            raise StopIteration

    def subsample(self, n, k):
        """ Draws n number of reads from self, then returns the kth read.
        Supply a random k for random sampling. """
        z = None
        for i, x in zip(range(n), self):
            if i == k:
                z = x
        if z:
            return z
        else:
            return None

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.file.close()

class Stats:
    """ Statistics for characterization of NGS reads       
    """
    def __init__(self):
        self.depth = defaultdict(int)
        self.nuc = defaultdict(lambda: defaultdict(int))
        self.qual = defaultdict(lambda: defaultdict(int))
        self.gc = defaultdict(int)
        self.kmers = Counter(defaultdict(int))
            
    def evaluate(self, seq, qual):
        """ Evaluate read object at each position, and fill in nuc and qual dictionaries """
        self.gc[gc(seq)] += 1
        for i, s, q in zip(range(len(seq)), seq, qual):
            i += 1
            self.depth[i] += 1
            self.nuc[i][s[0]] += 1
            self.qual[i][q[0]] += 1
            
    def kmercount(self, seq, k=5):
        """ Count all kmers of k length in seq and update kmer counter.
        """
        for kmer in window(seq, n=k):
            self.kmers[kmer] += 1

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass
    
            
def gc(seq):
    """ Return the GC content of as an int
    >>> x = tuple('TTTTTATGGAGGTATTGAGAACGTAAGATGTTTGGATAT')
    >>> gc(x)
    30
    """
    g = seq.count('G')
    c = seq.count('C')
    return int((g + c) / len(seq) * 100)
    
def padbases(bases):
    """For each base call in dictionary D, add an entry base:0 if key base does not exist."""
    def inner(D):
        for key in bases:
            if key not in D:
                D[key] = 0
    return inner
    
def percentile(D, percent):
    """
    modified from: http://stackoverflow.com/a/2753343/717419
    Find the percentile of a list of values.

    N - is a dictionary with key=numeric value and value=frequency.
    percent - a float value from 0.0 to 1.0.
    
    outlier removal: http://www.purplemath.com/modules/boxwhisk3.htm
    
    return the percentile of the values
    """
    N = sorted(D.keys()) ## dict keys
    P = [D[n] for n in N] ## dict values
    if not N:
        return None
    k = (sum(P)) * percent
    l = (sum(P)) * 0.25 ## lower quartile
    u = (sum(P)) * 0.75 ## upper quartile
    e = int()
    for n,p in zip(N,P): ## find percentile
        e += p
        if e >= k:
            z = n ## value at percentile
            break
    e = int()
    for n,p in zip(N,P): ## find upper quartile
        e += p
        if e >= u:
            uz = n ## value at quartile
            break
    e = int()
    for n,p in zip(N,P): ## find lower quartile
        e += p
        if e >= l:
            lz = n ## value at quartile
            break
    iqd = 1.5 * (uz - lz) ## 1.5 times the inter-quartile distance
    if (z) & (z < lz - iqd):
        return int(lz - iqd)
    elif (z) & (z > uz + iqd):
        return int(uz + iqd)
    elif z:
        return int(z)
    else:
        return N[-1]
    
def qualplot(positions, quantiles, filename, fig_kw):
    Q0 = [q[0] for q in quantiles]
    Q1 = [q[1] for q in quantiles]
    Q2 = [q[2] for q in quantiles]
    Q3 = [q[3] for q in quantiles]
    Q4 = [q[4] for q in quantiles]
    fig = plt.figure(**fig_kw)
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.85)
    ax.fill_between(positions, Q4, color=(1.,.8,0.))
    ax.fill_between(positions, Q3, color=(1.,.6,0.))
    ax.fill_between(positions, Q2, color=(1.,.4,0.))
    ax.fill_between(positions, Q1, color=(1.,.2,0.))
    ax.fill_between(positions, Q0, color=(1.,1.,1.))
    ax.plot(positions, Q4, color=(1.,.8,0.))
    ax.plot(positions, Q3, color=(1.,.6,0.))
    ax.plot(positions, Q2, color=(1.,.4,0.))
    ax.plot(positions, Q1, color=(1.,.2,0.))
    ax.plot(positions, Q2, color='black')
    x1,x2,y1,y2 = ax.axis()
    ax.axis((x1,x2,0,max(Q4)))
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width, box.height])
    ax.yaxis.grid(b=True, which='major', **{'color':'gray', 'linestyle':':'})
    ax.legend(('100-75%', '75-50%', '50-25%', '25-0%', 'Median'), bbox_to_anchor=(0, 0.25), 
                loc='center left')
    ax.set_title('Quality score percentiles')
    ax.set_xlabel('Cycle')
    ax.set_ylabel('Phred score')
    plt.savefig(filename + '_quals.png')

def qualdist(qualities, filename, fig_kw):
    fig = plt.figure(**fig_kw)
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.85)
    values = map(Counter, qualities)
    counts = Counter()
    for value in values:
        counts = counts + value
    ax.plot(tuple(counts.keys()), tuple(counts.values()), color='black')
    ax.fill_between(tuple(counts.keys()), tuple(counts.values()), color=(0.1,0.6,0.8))
    ax.yaxis.grid(b=True, which='major', **{'color':'gray', 'linestyle':':'})
    ax.set_axisbelow(True)
    x1,x2,y1,y2 = ax.axis()
    ax.axis((x1,max(tuple(counts.keys())),0,y2))
    ax.set_title('Quality score cumulative distribution')
    ax.set_xlabel('Phred score')
    ax.set_ylabel('Cumulative sum')
    plt.savefig(filename + '_qualdist.png')
    
def qualmap(qualities, filename, fig_kw):
    fig = plt.figure(**fig_kw)
    ax = fig.add_subplot(111)
    values = map(Counter, tuple(qualities.values()))
    counts = Counter()
    for value in values:
        counts = counts + value
    max_qual = max(tuple(counts.keys()))
    max_pos = max(tuple(qualities.keys()))
    heat_map = np.zeros((max_qual, max_pos))
    for p in range(max_pos):
        for q in range(max_qual):
            try:
                heat_map[q][p] = qualities[p+1][q+1]
            except KeyError:
                pass
    imax = ax.imshow(np.array(heat_map), cmap=plt.cm.gist_heat, origin='lower', interpolation='none', aspect='auto')
    ax.axhline(y=10, linestyle=':', color='gray')
    ax.axhline(y=20, linestyle=':', color='gray')
    ax.axhline(y=30, linestyle=':', color='gray')
    cbar = fig.colorbar(imax, orientation='horizontal', shrink=0.5)
    cbar_labels = [item.get_text() for item in cbar.ax.get_xticklabels()]
    cbar.ax.set_xticklabels(cbar_labels, rotation=45)
    cbar.ax.set_title('')
    ax.set_title('Quality score heatmap')
    ax.set_xlabel('Cycle')
    ax.set_ylabel('Sum of Phred qualities')
    plt.savefig(filename + '_qualmap.png')
    
def nucplot(positions, nucs, counts, filename, fig_kw):
    nuc_order = ['A','T','C','G','N','M','R','W','S','Y','K','V','H','D','B']
    max_depth = sum(tuple(counts[1].values()))
    cmap = mpl.cm.get_cmap(name='Set1')
    colors = [cmap(i) for i in np.linspace(0, 1, len(nuc_order))]
    mpl.rc('axes', color_cycle=colors)
    fig, axes = plt.subplots(nrows=1, subplot_kw={'axis_bgcolor':'white'}, **fig_kw)
    for pos, count in tuple(counts.items()):
        max_depth = sum(tuple(count.values()))
        for nuc in nucs:
            counts[pos][nuc] = float(count[nuc]) / max_depth * 100
    for nuc in nuc_order:
        if nuc in nucs:
            axes.plot(positions, [count[nuc] for count in counts.values()])
    # Shink current axis by 20%
    box = axes.get_position()
    axes.set_position([box.x0, box.y0, box.width, box.height])
    axes.yaxis.grid(b=True, which='major', **{'color':'gray', 'linestyle':':'})
    axes.set_axisbelow(True)
    axes.set_title('Base content')
    axes.set_xlabel('Cycle')
    axes.set_ylabel('Base content (% basecall)')     
    legend = axes.legend(tuple((n for n in nuc_order if n in nucs)), ncol=len(nucs),
                bbox_to_anchor=(0.5,0.25), loc='center', prop={'size':8})
    frame = legend.get_frame()
    frame.set_facecolor('white')
    for label in legend.get_texts():
        label.set_color('black')
    plt.savefig(filename + '_nucs.png')
    
def depthplot(positions, depths, filename, fig_kw):
    fig, axes = plt.subplots(nrows=1, **fig_kw)
    axes.plot(positions, depths, color=(0.1,0.6,0.8))
    axes.yaxis.grid(b=True, which='major', **{'color':'gray', 'linestyle':':'})
    axes.set_title('Cumulative read lengths distribution')
    axes.set_xlabel('Cycle')
    axes.set_ylabel('Number of reads at cycle')
    plt.savefig(filename + '_depth.png')
    
def gcplot(positions, counts, filename, fig_kw):
    fig, axes = plt.subplots(nrows=1, **fig_kw)
    gc = [c['G'] + c['C'] for c in counts]
    at = [c['A'] + c['T'] for c in counts]
    tots = [x + y for x,y in zip(gc, at)]
    gcs = [float(n) / m * 100 for n,m in zip(gc, tots)]
    axes.plot(positions, gcs, color=(0.1,0.6,0.8))
    x1,x2,y1,y2 = axes.axis()
    axes.axis((x1,x2,0,100))
    axes.yaxis.grid(b=True, which='major', **{'color':'gray', 'linestyle':':'})
    axes.set_axisbelow(True)
    axes.set_title('GC content distribution')
    axes.set_xlabel('Cycle')
    axes.set_ylabel('GC% at cycle')
    plt.savefig(filename + '_gc.png')
    
def gcdist(counts, filename, fig_kw):
    m = int(sum([k * v for k,v in zip(counts.keys(), counts.values())]) / sum(counts.values()))
    variances = [(k - m)**2 for k in counts.keys()]
    variance = int(sum([k * v for k,v in zip(variances, counts.values())]) / sum(counts.values()))
    sigma = math.sqrt(variance)
    x = np.linspace(0,100,100)
    fig, axes = plt.subplots(nrows=1, **fig_kw)
    axes.plot(tuple(counts.keys()), tuple(counts.values()), color='black')
    axes.fill_between(tuple(counts.keys()), tuple(counts.values()), color=(0.1,0.6,0.8))
    x1,x2,y1,y2 = axes.axis()
    axes.axis((x1,x2,0,y2))
    axes2 = axes.twinx()
    axes2.plot(x,mlab.normpdf(x,m,sigma), color='red')
    axes2.get_yaxis().set_visible(False)
    handles, labels = axes.get_legend_handles_labels()
    display = (0,1,2)
    a = plt.Line2D((0,1),(0,0), color=(0.1,0.6,0.8))
    b = plt.Line2D((0,1),(0,0), color='red')
    legend = axes.legend([handle for i,handle in enumerate(handles) if i in display]+[a,b],
                         [label for i,label in enumerate(labels) if i in display]+['Actual','Theoretical'],
                fancybox=True, shadow=True)
    axes.yaxis.grid(b=True, which='major', **{'color':'gray', 'linestyle':':'})
    axes.set_axisbelow(True)
    axes.set_title('Read GC content distribution')
    axes.set_xlabel('Mean read GC content (%)')
    axes.set_ylabel('Number of reads')
    plt.savefig(filename + '_gcdist.png')
    
def window(seq, n=2):
    """ Returns a sliding window (of width n) over data from the iterable
    s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ... """
    it = iter(seq)
    result = tuple(itertools.islice(it, n))
    if len(result) == n:
        yield ''.join(result)
    for elem in it:
        result = result[1:] + (elem,)
        yield ''.join(result)
        
if __name__ == "__main__":
    import doctest
    doctest.testmod()
