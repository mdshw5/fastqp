""" 
Classes and functions for quality assessment of FASTQ and SAM format NGS reads
"""

from __future__ import division
import sys
import os
import gzip
import itertools
import functools
import plots
import strings
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import random
from collections import defaultdict, Counter
        
class read(object):
    """
    A class to hold features from NGS reads. 
    """
    __slots__ = ['name', 'seq', 'strand', 'qual']
    def __init__(self, name='', seq='', strand='', qual=''):
        self.name = name
        self.seq = seq
        self.strand = strand
        self.qual = qual
        
    def __iter__(self):
        for seq, qual in itertools.izip(self.seq, self.qual):
            yield self.__class__(self.name, seq, self.strand, qual)
                
    def __getitem__(self, key):
        return self.__class__(self.name, self.seq[key], self.strand, self.qual[key])       

    def __repr__(self):
        return '\n'.join([self.name, self.seq, self.strand, self.qual])
            
    def __len__(self):
        return len(self.seq)

    def complement(self):
        """ Returns the compliment of self. This only affects the sequence slot. """
        dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
        compseq = ''.join(map(lambda x: dict[x], self.seq))
        return self.__class__(self.name, compseq, self.strand, self.qual)
        
    def pair(self):
        """ Return 0 1 or 2 depending on whether read is 1 or 2 in a pair, or unpaired (0). """
        n = self.name[-2:]
        if n[0] != '/':
            return 0
        else:
            return int(n[1])
            
    def gc(self):
        """ Return the GC content of self as an int """
        g = self.seq.count('G')
        c = self.seq.count('C')
        return int((g + c) / len(self) * 100)
    
    def cpg(self):
        """ Return the number of CpG sites in self.seq """
        return self.seq.count('CG')

class reader:
    """ 
    A class to read the name, sequence, strand and qualities from a fastq file
    """
    def __init__(self, filename, format='fastq'):
        name, ext = os.path.splitext(filename)
        if ext == '.gz':
            self.file = gzip.open(filename, 'rb')
        else:
            self.file = open(filename, 'r')
        assert format in ['fastq', 'sam']
        if format == 'fastq':
            self.__iter__ = self._iterfq
        elif format == 'sam':
            self.__iter__ = self._itersam
            
    def _iterfq(self):
        """ 
        Return :py:class:`read` object.
        """
        for i, line in enumerate(self.file):
            if i % 4 == 0:
                name = line.split()[0].strip()[1:]
            elif i % 4 == 1:
                sequence = line.strip()
            elif i % 4 == 2:
                strand = line.strip()
            elif i % 4 == 3:
                qualities = line.rstrip('\n\r')
                yield read(name, sequence, strand, qualities)

    def _itersam(self):
        """ 
        Return :py:class:`read` object.
        """
        for i, line in enumerate(self.file):
            if line[0] == '@':
                self.header.append(line.rstrip('\n\r'))
            else:
                fields = tuple(line.rstrip('\n\r').split('\t'))
                yield read(fields[0], fields[9], '+', fields[10])

    def subsample(self, n, k):
        """ Draws n number of reads from self, then returns the kth read. 
        Supply a random k for random sampling. """
        z = None
        for i, x in itertools.izip(xrange(n), self):
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

class stats:
    """ Statistics for characterization of GEMINI reads       
    """
    def __init__(self):
        self.depth = defaultdict(int)
        self.nuc = defaultdict(lambda: defaultdict(int))
        self.qual = defaultdict(lambda: defaultdict(int))
        self.gc = defaultdict(int)
        self.kmers = Counter(defaultdict(int))
            
    def evaluate(self, read):
        """ Evaluate read object at each position, and fill in nuc and qual dictionaries """
        self.gc[read.gc()] += 1
        self.kmerdict(read.seq, self.kmers, k=5)
        for i,r in enumerate(read):
            i += 1
            self.depth[i] += 1
            self.nuc[i][r.seq] += 1
            self.qual[i][r.qual] += 1
    
    @staticmethod            
    def kmerdict(string, D, k):
        """ Count all kmers of k length in string and return a dictionary.
        D: a defaultdict(int) to update
        >>> D = defaultdict(int)
        >>> stats.kmerdict('ACTGTGCATGTACTGTACGTGGA', D, 22)
        >>> D
        defaultdict(<type 'int'>, {'CTGTGCATGTACTGTACGTGGA': 1, 'ACTGTGCATGTACTGTACGTGG': 1})
        """
        for kmer in strings.window(string, n=k):
            D[kmer] += 1
                
    @staticmethod    
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
    
    @staticmethod
    def padbases(bases):
        """For each base call in dictionary D, add an entry base:0 if key base does not exist."""
        def inner(D):
            for key in bases:
                if key not in D:
                    D[key] = 0
        return inner
        
    def summarize(self, filename=None, figures=True):
        positions = [k for k in sorted(self.depth.iterkeys())]
        depths = [self.depth[k] for k in sorted(self.depth.iterkeys())]
        
        basecalls = [self.nuc[k].keys() for k in sorted(self.nuc.iterkeys())]
        bases = set(list(itertools.chain.from_iterable(basecalls)))
        _padbases = self.padbases(bases)
        map(_padbases, self.nuc.itervalues())
        nbasecalls = [ '\t'.join([str(v) for v in self.nuc[k].values()]) for k in sorted(self.nuc.iterkeys())]
        
        quantile_values = [0.05,0.25,0.5,0.75,0.95]
        quantiles = []
        for k,v in sorted(self.qual.iteritems()):
            for qual in v.keys(): ## replace ASCII quality with integer
                v[ord(str(qual)) - 33] = v.pop(qual)
            line = [self.percentile(v, p) for p in quantile_values]
            quantiles.append(line)
            
        sys.stdout.write(str(self.kmers.most_common(10)).strip('[]'))
            
        if not file:
            sys.stdout.write("{pos}\t{dep}\t{qual}\t{base}\n".format(pos='Pos',
                                                                         dep='Depth',
                                                                         base='\t'.join(bases),
                                                                         qual='\t'.join(map(str,quantile_values))))
            
            for i, position in enumerate(positions):
                sys.stdout.write("{pos}\t{dep:.1E}\t{qual}\t{base}\n".format(pos=position,
                                                            dep=depths[i],
                                                            base=nbasecalls[i],
                                                            qual='\t'.join(map(str, quantiles[i]))))
                
        elif filename:
            with open(filename + 'stats.txt', 'w') as out:
                out.write("{pos}\t{dep}\t{qual}\t{base}\n".format(pos='Pos',
                                                                             dep='Depth',
                                                                             base='\t'.join(bases),
                                                                             qual='\t'.join(map(str,quantile_values))))
                
                for i, position in enumerate(positions):
                    out.write("{pos}\t{dep:.1E}\t{qual}\t{base}\n".format(pos=position,
                                                                dep=depths[i],
                                                                base=nbasecalls[i],
                                                                qual='\t'.join(map(str, quantiles[i]))))
                    
            if figures:
                fig_kw = {'figsize':(10,6)}
                qualplot(positions, quantiles, filename, fig_kw)
                qualdist(self.qual.values(), filename, fig_kw)
                depthplot(positions, depths, filename, fig_kw)
                gcplot(positions, self.nuc.values(), filename, fig_kw)
                gcdist(self.gc, filename, fig_kw)
                nucplot(positions, bases, self.nuc.values(), filename, fig_kw)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass
        
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
    ax.plot(positions, Q0, color=(1.,1.,1.))
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
    ax.legend(('100%', '75%', '50%', '25%', '0%'), bbox_to_anchor=(1, 0.5), 
                fancybox=True, shadow=True, loc='center left')
    ax.set_title('Quality score percentiles by position in reads')
    ax.set_xlabel('position')
    ax.set_ylabel('Phred quality ')
    plt.savefig(filename + '_quals.png')

def qualdist(qualities, filename, fig_kw):
    fig = plt.figure(**fig_kw)
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.85)
    values = map(Counter, qualities)
    counts = Counter()
    for value in values:
        counts = counts + value
    ax.fill_between(counts.keys(), counts.values(), color=(0.1,0.6,0.8))
    ax.set_title('Quality score distribution')
    ax.set_xlabel('quality score')
    ax.set_ylabel('count')
    plt.savefig(filename + '_qualdist.png')
    
def nucplot(positions, nucs, counts, filename, fig_kw):
    max_depth = sum(counts[0].values())
    colors = [(.5,.5,.5),
              (1,0,0),
              (0,1,0),
              (0,0,1),
              (1,1,0),
              (0,1,1),
              (1,0,1),
              (1,1,1),
              (.5,0,0),
              (0,.5,0)]
    mpl.rc('axes', color_cycle=colors)
    fig, axes = plt.subplots(nrows=1, subplot_kw={'axis_bgcolor':'black'}, **fig_kw)
    for i,count in enumerate(counts):
        max_depth = sum(count.values())
        for nuc in nucs:
            counts[i][nuc] = float(count[nuc]) / max_depth * 100
    for nuc in nucs:
        axes.plot(positions, [count[nuc] for count in counts])
    # Shink current axis by 20%
    box = axes.get_position()
    axes.set_position([box.x0, box.y0, box.width * 0.9, box.height])
    axes.set_title('basecall% by position in read')
    axes.set_xlabel('position')
    axes.set_ylabel('basecall percent')
    legend = axes.legend(tuple(nucs), bbox_to_anchor=(1, 0.5), 
                fancybox=True, shadow=True, loc='center left')
    frame = legend.get_frame()
    frame.set_facecolor('black')
    for label in legend.get_texts():
        label.set_color('white')
    plt.savefig(filename + '_nucs.png')
    
def depthplot(positions, depths, filename, fig_kw):
    fig, axes = plt.subplots(nrows=1, **fig_kw)
    axes.fill_between(positions, depths, color=(0.1,0.6,0.8))
    axes.set_title('coverage of each position in read')
    axes.set_xlabel('position')
    axes.set_ylabel('depth')
    plt.savefig(filename + '_depth.png')
    
def gcplot(positions, counts, filename, fig_kw):
    fig, axes = plt.subplots(nrows=1, **fig_kw)
    gc = [c['G'] + c['C'] for c in counts]
    at = [c['A'] + c['T'] for c in counts]
    tots = [x + y for x,y in zip(gc, at)]
    gcs = [float(n) / m * 100 for n,m in zip(gc, tots)]
    axes.plot(positions, gcs)
    axes.set_title('GC content at each position in read')
    axes.set_xlabel('position')
    axes.set_ylabel('gc')
    plt.savefig(filename + '_gc.png')
    
def gcdist(counts, filename, fig_kw):
    fig, axes = plt.subplots(nrows=1, **fig_kw)
    axes.fill_between(counts.keys(), counts.values(), color=(0.1,0.6,0.8))
    axes.set_title('GC content distribution')
    axes.set_xlabel('GC')
    axes.set_ylabel('counts')
    plt.savefig(filename + '_gcdist.png')
        
if __name__ == "__main__":
    import doctest
    doctest.testmod()