"""
Classes and functions for quality assessment of FASTQ and SAM format NGS reads
"""

from __future__ import division
import sys
import os
import re
from six.moves import range, zip, map
from six import string_types
import math
import matplotlib as mpl
from itertools import groupby, islice
if sys.platform is not 'darwin':
    mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
from collections import defaultdict
try:
    from collections import Counter
except:
    from fastqp.backports import Counter
from subprocess import Popen, PIPE
from io import TextIOWrapper


class Gzip(object):
    """ Call system gzip and maintain interface compatibility with python
    gzip module """
    def __init__(self, filename, mode):
        self.stream, self.p = self.open(filename, mode)
        self.mode = mode
        self.filename = filename

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def next(self):
        return next(self.stream)

    def open(self, filename, mode):
        if 'r' in mode:
            self.fh = open(filename, 'rb', 0)
            p = Popen(['gzip', '-dc', filename], stdout=PIPE)
            if 'b' in mode:
                fh = p.stdout
            else:
                fh = TextIOWrapper(p.stdout)
        elif 'w' in mode:
            self.fh = open(filename, 'wb', 0)
            p = Popen(['gzip', '-c'], stdin=PIPE, stdout=self.fh)
            fh = p.stdin
        return (fh, p)

    def write(self, string):
        self.stream.write(string.encode('utf-8'))

    def read(self, string):
        self.stream.read(string)

    def close(self):
        self.__exit__()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.p.communicate()
        if self.fh:
            self.fh.close()


class Fastq(object):
    """
    A class to hold features from fastq reads.
    """
    def __init__(self, name='', seq='', strand='+', qual='', conv=None):
        self.name = name
        self.seq = seq
        self.strand = strand
        self.qual = qual
        self.conv = conv
        self.i = int()
        assert isinstance(name, string_types)
        assert isinstance(seq, string_types)
        assert isinstance(qual, string_types)

    def __iter__(self):
        return self

    def next(self):
        if self.i < len(self):
            value, self.i = self[self.i], self.i + 1
            return value
        else:
            raise StopIteration()

    def __getitem__(self, key):
        if self.conv:
            return self.__class__(self.name, self.seq[key], self.strand,
                                  self.qual[key], self.conv[key])
        else:
            return self.__class__(self.name, self.seq[key], self.strand,
                                  self.qual[key])

    def __next__(self):
        return self.next()

    def __repr__(self):
        return str(self)

    def __str__(self):
        if self.name[0] != '@':
            self.name = ''.join(['@', self.name])
        if self.conv:
            return '\n'.join(['{0}:YM:Z:{1}'.format(self.name, self.conv),
                             self.seq, self.strand, self.qual]) + '\n'
        else:
            return '\n'.join([self.name, self.seq, self.strand, self.qual]) + '\n'

    def __len__(self):
        return len(self.seq)

    def gc(self):
        """ Return the GC content of self as an int
        >>> x = Fastq(name='test', seq='TTTTTATGGAGGTATTGAGAACGTAAGATGTTTGGATAT', qual=' # # ##EC4<?4A<+EFB@GHC<9FAA+DDCAFFC=22')
        >>> x.gc()
        30
        """
        g = self.seq.count('G')
        c = self.seq.count('C')
        return int((g + c) / len(self) * 100)


class Sam(object):
    """ Store fields in each line of a SAM file, provided as a tuple. """
    __slots__ = ['qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual', 'tags', '_tags', '_cigars']

    def __init__(self, fields):
        self.qname = fields[0]
        self.flag = int(fields[1])
        self.rname = fields[2]
        self.pos = int(fields[3])
        self.mapq = int(fields[4])
        self.cigar = fields[5]
        self.rnext = fields[6]
        self.pnext = int(fields[7])
        self.tlen = int(fields[8])
        self.seq = fields[9]
        self.qual = fields[10]
        self.tags = None
        self._tags = fields[11:]
        self._cigars = None

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
        if not self.tags:
            self.tags = parse_sam_tags(self._tags)
        return '\t'.join((self.qname, str(self.flag), self.rname, str(self.pos),
                          str(self.mapq), str(self.cigar), self.rnext, str(self.pnext),
                          str(self.tlen), ''.join(self.seq), ''.join(self.qual)) + \
                          tuple(':'.join((tag, self.tags[tag][0], str(self.tags[tag][1]))) for tag in sorted(self.tags.keys()))) + '\n'
    def __repr__(self):
        return "Sam({0}:{1}:{2})".format(self.rname, self.pos, self.qname)

    def __len__(self):
        return sum(c[0] for c in self.cigars if c[1] in
                   ("M", "D", "N", "EQ", "X", "P"))

    def __getitem__(self, key):
        if not self.tags:
            self.tags = parse_sam_tags(self._tags)
        return self.tags[key][1]

    def __setitem__(self, key, value):
        if not self.tags:
            self.tags = parse_sam_tags(self._tags)
        self.tags[key] = value

    def cigar_split(self):
        """ CIGAR grouping function modified from:
        https://github.com/brentp/bwa-meth
        """
        if self.cigar == "*":
            yield (0, None)
            raise StopIteration
        cig_iter = groupby(self.cigar, lambda c: c.isdigit())
        for g, n in cig_iter:
            yield int("".join(n)), "".join(next(cig_iter)[1])

    @property
    def conv(self):
        return self['YM']

    @property
    def cigars(self):
        if not self._cigars:
            self._cigars = tuple(self.cigar_split())
        return self._cigars

    @property
    def mapped(self):
        return not (self.flag & 0x4)

    @property
    def secondary(self):
        return bool(self.flag & 0x100)

    @property
    def reverse(self):
        return bool(self.flag & 0x10)

    @property
    def duplicate(self):
        return bool(self.flag & 0x400)

    def gapped(self, string, gap_char='-'):
        """ Return string with all deletions wrt reference
         represented as gaps '-' and all insertions wrt reference
         removed.
        i: sequence index
        """
        seq = []
        i = 0
        for n, t in self.cigars:
            if t in ("M", "N", "EQ", "X", "P"):
                seq.extend(string[i:i+n])
                i += n
            elif t in ("D",):
                seq.extend(('-',) * n)
            elif t in ("I",):
                i += n
        return ''.join(seq)

    @property
    def coords(self):
        return range(self.pos, self.pos + len(self))


class FastqReader:
    """
    A class to read the name, sequence, strand and qualities from a fastq file
    """
    def __init__(self, f):
        name, ext = os.path.splitext(f.name)
        if ext == '.gz':
            self.file = Gzip(''.join([name, ext]), 'r')
        else:
            self.file = f

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def next(self):
        try:
            name = next(self.file).strip().split()[0]  # remove whitespace
            seq = next(self.file).strip()
            strand = next(self.file).strip()
            qual = next(self.file).strip()
            if name.count(':YM:Z:') > 0:
                tag, dtype, data = name.split(':')[-3:]
                name = ':'.join(name.split(':')[:-3])
                return Fastq(name=name, seq=seq, strand=strand, qual=qual, conv=data)
            else:
                return Fastq(name=name, seq=seq, strand=strand, qual=qual)
        except StopIteration:
            raise StopIteration

    def subsample(self, n):
        """ Draws every nth read from self. Returns Fastq. """
        n = n * 4
        for i, line in enumerate(self.file):
            if i % n == 0:
                name = line.strip().split()[0]
            elif i % n == 1:
                seq = line.strip()
            elif i % n == 2:
                strand = line.strip()
            elif i % n == 3:
                qual = line.strip()
                if name.count(':YM:Z:') > 0:
                    tag, dtype, data = name.split(':')[-3:]
                    name = ':'.join(name.split(':')[:-3])
                    yield Fastq(name=name, seq=seq, strand=strand, qual=qual, conv=data)
                else:
                    yield Fastq(name=name, seq=seq, strand=strand, qual=qual)

    def fileno(self):
        return self.file.fileno()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.file.close()


class Reader(object):
    """ Read SAM/BAM format file using iterator. """
    def __init__(self, f):
        name, ext = os.path.splitext(f.name)
        if ext == '.bam':
            BamReaderSamtools.__init__(self, f)
        else:
            SamReader.__init__(self, f)

    def next(self):
        try:
            line = next(self.file).rstrip('\n\r')
            return Sam(tuple(line.split('\t')))
        except StopIteration:
            raise StopIteration

    def __next__(self):
        return self.next()

    def __iter__(self):
        return self

    def subsample(self, n):
        """ Draws every nth read from self. Returns Sam. """
        for i, line in enumerate(self.file):
            if i % n == 0:
                yield Sam(tuple(line.rstrip().split('\t')))

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.file.close()


class SamReader(Reader):
    """ Read SAM format file using iterator. """
    def __init__(self, f):
        self.header = []
        self.file = f
        for line in self.file:
            if line[0] == '@':
                self.header.append(line.rstrip('\n\r'))
            else:
                break


class BamReaderSamtools(Reader):
    """ Read BAM format file using iterator. """
    def __init__(self, f):
        pline = ['samtools', 'view', '-H', f.name]
        try:
            p = Popen(pline, bufsize=-1, stdout=PIPE,
                      stderr=PIPE)
        except OSError:
            sys.stderr.write('Samtools must be installed for BAM file support!\n')
            sys.exit(1)
        self.header = [line.decode('utf-8').rstrip('\n\r') for line in p.stdout]
        p.wait()
        pline = ['samtools', 'view', f.name]
        self.p = Popen(pline, bufsize=-1, stdout=PIPE,
                       stderr=PIPE)
        self.file = TextIOWrapper(self.p.stdout)

    def __exit__(self, *args):
        self.p.wait()


class Stats:
    """ Counter for characterization of NGS reads
    """
    def __init__(self):
        self.depth = defaultdict(int)
        self.nuc = defaultdict(lambda: defaultdict(int))
        self.qual = defaultdict(lambda: defaultdict(int))
        self.gc = defaultdict(int)
        self.kmers = Counter(defaultdict(int))
        self.conv = defaultdict(lambda: defaultdict(int))

    def evaluate(self, seq, qual, conv=None):
        """ Evaluate read object at each position, and fill in nuc and qual dictionaries """
        self.gc[gc(seq)] += 1
        if conv:
            cpgs = cpg_map(seq)
        for i in range(1, len(seq) + 1):
            self.depth[i] += 1
            self.nuc[i][seq[i-1]] += 1
            self.qual[i][qual[i-1]] += 1
            if conv:
                if cpgs[i-1] != 'N':
                    self.conv[i][conv[i-1]] += 1

    def kmercount(self, seq, k=5):
        """ Count all kmers of k length in seq and update kmer counter.
        """
        for kmer in window(seq, n=k):
            self.kmers[kmer] += 1

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass


def bam_read_count(bamfile):
    """ Return a tuple of the number of mapped and unmapped reads in a bam file """
    p = Popen(['samtools', 'idxstats', bamfile], stdout=PIPE)
    mapped = 0
    unmapped = 0
    for line in p.stdout:
        rname, rlen, nm, nu = line.rstrip().split()
        mapped += int(nm)
        unmapped += int(nu)
    return (mapped, unmapped)


def parse_sam_tags(tagfields):
    """ Return a dictionary containing the tags """
    return dict((tag, (dtype, data)) for tag, dtype, data in (decode_tag(x) for x in tagfields))


def encode_tag(tag, data_type, data):
    """ Write a SAM tag in the format ``TAG:TYPE:data``
    >>> encode_tag('YM', 'Z', '#""9O"1@!J')
    'YM:Z:#""9O"1@!J'
    """
    value = ':'.join(list((tag.upper(), data_type.upper(), data)))
    return value


def decode_tag(tag):
    """ Parse a SAM format tag to a (TAG, TYPE, data) tuple.

    TYPE in A, i, f, Z, H, B

    >>> decode_tag('YM:Z:#""9O"1@!J')
    ('YM', 'Z', '#""9O"1@!J')
    >>> decode_tag('XS:i:5')
    ('XS', 'i', 5)
    """
    values = tag.split(':')
    if values[1] == 'i':
        values[2] = int(values[2])
    elif values[1] == 'f':
        values[2] = float(values[2])
    elif values[1] == 'H':
        raise(NotImplementedError, "Hex array SAM tags are currently not parsed.")
    elif values[1] == 'B':
        raise(NotImplementedError, "Byte array SAM tags are currently not parsed.")
    return tuple(values)


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
    N = sorted(D.keys())  # dict keys
    P = [D[n] for n in N]  # dict values
    if not N:
        return None
    k = (sum(P)) * percent
    l = (sum(P)) * 0.25  # lower quartile
    u = (sum(P)) * 0.75  # upper quartile
    e = int()
    for n,p in zip(N, P):  # find percentile
        e += p
        if e >= k:
            z = n  # value at percentile
            break
    e = int()
    for n,p in zip(N, P):  # find upper quartile
        e += p
        if e >= u:
            uz = n  # value at quartile
            break
    e = int()
    for n,p in zip(N, P):  # find lower quartile
        e += p
        if e >= l:
            lz = n  # value at quartile
            break
    iqd = 1.5 * (uz - lz)  # 1.5 times the inter-quartile distance
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
    axes.legend([handle for i,handle in enumerate(handles) if i in display]+[a,b],
                         [label for i,label in enumerate(labels) if i in display]+['Actual','Theoretical'])
    axes.yaxis.grid(b=True, which='major', **{'color':'gray', 'linestyle':':'})
    axes.set_axisbelow(True)
    axes.set_title('Read GC content distribution')
    axes.set_xlabel('Mean read GC content (%)')
    axes.set_ylabel('Number of reads')
    plt.savefig(filename + '_gcdist.png')


def mbiasplot(positions, conv_dict, filename, fig_kw):
    methyl_values = [(conv_dict[pos]['G'] + conv_dict[pos]['C']) / (conv_dict[pos]['R'] + conv_dict[pos]['Y'] + conv_dict[pos]['G'] + conv_dict[pos]['C']) for pos in positions]
    fig, axes = plt.subplots(nrows=1, **fig_kw)
    axes.plot(positions, methyl_values, color='red')
    x1,x2,y1,y2 = axes.axis()
    axes.axis((x1,positions[-1],0,1))
    axes.yaxis.grid(b=True, which='major', **{'color':'gray', 'linestyle':':'})
    axes.set_axisbelow(True)
    axes.set_title('Methylation bias (M-Bias)')
    axes.set_xlabel('Cycle')
    axes.set_ylabel('Bias')
    plt.savefig(filename + '_mbias.png')


def window(seq, n=2):
    """ Returns a sliding window (of width n) over data from the iterable
    s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ... """
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield ''.join(result)
    for elem in it:
        result = result[1:] + (elem,)
        yield ''.join(result)


def mean(s):
    return sum(s) / len(s)


def cpg_map(seq):
    """ Return tuple of C/G/N.

    >>> cpg_map('CGCGTAGCCG')
    'CGCGNNNNCG'
    """
    starts = (x.start() for x in re.finditer('CG', ''.join(['N', seq, 'N'])))
    cpgs = ['N'] * (len(seq) + 2)
    for start in starts:
        cpgs[start] = 'C'
        cpgs[start+1] = 'G'
    return ''.join(cpgs[1:-1])


if __name__ == "__main__":
    import doctest
    doctest.testmod()
