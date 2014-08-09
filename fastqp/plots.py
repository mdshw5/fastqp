import sys
import numpy as np
import math
import matplotlib as mpl
if sys.platform is not 'darwin':
    mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from six.moves import map
try:
    from collections import Counter
except:
    from fastqp.backports import Counter
import itertools
from collections import defaultdict

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
    nuc_percent = defaultdict(lambda: defaultdict(int))
    for pos, count in tuple(counts.items()):
        max_depth = sum(tuple(count.values()))
        for nuc in nucs:
            if max_depth > 0:
                nuc_percent[pos][nuc] = float(count[nuc]) / max_depth * 100
            else:
                nuc_percent[pos][nuc] = 0.

    for nuc in nuc_order:
        if nuc in nucs:
            axes.plot(positions, [nuc_percent[pos][nuc] for pos in positions])
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
    cycle_gc = [sum([counts[i]['C'], counts[i]['G']]) / sum([counts[i]['C'],
                                                              counts[i]['G'],
                                                              counts[i]['A'],
                                                              counts[i]['T']]) * 100 for i in positions]
    axes.plot(positions, cycle_gc, color=(0.1,0.6,0.8))
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
    x_vals, y_vals = zip(*sorted(counts.items(), key=lambda x: x[0]))
    axes.plot(x_vals, y_vals, color='black')
    axes.fill_between(x_vals, y_vals, color=(0.1,0.6,0.8))
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
    top_strand_c = conv_dict['C']
    bot_strand_c = conv_dict['G']
    methyl_values = [(top_strand_c[pos]['C'] + bot_strand_c[pos]['G']) / (top_strand_c[pos]['Y'] + \
                                                                          bot_strand_c[pos]['R'] + \
                                                                          top_strand_c[pos]['C'] + \
                                                                          bot_strand_c[pos]['G'] + 0.1) for pos in positions]
    fig, axes = plt.subplots(nrows=1, **fig_kw)
    axes.plot(positions, methyl_values, color='red')
    x1,x2,y1,y2 = axes.axis()
    axes.axis((x1,x2,0,1))
    axes.yaxis.grid(b=True, which='major', **{'color':'gray', 'linestyle':':'})
    axes.set_axisbelow(True)
    axes.set_title('Methylation bias (M-Bias)')
    axes.set_xlabel('Cycle')
    axes.set_ylabel('Bias')
    plt.savefig(filename + '_mbias.png')


def convplot(positions, conv_dict, filename, fig_kw):
    conv_values = [(conv_dict[pos]['R'] + conv_dict[pos]['Y']) / (conv_dict[pos]['R'] + conv_dict[pos]['Y'] + conv_dict[pos]['G'] + conv_dict[pos]['C']) for pos in positions]
    fig, axes = plt.subplots(nrows=1, **fig_kw)
    axes.plot(positions, conv_values, color='gray')
    x1,x2,y1,y2 = axes.axis()
    axes.axis((x1,x2,min(conv_values)-0.1,max(conv_values)+0.1))
    axes.yaxis.grid(b=True, which='major', **{'color':'gray', 'linestyle':':'})
    axes.set_axisbelow(True)
    axes.set_title('Cytosine conversion rate')
    axes.set_xlabel('Cycle')
    axes.set_ylabel('Rate')
    plt.savefig(filename + '_conversion.png')


def kmerplot(positions, counts, filename, fig_kw):
    all_kmers = [counts[k].keys() for k in sorted(counts.keys())]
    kmers = set(list(itertools.chain.from_iterable(all_kmers)))
    kmer_len = len(tuple(kmers)[0])
    kmer_sums = Counter(dict(zip(kmers, [sum([counts[pos].get(kmer, 0) for pos in positions]) for kmer in kmers])))
    top_kmers = [x[0] for x in kmer_sums.most_common(9)]
    cmap = mpl.cm.get_cmap(name='Set1')
    colors = [cmap(i) for i in np.linspace(0, 1, len(top_kmers))]
    mpl.rc('axes', color_cycle=colors)
    fig, axes = plt.subplots(nrows=1, subplot_kw={'axis_bgcolor':'white'}, **fig_kw)
    kmer_percent = defaultdict(lambda: defaultdict(int))
    for pos, count in tuple(counts.items()):
        for kmer in kmers:
            pos_count = sum(counts[pos].values())
            if pos_count > 0:
                kmer_percent[pos][kmer] = float(count.get(kmer, 0)) / pos_count * 100
            else:
                kmer_percent[pos][kmer] = 0.
    for kmer in top_kmers:
        axes.plot(positions, [kmer_percent[pos][kmer] for pos in positions])
    for kmer in [x for x in kmers if x not in top_kmers]:
        axes.plot(positions, [kmer_percent[pos][kmer] for pos in positions], color='0.4', linestyle='dotted')
    # Shink current axis by 20%
    box = axes.get_position()
    axes.set_position([box.x0, box.y0, box.width, box.height])
    axes.yaxis.grid(b=True, which='major', **{'color':'gray', 'linestyle':':'})
    axes.set_axisbelow(True)
    axes.set_title('Kmer ({0:n}) content'.format(kmer_len))
    axes.set_xlabel('Cycle')
    axes.set_ylabel('Kmer content (% kmer)')
    legend = axes.legend(top_kmers, ncol=int(len(top_kmers)/3),
                bbox_to_anchor=(0.5,0.25), loc='center', prop={'size':8})
    frame = legend.get_frame()
    frame.set_facecolor('white')
    for label in legend.get_texts():
        label.set_color('black')
    plt.savefig(filename + '_kmers.png')
