#!/usr/bin/python

import os
import sys
import argparse
import itertools
import random
import math
import time
from functools import reduce
from fastqp.fastqp import *

def mean(s):
    return sum(s) / len(s)

def run(args):
    """ read FASTQ or SAM and tabulate basic metrics """
    time_start = time.time()
    bsize = os.path.getsize(args.input)
    
    ## estimate the number of lines in args.input
    current_entry = int()
    sample_lengths = list()
    sample_binsizes = list()
    with open(args.input, 'r') as infile:
        while current_entry < 1000:
            if args.type == 'fastq':
                header = next(infile)
                seq = next(infile)
                strand = next(infile)
                qual = next(infile)
                seq_len = len(seq)
                line = header + seq + strand + qual
            elif args.type == 'sam':
                line = next(infile)
                read = Sam(tuple(line.rstrip('\n\r').split('\t')))
                seq_len = len(read)
            blen = bin(reduce(lambda x, y: 256*x+y, (ord(c) for c in line), 0))
            sample_binsizes.append(len(blen) / 8)
            sample_lengths.append(seq_len)
            current_entry += 1
    mean_bentry = mean(sample_binsizes)
    mean_len = mean(sample_lengths)
    est_nlines = int(bsize / mean_bentry)
    act_nlines = int()
    if not args.quiet:
        sys.stderr.write("At {bytes:.0f} bytes per read of {len:.0f} length we estimate {est:,} reads in input file.\n".format(bytes=mean_bentry,
                                                                                                               len=mean_len,
                                                                                                               est=est_nlines))                                                                                                               
    ## set up factor for sampling bin size                                                                                                                   
    if args.sample:
        n = args.sample
    else:
        nf = math.floor(est_nlines / 200000)
        if nf >=1:
            n = int(nf)
        else:
            n = 1
    if not args.quiet:
        sys.stderr.write("Bin size (-s) set to {binsize:n}.\n".format(binsize=n))        
    with Reader(args.input, format=args.type) as infile, Stats() as stats:
        percent_complete = 10
        reads = itertools.islice(infile, None, None, n)
        for read in reads:
            stats.evaluate(read)
            if not args.nokmer:
                stats.kmercount(read, args.kmer)
            if not args.quiet:
                if (act_nlines / est_nlines) * 100 >= percent_complete:
                    sys.stderr.write("Approximately {0:n}% complete at read {1:,} in {2}\n".format(percent_complete, act_nlines,
                                                                                                 time.strftime('%H:%M:%S',
                                                                                                               time.gmtime(time.time()-time_start))))
                    percent_complete += 10
            act_nlines += n
            
        figures=args.figures
        positions = [k for k in sorted(stats.depth.keys())]
        depths = [stats.depth[k] for k in sorted(stats.depth.keys())]
        
        basecalls = [stats.nuc[k].keys() for k in sorted(stats.nuc.keys())]
        bases = set(list(itertools.chain.from_iterable(basecalls)))
        map(padbases(bases), stats.nuc.values())
        nbasecalls = [ '\t'.join([str(v) for v in stats.nuc[k].values()]) for k in sorted(stats.nuc.keys())]
        
        quantile_values = [0.05,0.25,0.5,0.75,0.95]
        quantiles = []
        ## replace ASCII quality with integer
        for k,v in sorted(stats.qual.items()):
            for qual in tuple(v.keys()): ## py3 keys are iterator, so build a tuple to avoid recursion
                v[ord(str(qual)) - 33] = v.pop(qual)
            line = [percentile(v, p) for p in quantile_values]
            quantiles.append(line)
            
        if not args.output:
            sys.stdout.write("{pos}\t{dep}\t{qual}\t{base}\n".format(pos='Pos',
                                                                         dep='Depth',
                                                                         base='\t'.join(bases),
                                                                         qual='\t'.join(map(str,quantile_values))))
            
            for i, position in enumerate(positions):
                sys.stdout.write("{pos}\t{dep:.1E}\t{qual}\t{base}\n".format(pos=position,
                                                            dep=depths[i],
                                                            base=nbasecalls[i],
                                                            qual='\t'.join(map(str, quantiles[i]))))
            sys.stdout.write('kmer\tcount\n')
            for key, value in stats.kmers.most_common(10):
                sys.stdout.write(key + '\t' + str(value) + '\n')
                
        elif args.output:
            with open(args.output + '_stats.txt', 'w') as out:
                out.write("{pos}\t{dep}\t{qual}\t{base}\n".format(pos='Pos',
                                                                             dep='Depth',
                                                                             base='\t'.join(bases),
                                                                             qual='\t'.join(map(str,quantile_values))))
                
                for i, position in enumerate(positions):
                    out.write("{pos}\t{dep:n}\t{qual}\t{base}\n".format(pos=position,
                                                                dep=depths[i],
                                                                base=nbasecalls[i],
                                                                qual='\t'.join(map(str, quantiles[i]))))
                                                                
                out.write('kmer\tcount\n')
                for key, value in stats.kmers.most_common(10):
                    out.write(key + '\t' + str(value) + '\n')
                    
        if args.figures:
            fig_kw = {'figsize':(8,6)}
            basename = args.output if args.output else 'plot'
            qualplot(positions, quantiles, basename, fig_kw)
            qualdist(stats.qual.values(), basename, fig_kw)
            qualmap(stats.qual, basename, fig_kw)
            depthplot(positions, depths, basename, fig_kw)
            gcplot(positions, stats.nuc.values(), basename, fig_kw)
            gcdist(stats.gc, basename, fig_kw)
            nucplot(positions, bases, stats.nuc, basename, fig_kw)
        time_finish = time.time()
        elapsed = time_finish - time_start
        if not args.quiet:
            sys.stderr.write("There were approximately {counts:,} reads in the file. Analysis finished in {sec}.\n".format(counts=act_nlines,
                                                                                                                      sec=time.strftime('%H:%M:%S',
                                                                                                                          time.gmtime(elapsed))
                                                                                                                          ))
        
def main():
    parser = argparse.ArgumentParser(prog='fastqp', description="simple NGS read quality assessment using Python")
    parser.add_argument('input', type=str, help="input file(.gz))")
    parser.add_argument('-q', '--quiet', action="store_true", default=False, help="do not print any messages (default: %(default)s)")
    parser.add_argument('-s', '--sample', type=int, help='number of reads to bin for sampling (default: auto sample 200,000 reads)')
    parser.add_argument('-k', '--kmer', type=int, default=5, choices=range(2, 11), help='length of kmer for over-repesented kmer counts (default: %(default)s)')
    parser.add_argument('-o', '--output', type=str, help="base name for output files")
    parser.add_argument('-t', '--type', type=str, default='fastq', choices=('sam', 'fastq'), help="file type for input file (default: %(default)s)")
    parser.add_argument('-f', '--figures', action="store_true", default=True, help="produce figures (default: %(default)s)")
    parser.add_argument('--nokmer', action="store_true", default=False, help="do not count kmers (default: %(default)s)")
    
    args = parser.parse_args()
    run(args)

if __name__ == "__main__": 
    main()
