#!/usr/bin/python

import os
import sys
import argparse
import itertools
import math
import time
from fastqp import FastqReader, Reader, padbases, percentile, qualplot, \
    qualdist, qualmap, nucplot, depthplot, gcplot, gcdist, mbiasplot, \
    window, mean, cpg_map, bam_read_count, gc
from collections import defaultdict

def run(args):
    """ read FASTQ or SAM and tabulate basic metrics """
    time_start = time.time()
    bsize = os.path.getsize(args.input.name)

    est_counter = int()
    sample_lengths = list()
    sample_binsizes = list()
    act_nlines = int()
    name, ext = os.path.splitext(args.input.name)

    if ext not in ['.fastq', '.sam', '.bam', '.gz']:
        sys.exit("Input file must end in either .sam, .bam, .fastq, or .fastq.gz\n")
    # estimate the number of lines in args.input
    if ext in ['.fastq']:
        with FastqReader(open(args.input.name)) as fh:
            for read in fh:
                sample_lengths.append(len(read))
                sample_binsizes.append(len(str(read)))
                est_counter += 1
                if est_counter == 10000:
                    break
            mean_bentry = mean(sample_binsizes)
            mean_len = mean(sample_lengths)
            est_nlines = int(bsize / mean_bentry)
            if not args.quiet:
                sys.stderr.write("At {bytes:.0f} bytes per read of {len:.0f} length "
                "we estimate {est:,} reads in input file.\n".format(bytes=mean_bentry,
                                                                    len=mean_len,
                                                                    est=est_nlines))
    elif ext  == '.sam':
        with Reader(open(args.input.name)) as fh:
            for read in fh:
                sample_lengths.append(len(read))
                sample_binsizes.append(len(str(read)))
                est_counter += 1
                if est_counter == 10000:
                    break
            mean_bentry = mean(sample_binsizes)
            mean_len = mean(sample_lengths)
            est_nlines = int(bsize / mean_bentry)
            if not args.quiet:
                sys.stderr.write("At {bytes:.0f} bytes per read of {len:.0f} length "
                "we estimate {est:,} reads in input file.\n".format(bytes=mean_bentry,
                                                                    len=mean_len,
                                                                    est=est_nlines))
    elif ext == '.bam':
        est_nlines = sum(bam_read_count(args.input.name))
        if not args.quiet:
            sys.stderr.write("{est:,} reads in input file.\n".format(est=est_nlines))
    elif ext == '.gz':
        if args.binsize:
            n = args.binsize
        else:
            n = 1
        if not args.quiet:
            sys.stderr.write("Gzipped file detected, bin size (-s) set to {binsize:n}.\n".format(binsize=n))

    if ext != '.gz':
        # set up factor for sampling bin size
        if args.binsize:
            n = args.binsize
        else:
            nf = math.floor(est_nlines / args.nreads)
            if nf >= 1:
                n = int(nf)
            else:
                n = 1
        if not args.quiet:
            sys.stderr.write("Bin size (-s) set to {binsize:n}.\n".format(binsize=n))

    if ext in ['.sam', '.bam']:
        infile = Reader(args.input)
    else:
        infile = FastqReader(args.input)

    cycle_depth = defaultdict(int)
    cycle_nuc = defaultdict(lambda: defaultdict(int))
    cycle_qual = defaultdict(lambda: defaultdict(int))
    cycle_gc = defaultdict(int)
    cycle_kmers = defaultdict(lambda: defaultdict(int))
    cycle_conv = defaultdict(lambda: defaultdict(int))
    percent_complete = 10
    reads = infile.subsample(n)

    for read in reads:
        if ext in ['.sam', '.bam']:
            if args.aligned and not read.mapped:
                continue
            elif args.unaligned and read.mapped:
                continue
            if read.reverse:
                if args.mbias:
                    conv = read.conv[::-1]
                else:
                    conv = read.seq[::-1]
                seq = read.seq[::-1]
                qual = read.qual[::-1]
            else:
                if args.mbias:
                    conv = read.conv
                else:
                    conv = read.seq
                seq = read.seq
                qual = read.qual
        else:
            if args.mbias:
                conv = read.conv
            seq = read.seq
            qual = read.qual

        cycle_gc[gc(seq)] += 1
        cpgs = cpg_map(seq)
        for i, (s, q, c, p) in enumerate(zip(seq, qual, conv, cpgs)):
            cycle_depth[i+1] += 1
            cycle_nuc[i+1][s] += 1
            cycle_qual[i+1][q] += 1
            if p != 'N':
                cycle_conv[i+1][c] += 1

        if not args.nokmer:
            for i, kmer in enumerate(window(seq, n=3)):
                cycle_kmers[i][kmer] += 1

        if not args.quiet and ext != '.gz':
            if (act_nlines / est_nlines) * 100 >= percent_complete:
                sys.stderr.write("Approximately {0:n}% complete at "
                                 "read {1:,} in {2}\n".format(percent_complete,
                                                              act_nlines,
                                                              time.strftime('%H:%M:%S',
                                                                            time.gmtime(time.time()-time_start))))
                percent_complete += 10
        act_nlines += n

    positions = [k for k in sorted(cycle_depth.keys())]
    depths = [cycle_depth[k] for k in sorted(cycle_depth.keys())]

    basecalls = [cycle_nuc[k].keys() for k in sorted(cycle_nuc.keys())]
    bases = set(list(itertools.chain.from_iterable(basecalls)))
    nbasecalls = [ '\t'.join([str(cycle_nuc[p].get(k, 0)) for k in bases]) for p in sorted(cycle_nuc.keys())]
    map(padbases(bases), cycle_nuc.values())

    quantile_values = [0.05,0.25,0.5,0.75,0.95]
    quantiles = []
    ## replace ASCII quality with integer
    for k,v in sorted(cycle_qual.items()):
        for q in tuple(v.keys()): ## py3 keys are iterator, so build a tuple to avoid recursion
            v[ord(str(q)) - 33] = v.pop(q)
        line = [percentile(v, p) for p in quantile_values]
        quantiles.append(line)

    out = open(args.output + '_counts.txt', 'w')
    out.write("{pos}\t{dep}\t{qual}\t{base}\n".format(pos='Pos',
                                                                 dep='Depth',
                                                                 base='\t'.join(bases),
                                                                 qual='\t'.join(map(str, quantile_values))))

    for i, position in enumerate(positions):
        out.write("{pos}\t{dep:.1E}\t{qual}\t{base}\n".format(pos=position,
                                                    dep=depths[i],
                                                    base=nbasecalls[i],
                                                    qual='\t'.join(map(str, quantiles[i]))))

    out.close()

    if not args.nofigures:
        fig_kw = {'figsize':(8,6)}
        qualplot(positions, quantiles, args.output, fig_kw)
        qualdist(cycle_qual.values(), args.output, fig_kw)
        qualmap(cycle_qual, args.output, fig_kw)
        depthplot(positions, depths, args.output, fig_kw)
        gcplot(positions, cycle_nuc.values(), args.output, fig_kw)
        gcdist(cycle_gc, args.output, fig_kw)
        nucplot(positions, bases, cycle_nuc, args.output, fig_kw)
        if args.mbias:
            mbiasplot(positions, cycle_conv, args.output, fig_kw)
    time_finish = time.time()
    elapsed = time_finish - time_start
    if not args.quiet:
        sys.stderr.write("There were approximately {counts:,} reads in the file. Analysis finished in {sec}.\n".format(counts=act_nlines,
                                                                                                                  sec=time.strftime('%H:%M:%S',
                                                                                                                      time.gmtime(elapsed))
                                                                                                                      ))

def main():
    parser = argparse.ArgumentParser(prog='fastqp', description="simple NGS read quality assessment using Python")
    parser.add_argument('input', type=argparse.FileType('r'), help="input file (one of .sam, .bam, or .fastq(.gz) )")
    parser.add_argument('-q', '--quiet', action="store_true", default=False, help="do not print any messages (default: %(default)s)")
    parser.add_argument('-s', '--binsize', type=int, help='number of reads to bin for sampling (default: auto)')
    parser.add_argument('-n', '--nreads', type=int, default=2000000, help='number of reads sample from input (default: %(default)s)')
    parser.add_argument('-k', '--kmer', type=int, default=5, choices=range(2, 11), help='length of kmer for over-repesented kmer counts (default: %(default)s)')
    parser.add_argument('-o', '--output', type=str, default='plot', help="base name for output files (default: plot)")
    align_group = parser.add_mutually_exclusive_group()
    align_group.add_argument('--aligned', action="store_true", default=False, help="only aligned reads (default: %(default)s)")
    align_group.add_argument('--unaligned', action="store_true", default=False, help="only unaligned reads (default: %(default)s)")
    parser.add_argument('--nofigures', action="store_true", default=False, help="don't produce figures (default: %(default)s)")
    parser.add_argument('--nokmer', action="store_true", default=False, help="do not count kmers (default: %(default)s)")
    parser.add_argument('--mbias', action="store_true", default=False, help="make mbias plot for GEMINI reads (default: %(default)s)")

    args = parser.parse_args()
    run(args)

if __name__ == "__main__":
    main()
