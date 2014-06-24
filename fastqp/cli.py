#!/usr/bin/python

import os
import sys
import argparse
import itertools
import math
import time
from fastqp import *


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

    stats = Stats()
    percent_complete = 10
    reads = infile.subsample(n)

    for read in reads:
        if ext in ['.sam', '.bam']:
            if args.aligned and not read.mapped:
                continue
            if args.unaligned and read.mapped:
                continue
            if read.reverse:
                if args.mbias:
                    conv = read.conv[::-1]
                seq = read.seq[::-1]
                qual = read.qual[::-1]
            else:
                if args.mbias:
                    conv = read.conv
                seq = read.seq
                qual = read.qual
        else:
            if args.mbias:
                conv = read.conv
            seq = read.seq
            qual = read.qual
        if args.mbias:
            stats.evaluate(seq, qual, conv)
        else:
            stats.evaluate(seq, qual)

        if not args.nokmer:
            stats.kmercount(seq, k=args.kmer)

        if not args.quiet and ext != '.gz':
            if (act_nlines / est_nlines) * 100 >= percent_complete:
                sys.stderr.write("Approximately {0:n}% complete at "
                                 "read {1:,} in {2}\n".format(percent_complete,
                                                              act_nlines,
                                                              time.strftime('%H:%M:%S',
                                                                            time.gmtime(time.time()-time_start))))
                percent_complete += 10
        act_nlines += n

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
                                                                     base='\t'.join(stats.nuc[1].keys()),
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
                                                                         base='\t'.join(stats.nuc[1].keys()),
                                                                         qual='\t'.join(map(str,quantile_values))))

            for i, position in enumerate(positions):
                out.write("{pos}\t{dep:n}\t{qual}\t{base}\n".format(pos=position,
                                                            dep=depths[i],
                                                            base=nbasecalls[i],
                                                            qual='\t'.join(map(str, quantiles[i]))))

            out.write('kmer\tcount\n')
            for key, value in stats.kmers.most_common(10):
                out.write(key + '\t' + str(value) + '\n')

    if not args.nofigures:
        fig_kw = {'figsize':(8,6)}
        basename = args.output if args.output else 'plot'
        qualplot(positions, quantiles, basename, fig_kw)
        qualdist(stats.qual.values(), basename, fig_kw)
        qualmap(stats.qual, basename, fig_kw)
        depthplot(positions, depths, basename, fig_kw)
        gcplot(positions, stats.nuc.values(), basename, fig_kw)
        gcdist(stats.gc, basename, fig_kw)
        nucplot(positions, bases, stats.nuc, basename, fig_kw)
        if args.mbias:
            mbiasplot(positions, stats.conv, basename, fig_kw)
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
    parser.add_argument('-o', '--output', type=str, help="base name for output files (default: plot)")
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
