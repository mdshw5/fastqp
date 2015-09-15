#!/usr/bin/python
# pylint: disable=bad-builtin

import os
import sys
import argparse
import itertools
import math
import time
import shlex
from fastqp import FastqReader, padbases, percentile, mean, bam_read_count, gc, window
from fastqp.plots import qualplot, qualdist, qualmap, nucplot, depthplot, gcplot, gcdist, kmerplot, mismatchplot, adaptermerplot
from fastqp.adapters import all_adapter_sequences
from collections import defaultdict
from simplesam import Reader, Sam
from subprocess import Popen, PIPE
from scipy import stats
from operator import mul
from six.moves import reduce


def run(args):
    """ read FASTQ or SAM and tabulate basic metrics """
    time_start = time.time()
    if args.input.name != '<stdin>':
        bsize = os.path.getsize(args.input.name)

    est_counter = int()
    sample_lengths = list()
    sample_binsizes = list()
    act_nlines = int()
    name, ext = os.path.splitext(args.input.name)
    if (args.leftlimit > 0) and (args.rightlimit > 0):
        if args.rightlimit < args.leftlimit:
            sys.exit("Left limit must be less than right limit.\n")
    if args.type:
        ext = '.' + args.type
    if ext not in ['.fq','.fastq', '.sam', '.bam', '.gz'] and args.input.name != '<stdin>':
        sys.exit("Input file must end in either .sam, .bam, .fastq, or .fastq.gz\n")

    if args.name:
        sample_name = args.name
    else:
        sample_name = args.input.name
        
    # estimate the number of lines in args.input if we can
    if ext in ['.fastq','.fq']:
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
            est_nlines = None
            if not args.quiet:
                sys.stderr.write("Reading from gzipped file, bin size (-s) set to {binsize:n}.\n".format(binsize=n))
        else:
            sys.stderr.write("Gzipped file detected. Reading file to determine bin size (-s).\n")
            p1 = Popen(shlex.split('gzip -dc %s' % args.input.name), stdout=PIPE)
            p2 = Popen(shlex.split('wc -l'), stdin=p1.stdout, stdout=PIPE)
            est_nlines, _ = p2.communicate()
            est_nlines = int(est_nlines) // 4
            if not args.quiet:
                sys.stderr.write("{est:,} reads in input file.\n".format(est=est_nlines))
    elif name == '<stdin>':
        if args.binsize:
            n = args.binsize
        else:
            n = 1
        if not args.quiet:
            sys.stderr.write("Reading from <stdin>, bin size (-s) set to {binsize:n}.\n".format(binsize=n))
        est_nlines = None
    if est_nlines is not None:
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
        infile = FastqReader(args.input, ext=ext)

    read_len = defaultdict(int)
    cycle_nuc = defaultdict(lambda: defaultdict(int))
    cycle_qual = defaultdict(lambda: defaultdict(int))
    cycle_gc = defaultdict(int)
    cycle_kmers = defaultdict(lambda: defaultdict(int))
    cycle_mismatch = {'C': defaultdict(lambda: defaultdict(int)),
                      'G': defaultdict(lambda: defaultdict(int)),
                      'A': defaultdict(lambda: defaultdict(int)),
                      'T': defaultdict(lambda: defaultdict(int))}
    if args.count_duplicates:
        bloom_filter = ScalableBloomFilter(mode=ScalableBloomFilter.SMALL_SET_GROWTH)
    duplicates = 0
    percent_complete = 10
    reads = infile.subsample(n)

    if args.count_duplicates:
        from pybloom import ScalableBloomFilter

    for read in reads:
        if isinstance(read, Sam):
            if args.aligned_only and not read.mapped:
                continue
            elif args.unaligned_only and read.mapped:
                continue
            if read.reverse:
                seq = read.seq[::-1]
                qual = read.qual[::-1]
            else:
                seq = read.seq
                qual = read.qual
        else:
            seq = read.seq
            qual = read.qual

        # Set up limits
        if (args.leftlimit == 1) and (args.rightlimit < 0):
            pass
        elif (args.leftlimit >= 1) and (args.rightlimit > 0):
            try:
                seq = seq[args.leftlimit - 1:args.rightlimit]
                qual = qual[args.leftlimit - 1:args.rightlimit]
            except IndexError:
                act_nlines += n
                continue
        elif (args.leftlimit > 1) and (args.rightlimit < 0):
            try:
                seq = seq[args.leftlimit - 1:]
                qual = qual[args.leftlimit - 1:]
            except IndexError:
                act_nlines += n
                continue
        if len(seq) == 0:
            act_nlines += n
            continue
        cycle_gc[gc(seq)] += 1

        if args.count_duplicates:
            if seq in bloom_filter:
                duplicates += 1
            else:
                bloom_filter.add(seq)

        for i, (s, q) in enumerate(zip(seq, qual)):
            cycle_nuc[args.leftlimit + i][s] += 1
            cycle_qual[args.leftlimit + i][q] += 1
        read_len[len(qual)] += 1

        for i, kmer in enumerate(window(seq, n=args.kmer)):
            cycle_kmers[args.leftlimit+i][kmer] += 1

        if isinstance(read, Sam) and read.mapped:
            try:
                ref = read.parse_md()
                for i, (s, r) in enumerate(zip(seq, ref)):
                    if s != r:
                        try:
                            cycle_mismatch[r][args.leftlimit+i][s] += 1
                        except KeyError:
                            pass
            except KeyError:
                pass


        if est_nlines is not None:
            if (act_nlines / est_nlines) * 100 >= percent_complete:
                sys.stderr.write("Approximately {0:n}% complete at "
                                 "read {1:,} in {2}\n".format(percent_complete,
                                                              act_nlines,
                                                              time.strftime('%H:%M:%S',
                                                                            time.gmtime(time.time()-time_start))))
                percent_complete += 10
        act_nlines += n

    positions = [k for k in sorted(cycle_qual.keys())]
    depths = [read_len[k] for k in sorted(read_len.keys())]

    basecalls = [cycle_nuc[k].keys() for k in sorted(cycle_nuc.keys())]
    bases = set(list(itertools.chain.from_iterable(basecalls)))
    #nbasecalls = [ '\t'.join([str(cycle_nuc[p].get(k, 0)) for k in bases]) for p in sorted(cycle_nuc.keys())]
    map(padbases(bases), cycle_nuc.values())

    quantile_values = [0.05,0.25,0.5,0.75,0.95]
    quantiles = []
    ## replace ASCII quality with integer
    for _, v in sorted(cycle_qual.items()):
        for q in tuple(v.keys()): ## py3 keys are iterator, so build a tuple to avoid recursion
            v[ord(str(q)) - 33] = v.pop(q)
        line = [percentile(v, p) for p in quantile_values]
        quantiles.append(line)

    # build kmer set of known adapter sequences
    adapter_kmers = set()
    for adapter in all_adapter_sequences:
        for kmer in window(adapter, n=args.kmer):
            adapter_kmers.add(kmer)

    # test for nonuniform kmer profiles and calculate obs/exp
    observed_expected = dict()
    all_kmers = [cycle_kmers[k].keys() for k in sorted(cycle_kmers.keys())]
    kmers = set(list(itertools.chain.from_iterable(all_kmers)))
    bad_kmers = []
    sequenced_bases = sum((l * n for l, n in read_len.items()))
    priors = tuple(map(float, args.base_probs.split(',')))
    for kmer in kmers:
        kmer_counts = [(i, cycle_kmers[i][kmer]) for i in sorted(cycle_kmers.keys())]
        expected_fraction = reduce(mul, (p ** kmer.count(b) for b, p in zip(('A', 'T', 'C', 'G', 'N'), priors)), 1)
        expected = expected_fraction * sequenced_bases
        observed_expected[kmer] = sum((n for _, n in kmer_counts)) / expected
        slope, _, _, p_value, _ = stats.linregress(*zip(*kmer_counts))
        if abs(slope) > 2 and p_value < 0.05:
            bad_kmers.append((kmer, slope, p_value))
    bad_kmers = sorted(bad_kmers, key=lambda x: x[2])[:10]
    

    # see http://vita.had.co.nz/papers/tidy-data.pdf
    sys.stdout.write("{row}\t{column}\t{value}\n".format(row=sample_name,
                                                         column='reads',
                                                         value=act_nlines))
   
    for cycle, count in read_len.items():
        sys.stdout.write("{row}\t{column}\t{value:n}\n".format(row=sample_name,
                                                               column='cycle_' + str(cycle),
                                                               value=count))

    for i, position in enumerate(positions):
        sys.stdout.write("{row}\t{column}\t{value:n}\n".format(row=sample_name,
                                                               column=str(position) + '_q05',
                                                               value=quantiles[i][0]))
        sys.stdout.write("{row}\t{column}\t{value:n}\n".format(row=sample_name,
                                                               column=str(position) + '_q25',
                                                               value=quantiles[i][1]))
        sys.stdout.write("{row}\t{column}\t{value:n}\n".format(row=sample_name,
                                                               column=str(position) + '_q50',
                                                               value=quantiles[i][2]))
        sys.stdout.write("{row}\t{column}\t{value:n}\n".format(row=sample_name,
                                                               column=str(position) + '_q75',
                                                               value=quantiles[i][3]))
        sys.stdout.write("{row}\t{column}\t{value:n}\n".format(row=sample_name,
                                                               column=str(position) + '_q95',
                                                               value=quantiles[i][4]))
    for base in bases:
        for position in positions:
            sys.stdout.write("{row}\t{column}\t{value:n}\n".format(row=sample_name,
                                                                   column=base + str(position),
                                                                   value=cycle_nuc[position][base]))
    for i in range(101):
        sys.stdout.write("{row}\t{column}\t{value:n}\n".format(row=sample_name,
                                                               column='read_gc_' + str(i),
                                                               value=cycle_gc[i]))    

    for kmer, obs_exp in sorted(observed_expected.items(), key=lambda x: x[1]):
        sys.stdout.write("{row}\t{column}\t{value:n}\n".format(row=sample_name,
                                                               column=kmer,
                                                               value=obs_exp))    

    if args.count_duplicates:
        sys.stdout.write("0\tpct_duplicate\t{0:.2%}\n".format(duplicates/act_nlines))

    from zipfile import ZipFile
    with ZipFile(args.output + '.zip', mode='w') as zip_archive:
        fig_kw = {'figsize':(8, 6)}
        qualplot(positions, quantiles, zip_archive, fig_kw)
        median_qual = qualdist(cycle_qual.values(), zip_archive, fig_kw)
        qualmap(cycle_qual, zip_archive, fig_kw)
        depthplot(read_len, zip_archive, fig_kw)
        gcplot(positions, cycle_nuc, zip_archive, fig_kw)
        gcdist(cycle_gc, zip_archive, fig_kw)
        nucplot(positions, bases, cycle_nuc, zip_archive, fig_kw)
        kmerplot(positions, cycle_kmers, zip_archive, [fields[0] for fields in bad_kmers], fig_kw)
        adaptermerplot(positions, cycle_kmers, adapter_kmers, zip_archive, fig_kw)
        if isinstance(infile, Reader):
            mismatchplot(positions
                         , cycle_mismatch, zip_archive, fig_kw)
    time_finish = time.time()
    elapsed = time_finish - time_start
    if not args.quiet:
        sys.stderr.write("There were {counts:,} reads in the file. Analysis finished in {sec}.\n".format(counts=act_nlines,
                                                                                                                       sec=time.strftime('%H:%M:%S',
                                                                                                                                         time.gmtime(elapsed))
        ))
        if len(bad_kmers) > 0:
            for kmer in bad_kmers:
                sys.stderr.write("KmerWarning: kmer %s has a non-uniform profile (slope = %s, p = %s).\n" % (kmer))
        if median_qual < args.median_qual:
            sys.stderr.write("QualityWarning: median base quality score is %s.\n" % median_qual)


def main():
    parser = argparse.ArgumentParser(prog='fastqp', description="simple NGS read quality assessment using Python")
    parser.add_argument('input', type=argparse.FileType('r'), help="input file (one of .sam, .bam, .fq, or .fastq(.gz) or stdin (-))")
    parser.add_argument('-q', '--quiet', action="store_true", default=False, help="do not print any messages (default: %(default)s)")
    parser.add_argument('-s', '--binsize', type=int, help='number of reads to bin for sampling (default: auto)')
    parser.add_argument('-a', '--name', type=str, help='sample name identifier for text and graphics output (default: input file name)')
    parser.add_argument('-n', '--nreads', type=int, default=2000000, help='number of reads sample from input (default: %(default)s)')
    parser.add_argument('-p', '--base-probs', type=str, default='0.25,0.25,0.25,0.25,0.1', help='probabilites for observing A,T,C,G,N in reads (default: %(default)s)')    
    parser.add_argument('-k', '--kmer', type=int, default=5, choices=range(2, 8), help='length of kmer for over-repesented kmer counts (default: %(default)s)')
    parser.add_argument('-o', '--output', type=str, default='fastqp_figures', help="base name for output files (default: %(default)s)")
    parser.add_argument('-t', '--type', type=str, default=None, choices=['fastq', 'gz', 'sam', 'bam'], help="file type (default: auto)")
    parser.add_argument('-ll', '--leftlimit', type=int, default=1, help="leftmost cycle limit (default: %(default)s)")
    parser.add_argument('-rl', '--rightlimit', type=int, default=-1, help="rightmost cycle limit (-1 for none) (default: %(default)s)")
    parser.add_argument('-mq', '--median-qual', type=int, default=30, help="median quality threshold for failing QC (default: %(default)s)")
    align_group = parser.add_mutually_exclusive_group()
    align_group.add_argument('--aligned-only', action="store_true", default=False, help="only aligned reads (default: %(default)s)")
    align_group.add_argument('--unaligned-only', action="store_true", default=False, help="only unaligned reads (default: %(default)s)")
    parser.add_argument('-d', '--count-duplicates', action="store_true", default=False, help="calculate sequence duplication rate (default: %(default)s)")

    args = parser.parse_args()
    run(args)

if __name__ == "__main__":
    main()
