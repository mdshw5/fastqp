fastqp
======
[![Build Status](https://travis-ci.org/mdshw5/fastqp.svg?)](https://travis-ci.org/mdshw5/fastqp)
[![PyPI](https://img.shields.io/pypi/v/fastqp.svg?)](https://pypi.python.org/pypi/fastqp)

Simple FASTQ and SAM read quality assessment and plotting using Python.

Features
--------

- Requires only Python with Numpy and Matplotlib libraries
- Works with (gzipped) FASTQ, SAM, and BAM formatted reads
- Tabular output statistics so you can create your own graphs
- A useful set of default graphics rivaling comparable QC packages
- Counts *all* IPUAC ambiguous nucleotide codes (NMWSKRY) if present in sequences
- Downsamples input files to around 2,000,000 reads (user adjustable)
- Allows a 5' and 3' (left and right) cycle limit for graphics generation

Requirements
------------

Tested on Python 2.6, 2.7, 3.2, 3.3, 3.4

Tested on Mac OS 10.9 and Linux 2.6.18

Installation
------------

    pip install fastqp

Note: BAM file support requires [samtools](http://samtools.sourceforge.net)

Usage
-----

    simple NGS read quality assessment using Python

    positional arguments:
      input                 input file (one of .sam, .bam, or .fastq(.gz) or stdin
                            (-))

    optional arguments:
      -h, --help            show this help message and exit
      -q, --quiet           do not print any messages (default: False)
      -s BINSIZE, --binsize BINSIZE
                            number of reads to bin for sampling (default: auto)
      -n NREADS, --nreads NREADS
                            number of reads sample from input (default: 2000000)
      -k {2,3,4,5,6,7,8,9,10}, --kmer {2,3,4,5,6,7,8,9,10}
                            length of kmer for over-repesented kmer counts
                            (default: 5)
      -o OUTPUT, --output OUTPUT
                            base name for output files (default: plot)
      -ll LEFTLIMIT, --leftlimit LEFTLIMIT
                            leftmost cycle limit (default: 1)
      -rl RIGHTLIMIT, --rightlimit RIGHTLIMIT
                            rightmost cycle limit (-1 for none) (default: -1)
      --aligned             only aligned reads (default: False)
      --unaligned           only unaligned reads (default: False)
      --nofigures           don't produce figures (default: False)
      --nokmer              do not count kmers (default: False)
      --gemstone            reads have convolution string (default: False)

Changes
-------

New in 0.1.5:

- Added `.fq` as acceptable file extension. (Thanks @danielecook)
- Added cycle-specific kmer plots

Examples
--------

![quality percentiles](https://raw.github.com/mdshw5/fastqp/master/examples/example_quals.png)

![gc plot](https://raw.github.com/mdshw5/fastqp/master/examples/example_gc.png)

![gc distribution](https://raw.github.com/mdshw5/fastqp/master/examples/example_gcdist.png)

![nucleotide plot](https://raw.github.com/mdshw5/fastqp/master/examples/example_nucs.png)

![depth plot](https://raw.github.com/mdshw5/fastqp/master/examples/example_depth.png)

![quality heatmap](https://raw.github.com/mdshw5/fastqp/master/examples/example_qualmap.png)

![quality distribution](https://raw.github.com/mdshw5/fastqp/master/examples/example_qualdist.png)

![kmer distribution](https://raw.github.com/mdshw5/fastqp/master/examples/example_kmers.png)


Acknowledgements
----------------
This project is freely licensed by the author, [Matthew Shirley](http://mattshirley.com), and was completed under the mentorship
and financial support of Drs. [Sarah Wheelan](http://sjwheelan.som.jhmi.edu) and [Vasan Yegnasubramanian](http://yegnalab.onc.jhmi.edu) at
the Sidney Kimmel Comprehensive Cancer Center in the Department of Oncology.
