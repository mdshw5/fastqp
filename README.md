fastqp
======

simple NGS read quality assessment using Python

Requirements
------------
Python 2.7
    - matplotlib
    
Installation
------------
`python setup.py install`
    
Usage
-----

    usage: fastqp [-h] [-v] [-s SAMPLE] [-k {2,3,4,5,6,7,8,9,10}] [-o OUTPUT] [-f]
                  [--nokmer]
                  input
    
    simple NGS read quality assessment using Python
    
    positional arguments:
      input                 input file (FASTQ or SAM)
    
    optional arguments:
      -h, --help            show this help message and exit
      -v, --verbose         verbose output
      -s SAMPLE, --sample SAMPLE
                            number of reads to sample from
      -k {2,3,4,5,6,7,8,9,10}, --kmer {2,3,4,5,6,7,8,9,10}
                            length of kmer for over-repesented kmer counts
      -o OUTPUT, --output OUTPUT
                            base name for output files
      -f, --figures         produce figures
      --nokmer              do not count kmers