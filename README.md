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
      
Examples
--------

![quality percentiles](https://raw.github.com/mdshw5/fastqp/master/examples/example_quals.png)

![gc plot](https://raw.github.com/mdshw5/fastqp/master/examples/example_gc.png)

![gc distribution](https://raw.github.com/mdshw5/fastqp/master/examples/example_gcdist.png)

![nucleotide plot](https://raw.github.com/mdshw5/fastqp/master/examples/example_nucs.png)

![depth plot](https://raw.github.com/mdshw5/fastqp/master/examples/example_depth.png)

![quality distribution](https://raw.github.com/mdshw5/fastqp/master/examples/example_qualdist.png)

### Tabular output

    Pos	Depth	0.05	0.25	0.5	0.75	0.95	A	C	T	G	N
    1	4.0E+04	30	31	33	34	34	11750	8400	10644	9194	12
    2	4.0E+04	30	31	34	34	34	12589	7843	11610	7958	0
    3	4.0E+04	30	31	34	34	34	11808	7788	12924	7480	0
    4	4.0E+04	33	35	37	37	37	11498	8445	12273	7784	0
    5	4.0E+04	32	35	37	37	37	11637	7696	11703	8890	74
    6	4.0E+04	32	35	37	37	37	11791	7954	11303	8952	0
    7	4.0E+04	33	35	37	37	37	12552	7821	11594	8033	0
    8	4.0E+04	33	35	37	37	37	12554	7829	12086	7531	0
    9	4.0E+04	36	38	39	39	39	11291	7902	12079	8728	0
    10	4.0E+04	36	38	39	39	39	12467	7831	11690	8012	0
    11	4.0E+04	36	38	39	39	39	11534	8008	11494	8964	0
    12	4.0E+04	36	38	39	39	39	11517	8748	11677	8058	0
    13	4.0E+04	36	38	39	39	39	12225	8182	11744	7849	0
    14	4.0E+04	32	38	40	41	41	11626	9020	11535	7819	0
    15	4.0E+04	36	39	40	41	41	12546	7842	11627	7985	0
    16	4.0E+04	36	39	40	41	41	11456	8867	11616	8061	0
    17	4.0E+04	32	39	40	41	41	11445	8029	11557	8969	0
    18	4.0E+04	32	38	40	41	41	11540	7941	12639	7880	0
    19	4.0E+04	36	39	40	41	41	11422	8838	11794	7946	0
    20	4.0E+04	32	39	40	41	41	11485	8037	12497	7964	17
    .........................
    kmer	count
    AAAAA	25875
    TTTTT	22832
    ATTTT	14088
    AAAAT	13886
    TATTT	12018
    AAATA	11903
    AGAAA	11830
    TTTTA	11826
    TTTCT	11761
    TAAAA	11567
    
Acknowledgements
----------------
This project is freely licensed by the author, [Matthew Shirley](http://mattshirley.com), and was completed under the mentorship 
and financial support of Drs. [Sarah Wheelan](http://sjwheelan.som.jhmi.edu) and [Vasan Yegnasubramanian](http://yegnalab.onc.jhmi.edu) at 
the Sidney Kimmel Comprehensive Cancer Center in the Department of Oncology.