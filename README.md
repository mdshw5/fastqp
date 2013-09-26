fastqp
======

Simple FASTQ and SAM read quality assessment and plotting using Python. 

Requirements
------------

Python 2.7

    - numpy
    - matplotlib
    
Tested on Mac OS 10.8 and Linux 2.6.18
    
Installation
------------

    sudo pip install numpy matplotlib https://github.com/mdshw5/fastqp/archive/master.zip
    
Usage
-----

    usage: fastqp [-h] [-q] [-s SAMPLE] [-k {2,3,4,5,6,7,8,9,10}] [-o OUTPUT]
                  [-t {sam,fastq}] [-f] [--nokmer]
                  input
    
    simple NGS read quality assessment using Python
    
    positional arguments:
      input                 input file(.gz))
    
    optional arguments:
      -h, --help            show this help message and exit
      -q, --quiet           do not print any messages
      -s SAMPLE, --sample SAMPLE
                            number of reads to bin for sampling
      -k {2,3,4,5,6,7,8,9,10}, --kmer {2,3,4,5,6,7,8,9,10}
                            length of kmer for over-repesented kmer counts
      -o OUTPUT, --output OUTPUT
                            base name for output files
      -t {sam,fastq}, --type {sam,fastq}
                            file type for input file (default=fastq)
      -f, --figures         produce figures
      --nokmer              do not count kmers
    
    Note: fastqp randomly samples ~200,000 reads from the input file by default.
    To change the number of reads sample, specify the number of reads to bin for
    sampling. For example, '-s 100' will sample 1 in 100 reads. To evaluate the
    entire file set '-s 1'.
      
Examples
--------

![quality percentiles](https://raw.github.com/mdshw5/fastqp/master/examples/example_quals.png)

![gc plot](https://raw.github.com/mdshw5/fastqp/master/examples/example_gc.png)

![gc distribution](https://raw.github.com/mdshw5/fastqp/master/examples/example_gcdist.png)

![nucleotide plot](https://raw.github.com/mdshw5/fastqp/master/examples/example_nucs.png)

![depth plot](https://raw.github.com/mdshw5/fastqp/master/examples/example_depth.png)

![quality heatmap](https://raw.github.com/mdshw5/fastqp/master/examples/example_qualmap.png)

![quality distribution](https://raw.github.com/mdshw5/fastqp/master/examples/example_qualdist.png)

### Tabular output

    Pos	Depth	0.05	0.25	0.5	0.75	0.95	A	C	T	G	N
    1	2.3E+04	2	27	31	31	33	6087	4377	6552	3736	2228
    2	2.3E+04	2	25	31	31	34	7246	3254	6063	4569	1848
    3	2.3E+04	16	26	30	31	34	7782	3858	6084	4787	469
    4	2.3E+04	27	32	35	35	37	7534	4064	6128	5249	5
    5	2.3E+04	27	32	35	35	37	6957	4479	6696	4845	3
    6	2.3E+04	27	32	35	35	37	7162	4517	6510	4782	9
    7	2.3E+04	26	32	35	36	37	6969	4580	6747	4676	8
    8	2.3E+04	26	32	35	36	37	7021	4613	6726	4610	10
    9	2.3E+04	25	33	35	38	39	6955	4588	6955	4474	8
    10	2.3E+04	23	32	37	38	39	6978	4654	6801	4541	6
    11	2.3E+04	25	33	35	38	39	6968	4581	6872	4553	6
    12	2.3E+04	23	32	35	38	39	7052	4579	6817	4527	5
    13	2.3E+04	25	33	35	38	39	7083	4607	6915	4372	3
    14	2.3E+04	24	33	37	39	41	7021	4616	6776	4565	2
    15	2.3E+04	24	33	37	39	41	7092	4566	6840	4480	2
    16	2.3E+04	24	33	37	39	41	6929	4594	6864	4589	4
    17	2.3E+04	24	33	37	39	41	6942	4571	6860	4604	3
    18	2.3E+04	24	33	37	39	41	6999	4605	6836	4537	3
    19	2.3E+04	21	32	37	39	41	6952	4722	6810	4493	3
    20	2.3E+04	21	32	37	39	41	7000	4633	6849	4494	4
    21	2.3E+04	21	32	37	39	41	6909	4659	6793	4617	2
    22	2.3E+04	21	32	37	39	41	7031	4594	6786	4567	2
    23	2.3E+04	21	32	37	39	41	6968	4642	6787	4583	0
    24	2.3E+04	21	32	37	39	41	6922	4656	6779	4622	1
    25	2.3E+04	21	32	37	39	41	6948	4628	6690	4713	1
    26	2.3E+04	18	32	37	39	41	6942	4593	6759	4685	1
    27	2.3E+04	18	32	37	39	40	6900	4704	6791	4582	1
    28	2.3E+04	18	32	36	39	40	7005	4488	6822	4654	4
    29	2.3E+04	21	32	36	39	40	6891	4655	6751	4667	1
    30	2.3E+04	16	32	36	39	40	7062	4586	6770	4540	3
    .........................
    kmer	count
    GTACA	17874
    ATAAA	12880
    CTGCT	12772
    TGCTG	12768
    ACTGC	12648
    CTGTA	12207
    GCTGT	11781
    TGTAC	11521
    AAAGT	11159
    TAAAG	10748
    
Acknowledgements
----------------
This project is freely licensed by the author, [Matthew Shirley](http://mattshirley.com), and was completed under the mentorship 
and financial support of Drs. [Sarah Wheelan](http://sjwheelan.som.jhmi.edu) and [Vasan Yegnasubramanian](http://yegnalab.onc.jhmi.edu) at 
the Sidney Kimmel Comprehensive Cancer Center in the Department of Oncology.