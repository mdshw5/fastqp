fastqp
======

Simple FASTQ and SAM read quality assessment and plotting using Python.

Features
--------

- Requires only Python with Numpy and Matplotlib libraries
- Works with (gzipped) FASTQ, SAM, and BAM formatted reads
- Tabular output statistics so you can create your own graphs
- A useful set of default graphics rivaling comparable QC packages
- Counts *all* IPUAC ambiguous nucleotide codes (NMWSKRY) if present in sequences
- Downsamples input files to around 2,000,000 reads (user adjustable)

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
      input                 input file (one of .sam, .bam, or .fastq(.gz) )

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
      --nofigures           don't produce figures (default: False)
      --nokmer              do not count kmers (default: False)
      --mbias               make mbias plot for GEMINI reads (default: False)

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

    Pos	Depth	0.05	0.25	0.5	0.75	0.95	S	T	W	K	M	N	A	C	G
    1	1.7E+04	28	31	31	33	34	3534	34	11	1396	3576	3837	4264	28	37
    2	1.7E+04	28	31	31	33	34	4151	31	6	8	1114	4529	2869	3956	53
    3	1.7E+04	27	31	31	34	34	26	14	13	11	301	5063	3355	4031	3903
    4	1.7E+04	32	35	37	37	37	37	8	7	6	4572	3693	3946	4448
    5	1.7E+04	32	35	37	37	37	4271	7	10	9	4378	3883	4145	14
    6	1.7E+04	32	35	37	37	37	17	12	13	7	4401	4146	3996	4125
    7	1.7E+04	32	35	37	37	37	3990	9	9	2	4451	3859	4377	6	14
    8	1.7E+04	32	35	37	37	37	9	13	12	11	2	4472	4042	4331	3825
    9	1.7E+04	32	37	39	39	39	3744	12	9	2	4274	4214	4433	13	16
    10	1.7E+04	32	37	39	39	39	18	4	8	9	1	4350	4140	4339	3848
    11	1.7E+04	32	37	39	39	39	3898	8	11	3	1	4410	4022	4351	13
    12	1.7E+04	32	37	39	39	39	10	7	8	12	1	4195	4117	4463	3904
    13	1.7E+04	32	37	39	39	39	12	6	9	4309	4016	4293	8	4064
    14	1.7E+04	32	38	39	40	41	29	9	16	5	4404	3907	4381	3966
    15	1.7E+04	32	38	39	40	41	12	5	2	6	4600	3990	4271	3831
    16	1.7E+04	32	38	39	40	41	13	5	11	4474	3950	4335	14	3915
    17	1.7E+04	32	38	39	40	41	18	5	11	8	4392	4014	4266	4003
    18	1.7E+04	32	38	39	40	41	3938	6	11	4279	4077	4376	8	22
    19	1.7E+04	32	38	40	40	41	3879	3	11	7	4314	4160	4326	17
    20	1.7E+04	32	38	39	40	41	3898	5	12	4436	4111	4217	12	26
    21	1.7E+04	32	38	39	40	41	20	19	10	17	1	4286	4091	4371	3902
    22	1.7E+04	32	38	39	40	41	18	6	16	4360	4002	4269	6	4040
    23	1.7E+04	32	37	39	40	41	22	15	14	4258	4064	4362	13	3969
    24	1.7E+04	32	37	39	40	41	21	14	14	13	4406	3984	4266	3999
    25	1.7E+04	32	37	39	40	41	19	12	13	4270	4062	4244	17	4080
    26	1.7E+04	32	37	39	40	41	4100	19	25	14	4437	3851	4244	27
    27	1.7E+04	32	37	39	40	41	36	11	13	4508	3813	4291	13	4030
    28	1.7E+04	32	37	39	40	41	22	11	23	4234	4220	4227	9	3956
    29	1.7E+04	32	37	39	40	41	22	20	11	11	4289	3993	4265	4077
    30	1.7E+04	32	37	39	40	41	59	30	20	18	4321	4247	4200	3782
    31	1.7E+04	30	37	39	40	41	41	15	19	13	4312	4076	4090	4106
    32	1.7E+04	30	37	39	40	41	40	14	18	4395	4102	4206	14	3877
    33	1.7E+04	30	37	39	40	41	34	26	14	18	4274	4056	4286	3951
    34	1.7E+04	30	37	39	40	41	32	13	20	4286	4048	4069	16	4166
    35	1.7E+04	30	37	39	40	41	34	21	18	13	4448	3907	4229	3951
    36	1.7E+04	30	37	39	40	41	44	24	21	11	4298	3900	4230	4065
    37	1.7E+04	30	37	39	40	41	45	20	39	4	4225	3961	4121	12	4131
    38	1.7E+04	30	37	39	40	41	39	29	33	4375	3990	4147	28	3875
    39	1.6E+04	30	37	39	40	41	3938	22	30	4233	3904	4282	32	42
    40	1.6E+04	30	37	39	40	41	37	25	15	15	4375	3941	4074	3952
    41	1.6E+04	30	37	39	40	41	3968	20	12	13	4273	3830	4231	36
    42	1.6E+04	30	36	39	40	41	29	20	28	4181	3798	4364	25	3885
    43	1.6E+04	30	36	39	40	41	40	17	26	24	4045	4016	4090	41	3984
    44	1.6E+04	30	36	39	40	41	35	24	23	4255	3763	4205	22	3902
    45	1.6E+04	30	36	39	40	41	36	32	24	15	4149	3952	4129	3846
    46	1.6E+04	30	36	38	40	41	3872	28	22	24	4164	3896	4078	40
    47	1.6E+04	30	36	38	40	41	44	24	20	14	4416	3711	4026	3802
    48	1.6E+04	30	36	38	40	41	41	22	14	20	4205	3765	4112	3826
    49	1.6E+04	28	36	38	40	41	3984	21	28	3979	3800	4044	27	52
    50	1.6E+04	28	36	38	40	41	37	27	15	18	3991	3758	4093	3918
    51	1.6E+04	28	35	38	39	41	37	14	27	4095	3826	4025	18	3741
    52	1.6E+04	28	35	38	39	41	3743	21	22	4041	3824	3992	32	45
    53	1.6E+04	28	35	38	40	41	3830	18	19	4049	3709	3955	20	49
    54	1.6E+04	28	36	38	40	41	3702	16	20	3982	3772	4032	18	44
    55	1.5E+04	28	36	38	40	41	3716	13	28	3914	3763	3991	20	52
    56	1.5E+04	28	36	38	40	41	3626	21	16	4004	3676	4001	28	42
    57	1.5E+04	28	36	38	40	41	3706	14	22	3993	3573	3966	19	48
    58	1.5E+04	28	35	38	40	41	3629	34	23	16	4	4011	3526	3948	59
    59	1.5E+04	28	35	38	40	41	3751	27	26	16	3792	3602	3920	48
    60	1.5E+04	28	35	38	40	41	52	21	22	11	4023	3550	3860	3541
    61	1.5E+04	28	35	38	40	41	3537	33	17	13	1	3689	3743	3899	46
    62	1.5E+04	28	35	38	40	41	3568	28	19	9	3970	3438	3786	57
    63	1.5E+04	28	35	38	40	41	55	33	24	23	3775	3468	3977	3424
    64	1.5E+04	27	34	38	40	41	3593	34	19	22	3764	3482	3691	41
    65	1.5E+04	26	34	37	39	41	3501	21	24	3817	3423	3641	30	62
    66	1.4E+04	26	34	37	39	41	46	32	23	18	3668	3318	3802	3481
    67	1.4E+04	26	34	37	39	41	3351	34	24	21	3601	3530	3609	53
    68	1.4E+04	26	34	37	39	41	44	38	20	22	3593	3342	3555	3435
    69	1.4E+04	26	34	36	39	41	3314	34	23	12	3566	3304	3607	50
    70	1.4E+04	26	34	36	39	41	3227	29	27	21	3490	3328	3544	49
    71	1.4E+04	26	33	36	39	40	3291	24	34	3423	3204	3455	26	64
    72	1.3E+04	26	33	36	39	40	3101	12	25	1	3437	3128	3547	34	44
    73	1.3E+04	26	33	36	38	40	3164	35	24	24	3389	3160	3329	47
    74	1.3E+04	26	33	36	38	40	3090	14	41	3334	3062	3352	18	86
    75	1.3E+04	25	33	35	38	40	70	44	27	21	4	3275	2914	3436	3001
    76	1.3E+04	22	31	34	36	39	75	16	38	1	3272	2971	3218	27	2946
    77	1.2E+04	23	31	34	36	39	2909	14	34	3139	2871	3304	23	76
    78	1.2E+04	24	32	35	37	39	67	30	31	16	1	3065	2934	3202	2836
    79	1.2E+04	25	32	35	37	39	67	36	26	14	1	2956	2904	3243	2734
    80	1.2E+04	26	32	35	36	39	74	28	23	16	3008	2904	3070	2658
    81	1.2E+04	24	32	35	36	39	2548	23	36	18	1	3043	2803	3011	97
    82	1.1E+04	24	32	35	36	39	2501	26	35	2855	2794	3074	35	88
    83	1.1E+04	24	32	35	36	39	2505	46	26	19	2777	2794	2956	75
    84	1.1E+04	24	32	34	36	39	106	21	32	2592	2731	3010	43	2480
    85	1.1E+04	24	32	34	36	39	2473	25	31	2563	2698	2907	64	96
    86	1.1E+04	22	32	34	36	39	111	64	33	23	2408	2669	2934	2391
    87	1.0E+04	22	32	34	35	38	101	89	33	30	2350	2593	2933	2332
    88	1.0E+04	22	32	34	35	37	2313	126	29	29	2285	2434	2994	103
    89	1.0E+04	20	32	34	35	37	165	180	56	19	2315	2339	2839	2257
    90	1.0E+04	20	32	34	35	37	204	221	60	33	2214	2366	2755	2148
    91	9.9E+03	18	32	34	35	37	2029	241	85	23	2216	2339	2702	224
    92	9.7E+03	18	31	34	35	37	2012	294	135	25	2199	2275	2511	272
    93	9.6E+03	18	32	34	35	36	1914	315	157	28	3	2082	2112	2663	290
    94	9.4E+03	16	31	34	35	36	334	326	212	26	1954	2036	2504	2016
    95	9.1E+03	25	31	34	35	36	334	24	246	2117	1969	2408	35	1992
    96	8.9E+03	2	31	34	35	36	1853	29	277	17	2343	1887	2126	346
    97	8.6E+03	2	31	34	35	36	55	21	300	2241	2012	2071	21	1907
    98	8.4E+03	2	31	34	35	36	48	37	308	10	1	1940	1939	2204	1879
    99	8.0E+03	2	31	34	35	36	2070	1920	2101	1943
    100	7.7E+03	2	31	34	35	36	2107	2046	2014	1558
    kmer	count
    TGCTG	5641
    CTGCT	5623
    ACTGC	4980
    GCTGT	4819
    CTGTA	4700
    TGTAC	4638
    TACTG	3068
    CACTG	3046
    TGACC	2739
    CTTTA	2693

Acknowledgements
----------------
This project is freely licensed by the author, [Matthew Shirley](http://mattshirley.com), and was completed under the mentorship
and financial support of Drs. [Sarah Wheelan](http://sjwheelan.som.jhmi.edu) and [Vasan Yegnasubramanian](http://yegnalab.onc.jhmi.edu) at
the Sidney Kimmel Comprehensive Cancer Center in the Department of Oncology.
