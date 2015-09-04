"""
Oligonucleotide sequences copyright 2007-2013 Illumina, Inc. All rights reserved. 
Derivative works created by Illumina customers are authorized for use with 
Illumina instruments and products only. All other uses are strictly prohibited.
"""

from itertools import chain


class IndexedAdapter(object):
    def __init__(self, name, p1, p2=None, index_sequences=[None]):
        self.name = name
        self.p1 = p1
        self.p2 = p2
        self.index_sequences = index_sequences
        self.adapters = []
        if p2 is None:
            assert index_sequences == [None]
            self.adapters.append(p1)
        elif index_sequences == [None]:
            assert p1 is not None and p2 is not None
            self.adapters.append(''.join([p1, p2]))
        elif index_sequences is not None:
            assert p1 is not None and p2 is not None
            for index in index_sequences:
                self.adapters.append(''.join([p1, index, p2]))

    def __repr__(self):
        try:
            p1 = self.p1[:5]
        except (IndexError, TypeError):
            p1 = self.p1
        try:
            p2 = self.p2[-5:]
        except (IndexError, TypeError):
            p2 = self.p2
        try:
            index = '%s, %s' % (self.index_sequences[0], self.index_sequences[1])
        except (IndexError, TypeError):
            index = '%s' % self.index_sequences[0]
        if len(self.index_sequences) > 2:
            index += '...'
        return "%s: %s...%s (%s)" % (self.name, p1, p2, index)


# TODO: add small RNA prep kit, synthetic longread kit
    
truseq_ht_i5 = IndexedAdapter('TruSeq Stranded RNA/DNA HT i5',
                              'AATGATACGGCGACCACCGAGATCTACAC',
                              'ACACTCTTTCCCTACACGACGCTCTTCCGATCT',
                              ('TATAGCCT', 'ATAGAGGC', 'CCTATCCT', 'GGCTCTGA',
                               'AGGCGAAG', 'TAATCTTA', 'CAGGACGT', 'GTACTGAC',
                               'AGGCTATA', 'GCCTCTAT', 'AGGATAGG', 'TCAGAGCC',
                               'CTTCGCCT', 'TAAGATTA', 'ACGTCCTG', 'GTCAGTAC'))

truseq_ht_i7 = IndexedAdapter('TruSeq Stranded RNA/DNA HT i7',
                              'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC',
                              'ATCTCGTATGCCGTCTTCTGCTTG',
                              ('ATTACTCG', 'TCCGGAGA', 'CGCTCATT', 'GAGATTCC', 'ATTCAGAA', 'GAATTCGT',
                               'CTGAAGCT', 'TAATGCGC', 'CGGCTATG', 'TCCGCGAA', 'TCTCGCGC', 'AGCGATAG'))

truseq_universal = IndexedAdapter('TruSeq Universal Adapter',
                                  'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT')

truseq_lt = IndexedAdapter('TruSeq v1/v2/LT Adapters',
                           'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC',
                           'ATCTCGTATGCCGTCTTCTGCTTG',
                           ('ATCACG', 'CGATGT', 'TTAGGC', 'TGACCA', 'ACAGTG', 'GCCAAT',
                            'CAGATC', 'ACTTGA', 'GATCAG', 'TAGCTT', 'GGCTAC', 'CTTGTA',
                            'AGTCAA', 'AGTTCC', 'ATGTCA', 'CCGTCC', 'GTCCGC', 'GTGAAA',
                            'GTGGCC', 'GTTTCG', 'CGTACG', 'GAGTGG', 'ACTGAT', 'ATTCCT'))

nextera_dna_p1 = IndexedAdapter('Nextera DNA kit/1',
                                'AATGATACGGCGACCACCGAGATCTACACGCCTCCCTCGCGCCATCAG')
nextera_dna_p2 = IndexedAdapter('Nextera DNA kit/2',
                                'CAAGCAGAAGACGGCATACGAGATCGGTCTGCCTTGCCAGCCCGCTCAG')
nextera_tpn1 = IndexedAdapter('Nextera DNA kit tpn/1',
                                'GCCTCCCTCGCGCCATCAGAGATGTGTATAAGAGACAG')
nextera_tpn2 = IndexedAdapter('Nextera DNA kit tpn/2',
                                'GCCTTGCCAGCCCGCTCAGAGATGTGTATAAGAGACAG')


all_adapter_sequences = tuple(chain(*(getattr(kit, 'adapters') for kit in (truseq_lt, truseq_ht_i5, truseq_ht_i7,
                                                                          truseq_universal, nextera_dna_p1,
                                                                          nextera_dna_p2, nextera_tpn1,
                                                                          nextera_tpn2))))
