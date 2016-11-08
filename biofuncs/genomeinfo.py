#! /bin/env python3
import os, sys, re
import numpy as np
import pandas as pd
import sqlalchemy as sa
import itertools
import functools
import matplotlib.pyplot as plt
from Bio import SeqIO, Seq
# ------------------
from biofuncs import gfftable

# ------------------------------------------------------------------------------
# Class

class GenomeInfo():
    # {{{

    def __init__(self, species, annofile, genomefile):
        self.species = species
        colnames = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'] # column names of gff and gtf are the same.
        def get_file_format(filename):
            if filename.split('.')[1] == 'gtf':
                return 'gtf'
            elif filename.split('.')[1] == 'gff':
                return 'gff'
            else:
                raise ValueError('Annotation file format might be wrong.')
        self.anno = gff_table(
            annofile,
            get_file_format(annofile),
            colnames
        )                       # link attr columns with annotation.
        self.genome = list(
            SeqIO.parse(genomefile, 'fasta')
        )
    def __repr__(self):
        return '<GenomeInfo [{0:s}]>'.format(self.species)
    def get_sgrna(self):
        # return DataFrame contains possible sgRNAs.
        if not hasattr(self, 'sgrna'):
            ngg = re.compile(
                '([atgcATGC]{20})([atgcATGC](GG|gg|Gg|gG))'
            )
            ccn = re.compile(
                '((CC|cc|Cc|cC)[atgcATGC])([atgcATGC]{20})'
            )
            columns = ['seqname', 'start', 'cut', 'end', 'sgrna', 'pam']
            sgrna = list()
            for chromosome in self.genome:
                sglist = [
                    {
                        'seqname': chromosome.id,
                        'start': x.start(),
                        'cut': x.end() - 6,
                        'end': x.end() - 3,
                        'sgrna': x.group(1),
                        'pam': x.group(2)
                    }
                    for x in ngg.finditer(str(chromosome.seq))
                ]
                sglist.extend(
                    {
                        'seqname': chromosome.id,
                        'start': x.start() + 3,
                        'cut': x.start() + 6,
                        'end': x.end(),
                        'sgrna': Seq.reverse_complement(x.group(3)),
                        'pam': Seq.reverse_complement(x.group(1))
                    }
                    for x in ccn.finditer(str(chromosome.seq))
                )
                sgrna.append(
                    pd.DataFrame(
                        sglist,
                        columns = columns
                    )
                )
            self.sgrna = pd.concat(sgrna, axis = 0, ignore_index = True)
        return self.sgrna
