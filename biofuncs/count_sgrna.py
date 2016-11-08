#! /bin/env python3
import numpy as np
import pandas as pd
import collections as cl
import os, sys, re
from Bio import SeqIO

# ------------------------------------------------------------------------------
# Functions

def count_sgrna(pattern     = 'ACCG([ATGC]{19,23})GTTTA',
          read        = '',
          pair        = '',
          library     = '',
          readthrough = False):
    '''count is used for counting the sgRNA reads from fastq sequence data.'''
    # Checking file exists
    files = {'one': read}
    if pair != '':
        files['two'] = pair
    if sum([os.path.exists(files[fq]) for fq in files]) < len(files):
        raise ValueError(
            'Fastq file not exists. Pleae check the file.'
        )
    if library == '' or not os.path.exists(library):
        raise ValueError(
            'Library should be a tab or comma separate file with the first row as header: sgrna, gene, sequence.'
        )
    with open(library) as fi:
        firstrow = fi.readline()
        if firstrow.find('\t') != -1:
            libformat = 'tsv'
        elif firstrow.find(',') != -1:
            libformat = 'csv'
        else:
            raise ValueError(
                'Library should be a tab or comma separate file with the first row as header: sgrna, gene, sequence.'
            )
    # Read data
    with open(library) as fi:   # read library
        if libformat == 'tsv':
            lib = pd.read_table(
                fi,
                header  = 0,
                names = ['sgrna', 'gene', 'sequence']
            )
        elif libformat == 'csv':
            lib = pd.read_csv(
                fi,
                header  = 0,
                names = ['sgrna', 'gene', 'sequence']
            )
    lib.index = lib['sequence']
    data = dict()
    for name, path in files.items(): # read fastq
        with open(path) as fq:
            data[name] = pd.Series(
                str(x.seq) for x in SeqIO.parse(fq, 'fastq')
            )
    if pair == '':
        # find the sgrna in single end reads
        reads = data['one'].str.findall(pattern)
    else:
        # find the sgrna in pair end reads
        reads = pd.concat(
            [data['one'].str.findall(pattern),
             data['two'].str.findall(pattern)],
            axis = 1
        )
        reads = reads[
            np.logical_xor(
                reads[0].map(len) > 0,
                reads[1].map(len) > 0
            )
        ]                       # remove pair end read both having the sgrna
        reads = pd.concat(
            [reads[0], reads[1]],
            axis = 0
        )
    reads = reads[reads.map(len) > 0].map(lambda x: x[0])
    count = reads.value_counts() # count sgrna reads
    count.name = 'read_count'
    result = lib.join(count)    # map sgrna counts to sgrna library
    return result
    
    
    
                
    


