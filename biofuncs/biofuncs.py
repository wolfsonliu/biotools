#! /bin/env python3
import os, sys, re
import numpy as np
import pandas as pd
import collections as cl
from Bio import SeqIO

# ------------------------------------------------------------------------------
# Functions

# ------------------
# gfftable
def gfftable(filename, fileformat, columns):
    '''Parse gff3 or gtf file into pandas DataFrame.'''
    # {{{

    if fileformat == 'gtf': # make different separation according to fileformat.
        sep1, sep2 = ' "', '";'
    elif fileformat == 'gff':
        sep1, sep2 = '=', ';'
    data = pd.read_table(
        filename,
        sep     = '\t',
        header  = None,
        names   = columns,
        index_col = False,
        comment = '#'
    )                           # read gtf/gff data
    attr_split = data[columns[-1]].apply(
        lambda row: [x.split(sep1)
                     for x in row.split(sep2) if len(x) > 2]
    )                      # split the attributes in attribute columns
    attr_dict = attr_split.apply(
        lambda row: dict(
            zip(
                [x[0].strip() for x in row],
                [x[1].strip() for x in row]
            )
        )
    )                           # make splited attributes into dicts
    attr_columns = attr_dict.apply(
        lambda row: list(row.keys())
    ).tolist()                  # get attr columns names
    attr_names = list(
        set([
            attr_columns[i][j] for i in range(len(attr_columns))
            for j in range(len(attr_columns[i]))
        ])
    )                           # get attr columns names
    attr = pd.DataFrame(
        dict(
            zip(
                attr_names,
                [
                    pd.Series(
                        [x[attr_name] if attr_name in x else np.NaN
                         for x in attr_dict]
                    ) for attr_name in attr_names
                ]
            )
        )
    )                           # make attr columns
    data = data.join(
        attr
    )                           # link attr columns with annotation.
    return data

    # }}}

# ------------------
# count_sgrna    
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
