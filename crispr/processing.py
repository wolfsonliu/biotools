#! /bin/env python3

# used to cound sgrna in pair end sequencing.

import os, sys, re
from biofuncs import count_sgrna

base = './F'

filepath = {
    'one1'   : os.path.join(base, 'PP10_WWS_CZZ-1_lane1_AGTCAA_L001_R1_001.fastq'),
    'one2'   : os.path.join(base, 'PP10_WWS_CZZ-1_lane2_AGTCAA_L002_R1_001.fastq'),
    'two1'   : os.path.join(base, 'PP10_WWS_CZZ-2_lane1_AGTTCC_L001_R1_001.fastq'),
    'two2'   : os.path.join(base, 'PP10_WWS_CZZ-2_lane2_AGTTCC_L002_R1_001.fastq'),
    'three1' : os.path.join(base, 'PP10_WWS_CZZ-3_lane1_ATGTCA_L001_R1_001.fastq'),
    'three2' : os.path.join(base, 'PP10_WWS_CZZ-3_lane2_ATGTCA_L002_R1_001.fastq'),
    'four1'  : os.path.join(base, 'PP10_WWS_CZZ-4_lane1_CCGTCC_L001_R1_001.fastq'),
    'four2'  : os.path.join(base, 'PP10_WWS_CZZ-4_lane2_CCGTCC_L002_R1_001.fastq')
}

lib = os.path.join(base, 'library.csv')

count = dict()
for condition in filepath.keys():
    data = count_sgrna(
        pattern = '^([ATGC]{19,23})GTTTA',
        read    = filepath[condition],
        library = lib
    )
    print('_'.join([os.path.basename(filepath[condition]).split('[.]')[0], 'result.csv']))
    data.to_csv(
        os.path.join(
            base,
            '_'.join([os.path.basename(filepath[condition]).split('[.]')[0], 'result.csv'])
        ),
        index = False
    )
