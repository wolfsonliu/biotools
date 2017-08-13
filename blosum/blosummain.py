#! /usr/bin/python3
# -*-coding:utf-8-*-

import blosumnonumpy as blosum

#aa = np.array(['C', 'S', 'T', 'P', 'A', 'G', 'N', 'D', 'E', 'Q', 'H', 'R', 'K', 'M', 'I', 'L', 'V', 'F', 'Y', 'W'], str)  #amino acid 20
aashort = ['A', 'I', 'L', 'S', 'T', 'V']

seq1 = ['A', 'A', 'I']
seq2 = ['S', 'A', 'L']
seq3 = ['T', 'A', 'L']
seq4 = ['T', 'A', 'V']
seq5 = ['A', 'A', 'L']
sequences = [seq1, seq2, seq3, seq4, seq5]

blosum.blosum(sequences, aashort)
