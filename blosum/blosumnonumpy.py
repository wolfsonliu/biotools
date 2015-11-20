#! /usr/bin/python3

#####Calculate the BLOSUM matrics without numpy
#####Author: Wolfson
#####Date: 20151017


import math

#####array and matrix functions
def py_array(m, i, j):
    return [[m for col in range(j)] for row in range(i)]

#####print results function
def print_sample_sequences(sequences):
    for seq in sequences:
        print(seq)

def print_matrix(mat, aa):
    print(' \t', end = '')
    for i in range(len(aa)):
        print(aa[i], end = '\t')
    print()
    for i in range(len(aa)):
        print(aa[i], end = '\t')
        for j in range(len(mat[i])):
            print('{0:.2}'.format(mat[i][j]), end = '\t')
        print()



#####blosum function

def blosum(seqs, aa, prt = True):
    ###get basic information
    aalen = len(aa)  #amino acids number
    seqlen = len(seqs[1])  #get the length of seq
    seqnum = len(seqs)  #get the number of sequences
    ###calculate 
    n = seqnum  #number of sequences
    w = seqlen  #length of sequences
    t = w*n*(n-1)/2  #total frequency
    ###initiate matrices
    c = py_array(0.0, aalen, aalen)  #initiate the cij frequency matrix
    q = py_array(0.0, aalen, aalen)  #initiate the qij frequence matrix
    e = py_array(0.0, aalen, aalen)  #initiate the eij matrix
    s = py_array(0.0, aalen, aalen)  #initiate the sij matrix
    blosum = py_array(' ', aalen, aalen)  #initiate the blosum matrix
    ###calculate matrices
    for i in range(aalen):
        for j in range(aalen):
            cij = [0] * seqlen
            for k in list(range(seqlen)):
                ni = 0.0
                nj = 0.0
                for m in list(range(seqnum)):
                    if aa[i] == seqs[m][k]:
                        ni += 1.0
                    if aa[j] == seqs[m][k]:
                        nj += 1.0
            #print(ni, nj)
                if i == j:
                    cij[k] = ni * (ni - 1) / 2.0
                elif i != j:
                    cij[k] = ni * nj
            #print(cij[k])
            c[i][j] = sum(cij)
            q[i][j] = c[i][j] / float(t)
    #print_matrix(c,aa)
    #print_matrix(q,aa)

    p1 = [sum([c[i][j] for i in range(aalen)]) for j in range(aalen)]
    #print(p1)
    p = [x / (float(t) * 2) for x in [p1[i] + c[i][i] for i in range(aalen)]]  #calculate the pi array
    #print(p)
    ###calculate blosum
    for i in range(aalen):
        for j in range(aalen):
            if i == j : e[i][j] = p[i] * p[i]
            elif i != j : e[i][j] = 2 * p[i] * p[j]
            if q[i][j] == 0: continue
            s[i][j] = math.log(q[i][j] / e[i][j], 2)
            blosum[i][j] = str(round(s[i][j] * 2))
    ###print matrices
    if prt == True:
        print("The c matrix is:")
        print_matrix(c, aa)
        print("The q matrix is:")
        print_matrix(q, aa)
        print("The e matrix is:")
        print_matrix(e, aa)
        print("The s matrix is:")
        print_matrix(s, aa)
        print("The BLOSUM is:")
        print_matrix(blosum, aa)
    return blosum


#####EOF
        

