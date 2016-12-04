#! /bin/env python2

################################################################################
## Author:   Wolfson                                                          ##
## Date:     20160217                                                         ##
## Modified: 20160218                                                         ##
## Description:                                                               ##
##    get the guide region for Crisper/Cas9 sgRNA                             ##
################################################################################

################################################################################
## Information                                                                ##
################################################################################
## cut: exon: -16 or 16 intron: 4 or -4                                       ##
## NGG: exon: -13 or 13 intron: 7 or -7                                       ##
## GG:  exon: -12 or 12 intron: 8 or -8                                       ##
## need the calculation of cutting site                                       ##
## final data:                                                                ##
##     * transcript (gene name)                                               ##
##     * element of cutting site                                              ##
##     * cutting site to the edge                                             ##
##     * 30 bp sequence (24 before PAM, 3 after PAM)                          ##
## algorithm:                                                                 ##
##     * get all the guide RNA.                                               ## 
##     * get the information of gene.                                         ##
##     * judge the place of PAM, save right ones.                             ##
## input:                                                                     ##
##     * gene fasta                                                           ##
## output:                                                                    ##
##     * tsv file with columns                                                ##
##         * id: gene id                                                      ##
##         * num: num of sequence for one gene                                ##
##         * pam: site of pam                                                 ##
##         * cut: site of cut                                                 ##
##         * distance: distance of cut site to the junction of exon and       ##
##               intron                                                       ##
##         * exon/intron: if cut inside exon, e; if cut inside intron, i      ##
##         * seq: sequence, with exon base in upper, intron base in lower     ##
##         * antisense: 1 for antisense                                       ##
################################################################################

from Bio import SeqIO
import re
import argparse

# arg pase
parser = argparse.ArgumentParser(description='Get the guide sequence for sgRNA.')

parser.add_argument('-i', help = "input file name, fasta",
                    dest = "infilename",
                    required = True)
parser.add_argument('-o', help = "output file name",
                    dest = "outfilename",
                    default = "./result.tsv")

args = parser.parse_args()



class SgRNA:
    def __init__(self, gene, num, pam, cut):
        self.id            = gene
        self.num           = num
        self.pam_site      = pam
        self.cut_site      = cut
        self.cut_d2j       = 0
        self.cut_place     = ""
        self.guide         = ""
        self.antisense     = ""
    
def isAccept(site, start, end):
    # judge site bigger than start and smaller than end.
    for n in range(len(start)):
        if site >= start[n] and site <= end[n]:
            return True
    return False

def distanceToJunction(site, start, end, junction):
    # return the distance of site to the nearest junction of exon and intron.
    for n in range(len(start)):
        if site >= start[n] and site <= end[n]:
            return site - junction[n]
        
records = SeqIO.index(args.infilename, "fasta")
# index fasta file


gene_names = list()

with open(args.infilename, "r") as fasta:
    for line in fasta:
        if line.startswith(">"):
            gene_names.append(line.replace(">","").split(" ")[0])

result = list()

for gene in gene_names:
    junction_pattern_ei = re.compile("[A-Z][a-z]")    # exon before intron
    junction_pattern_ie = re.compile("[a-z][A-Z]")    # intron before exon
    pam_pattern    = re.compile("[A-Za-z][Gg][Gg]")
    num                 = 0    # the number for the guide sequence of one gene
    for is_reverse in range(2):
        pam_start_site      = -1    # should exist in is_reverse for, or the later one would have problem.
        sequence = str(records[gene].seq)
        if is_reverse == 1:
            sequence = str(records[gene].reverse_complement().seq)
            # get the reversed sequence
        accept_start        = []    # accept region site floor for cutting site
        accept_end          = []    # accept region site roof for cutting site
        junction            = []
        tmpsite1            = -1
        tmpsite2            = -1
        for n in range(len(junction_pattern_ei.findall(sequence))):
        # store the accept region for cutting sites.
            tmpsite1 = junction_pattern_ei.search(sequence, tmpsite1 + 1).start()
            if tmpsite1 - 15 >= 0:
                junction.append(tmpsite1 + 1)
                accept_start.append(tmpsite1 - 15)
                accept_end.append(tmpsite1 + 4 + 1)
            tmpsite2 = junction_pattern_ie.search(sequence, tmpsite2 + 1).start()
            if tmpsite2 + 16 + 1 <= len(sequence) - 1:
                junction.append(tmpsite2 + 1)
                accept_start.append(tmpsite2 - 3)
                accept_end.append(tmpsite2 + 16 + 1)
        if len(accept_start) != len(accept_end) or len(accept_start) != len(junction) or len(accept_end) != len(junction):
            print "accept start end and junction not the same length"
            exit
        while True:
            
            if pam_pattern.search(sequence, pam_start_site + 1) is None:
            # end the while loop if at the end.
                break
            pam_start_site = pam_pattern.search(sequence, pam_start_site + 1).start()
            if pam_start_site - 24 < 0:
                continue
            cutting_site = pam_start_site - 3
            if isAccept(cutting_site, accept_start, accept_end):
            # cutting site inside acceptted region.
                result_record          = SgRNA(gene, num, pam_start_site, cutting_site)
                num                    += 1    # number increase
                result_record.cut_d2j  = distanceToJunction(cutting_site, accept_start, accept_end, junction)
                if sequence[result_record.cut_site] in "ATCG":
                    result_record.cut_place = "e"
                else:
                    result_record.cut_place = "i"
                result_record.guide = sequence[pam_start_site - 24 : pam_start_site + 6]
                result_record.antisense = str(is_reverse)
                result.append("\t".join((result_record.id, 
                                         str(result_record.num),
                                         str(result_record.pam_site),
                                         str(result_record.cut_site),
                                         str(result_record.cut_d2j),
                                         result_record.cut_place,
                                         result_record.guide,
                                         result_record.antisense)))
            

with open(args.outfilename, "w") as output:
    output.write("\t".join(("id", "num", "pam", "cut", "distance", "exon/intron", "seq", "antisense\n")))
    for line in result:
        output.write(line + "\n")



                
