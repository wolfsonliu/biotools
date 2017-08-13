#! /usr/bin/python3
#-*-coding:utf-8-*-
#####main file to calculate the distribution

import argparse
import time
import os
import bincount

workdir = os.popen("pwd").read()[:-1] + '/'
#####argument settings
parser = argparse.ArgumentParser(prog = "countbin", #using %(prog)s to get the program name later.
                                 description = "Count numbers of reads in bins"
#parser.add_argument("-v", "--verbosity", type=int, choices=[0, 1, 2],default=0
#                    help="increase output verbosity")
#action="store_true" "count" 
#group = parser.add_mutually_exclusive_group()
#group.add_argument("-v", "--verbose", action="store_true")
#group.add_argument("-q", "--quiet", action="store_true")
parser.add_argument("--situation",
                    type = str, 
                    default = "countbin.out",
                    help = "run situation file, default to be countbin.out.")
parser.add_argument("-S", "--sam",
                    type = str,
                    help = "setting input sam file.")
parser.add_argument("-G", "--genome",
                    type = str,
                    default = "genome.csv",
                    help = "setting input genome file, default is genome.csv")
parser.add_argument("-d", "--directory",
                    type = str,
                    default = workdir,
                    help = "setting work directory, default is current working directory.")
parser.add_argument("-b", "--bin",
                    type = str,
                    default = '1000',
                    help = "setting bin size, use ',' to seperate different bin size without space, default is 1000.")
args = parser.parse_args()

#####program initiation
situationfile = open(args.situation, 'w+')  #open situation file to output running situation

binsize = [int(i) for i in args.bin.split(',')]  #transform type to int

####setting work directory
if args.sam[0:1] == './':
    workdirectory = ''
elif args.sam[0] == '/':
    workdirectory = ''
else:
    workdirectory = args.directory

dist = {} #initiate dist dict

#####run count bin reads numbers
for bin in binsize:
    samfile = workdirectory + args.sam
    genomefile = workdirectory + args.genome
    outfile = args.directory + str(bin) + ".csv"
    dist[str(bin)] = bincount.Distribution('genome.csv',bin = bin)  #initiate genome setting
    dist[str(bin)].calculate(samfile)  #count bin reads
    time.sleep(60)  #wait
    dist[str(bin)].writefile(outfile)  #output result
    situationfile.write('{0} is over.\n'.format(bin))  #situation file output
    time.sleep(10)

situationfile.close()

#####EOF
