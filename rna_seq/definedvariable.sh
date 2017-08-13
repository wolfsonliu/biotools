#! /bin/bash
################################################################################
## Author:      Wolfson                                                       ##
## Date:        2016-06-12                                                    ##
## Modified:    2016-06-13                                                    ##
## Description: variables defined in RNA seq                                  ##
################################################################################

# Reference
Homo_sapiens_fa=/home/liuzh/lustre/Code/yangyl/data/genome.fa
Homo_sapiens_index_dir=/home/liuzh/lustre/Code/yangyl/data
Homo_sapiens_index=genome
Homo_sapiens_gtf=/home/liuzh/lustre/Code/yangyl/data/gene.pc.gtf

# Data directory
data_dir=/home/liuzh/lustre/Data/yangyl
clean_data_dir=/home/liuzh/lustre/Data/yangyl/clean_data
tophat_result_dir=${data_dir}/tophat_result2
cuffdiff1_result_dir=${data_dir}/cuffdiff_result1 # using bam
cuffdiff2_result_dir=${data_dir}/cuffdiff_result2 # using cxb from cuffquant
cuffdiff3_result_dir=${data_dir}/cuffdiff_result3 # time series using bam
cuffdiff4_result_dir=${data_dir}/cuffdiff_result4 # time series using cxb
cuffquant2_result_dir=${data_dir}/cuffquant_result2 # used by cuffdiff
cuffnorm_result_dir=${data_dir}/cuffnorm_result2
log_dir=/home/liuzh/lustre/Code/yangyl/log

# Setting parallel
cpu_num=2
cpu_big_num=8