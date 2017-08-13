#! /bin/bash
################################################################################
## Author:      Wolfson                                                       ##
## Date:        2016-06-12                                                    ##
## Modified:    2016-06-13                                                    ##
## Description: functions defined in RNA seq                                  ##
################################################################################

function createdirs
{
    ## Check dir existance and create dir if not exists.
    for file_path in $@; do
        if [ ! -e "${file_path}" ]; then
        mkdir -p ${file_path}
        fi
    done
}

function ngs_mapping_tophat
{
    ## Run tophat and cuffquant
    ########## tophat
    # $1 is the pair-end 1
    # $2 is the pair-end 2
    # $3 is the tophat_result_dir for sample
    # $4 reference index dir
    # $5 reference name
    # $6 reference gtf
    # $7 parallel
    ########## cuffquant
    # $8 is the cuffquant_result_dir for sample
    # $9 reference fa
    trap "" HUP
    local now_date=`date +%y%m%d`
    local output_bam_file=${3}/accepted_hits.bam
    tophat -p ${7:-1} -G ${6} -o ${3} ${4}/${5} ${1} ${2}
    wait
    cuffquant -p ${7:-1} -o ${8} -u -b ${9} ${6} ${output_bam_file}
}

