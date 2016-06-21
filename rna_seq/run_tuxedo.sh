#! /bin/bash
################################################################################
## Author:      Wolfson                                                       ##
## Date:        2016-06-12                                                    ##
## Modified:                                                                  ##
## Description: use tophat2 to do reads mapping for project yangyl            ##
################################################################################

trap "" HUP
# no hup

################################################################################
##                             Source files                                   ##
################################################################################
source /home/liuzh/lustre/Code/yangyl/shell/definedfunction.sh
source /home/liuzh/lustre/Code/yangyl/shell/definedvariable.sh


################################################################################
##                       Create dir if not exist                              ##
################################################################################
createdirs \
    ${tophat_result_dir} \
    ${cuffdiff1_result_dir} \
    ${cuffdiff2_result_dir} \
    ${cuffdiff3_result_dir} \
    ${cuffdiff4_result_dir} \
    ${cuffquant2_result_dir} \
    ${log_dir}


################################################################################
##                 Choose cell line according to host                         ##
################################################################################
if [ "$(hostname)" = "knight" ]; then
    cell_line="jeko"
elif [ "$(hostname)" = "hi-c" ]; then
    cell_line="raji"
fi


################################################################################
##                            Tophat mapping                                  ##
################################################################################
declear -a tophat_array 
# used to store tophat process id

i_tophat=0

prefixnames=$(\
    ls ${clean_data_dir} | \
    grep ${cell_line} | \
    sed -r s/_[[:alnum:][:punct:]]*// | \
    uniq)

for samplename in ${prefixnames}; do
    createdirs \
        ${tophat_result_dir}/${samplename} \
        ${cuffquant2_result_dir}/${samplename}
    ngs_mapping_tophat \
        ${clean_data_dir}/${samplename}_1.clean.fq \
        ${clean_data_dir}/${samplename}_2.clean.fq \
        ${tophat_result_dir}/${samplename} \
        ${Homo_sapiens_index_dir} \
        ${Homo_sapiens_index} \
        ${Homo_sapiens_gtf} \
        ${cpu_num} \
        ${cuffquant2_result_dir}/${samplename} \
        ${Homo_sapiens_fa} >\
            ${log_dir}/tophat_${samplename}_$(date +%y%m%d%H%M).log &
        tophat_array[${i_tophat}]=$!
        i_tophat=$(echo ${i_tophat}+1 | bc)    
done


for ii in $(seq ${i_tophat}); do
   wait ${tophat_array[$(echo ${ii}-1 | bc)]}
   # wait until all tophat running over
done


# ################################################################################
# ##              Cuffdiff for differential expressed genes                     ##
# ################################################################################
declear -a cuffdiff_array 
# used to store cuffdiff process id

i_cuffdiff=0

for timeone in 0 6 24 48; do
    for timetwo in 6 24 48; do
        if [ ${timeone} -lt ${timetwo} ]; then
            cuffdiff \
                -o ${cuffdiff1_result_dir}/${cell_line}_${timeone}_${timetwo} \
                -b ${Homo_sapiens_fa} \
                -p ${cpu_num} -L H${timeone},H${timetwo} -u \
                ${Homo_sapiens_gtf} \
                ${tophat_result_dir}/${cell_line}-${timeone}-1/accepted_hits.bam,${tophat_result_dir}/${cell_line}-${timeone}-2/accepted_hits.bam,${tophat_result_dir}/${cell_line}-${timeone}-3/accepted_hits.bam \
                ${tophat_result_dir}/${cell_line}-${timetwo}-1/accepted_hits.bam,${tophat_result_dir}/${cell_line}-${timetwo}-2/accepted_hits.bam,${tophat_result_dir}/${cell_line}-${timetwo}-3/accepted_hits.bam &
            # using bam
            cuffdiff_array[${i_cuffdiff}]=$!
            i_cuffdiff=$(echo ${i_cuffdiff}+1 | bc)    

            cuffdiff \
                -o ${cuffdiff2_result_dir}/${cell_line}_${timeone}_${timetwo} \
                -b ${Homo_sapiens_fa} \
                -p ${cpu_num} -L H${timeone},H${timetwo} -u \
                ${Homo_sapiens_gtf} \
                ${cuffquant2_result_dir}/${cell_line}-${timeone}-1/abundances.cxb,${cuffquant2_result_dir}/${cell_line}-${timeone}-2/abundances.cxb,${cuffquant2_result_dir}/${cell_line}-${timeone}-3/abundances.cxb \
                ${cuffquant2_result_dir}/${cell_line}-${timetwo}-1/abundances.cxb,${cuffquant2_result_dir}/${cell_line}-${timetwo}-2/abundances.cxb,${cuffquant2_result_dir}/${cell_line}-${timetwo}-3/abundances.cxb &
            # using cxb from cuffquant
            cuffdiff_array[${i_cuffdiff}]=$!
            i_cuffdiff=$(echo ${i_cuffdiff}+1 | bc)    
        fi
    done
done

cuffdiff \
    -o ${cuffdiff3_result_dir}/${cell_line} \
    -b ${Homo_sapiens_fa} \
    -p ${cpu_big_num} -T -L H0,H6,H24,H48 -u \
    ${Homo_sapiens_gtf} \
    ${tophat_result_dir}/${cell_line}-0-1/accepted_hits.bam,${tophat_result_dir}/${cell_line}-0-2/accepted_hits.bam,${tophat_result_dir}/${cell_line}-0-3/accepted_hits.bam \
    ${tophat_result_dir}/${cell_line}-6-1/accepted_hits.bam,${tophat_result_dir}/${cell_line}-6-2/accepted_hits.bam,${tophat_result_dir}/${cell_line}-6-3/accepted_hits.bam \
    ${tophat_result_dir}/${cell_line}-24-1/accepted_hits.bam,${tophat_result_dir}/${cell_line}-24-2/accepted_hits.bam,${tophat_result_dir}/${cell_line}-24-3/accepted_hits.bam \
    ${tophat_result_dir}/${cell_line}-48-1/accepted_hits.bam,${tophat_result_dir}/${cell_line}-48-2/accepted_hits.bam,${tophat_result_dir}/${cell_line}-48-3/accepted_hits.bam &
# time series using bam
cuffdiff_array[${i_cuffdiff}]=$!
i_cuffdiff=$(echo ${i_cuffdiff}+1 | bc)    


cuffdiff \
    -o ${cuffdiff4_result_dir}/${cell_line} \
    -b ${Homo_sapiens_fa} \
    -p ${cpu_big_num} -T -L H0,H6,H24,H48 -u \
    ${Homo_sapiens_gtf} \
    ${cuffquant2_result_dir}/${cell_line}-0-1/abundances.cxb,${cuffquant2_result_dir}/${cell_line}-0-2/abundances.cxb,${cuffquant2_result_dir}/${cell_line}-0-3/abundances.cxb \
    ${cuffquant2_result_dir}/${cell_line}-6-1/abundances.cxb,${cuffquant2_result_dir}/${cell_line}-6-2/abundances.cxb,${cuffquant2_result_dir}/${cell_line}-6-3/abundances.cxb \
    ${cuffquant2_result_dir}/${cell_line}-24-1/abundances.cxb,${cuffquant2_result_dir}/${cell_line}-24-2/abundances.cxb,${cuffquant2_result_dir}/${cell_line}-24-3/abundances.cxb \
    ${cuffquant2_result_dir}/${cell_line}-48-1/abundances.cxb,${cuffquant2_result_dir}/${cell_line}-48-2/abundances.cxb,${cuffquant2_result_dir}/${cell_line}-48-3/abundances.cxb &
# time series using cxb
cuffdiff_array[${i_cuffdiff}]=$!
i_cuffdiff=$(echo ${i_cuffdiff}+1 | bc)    


for ii in $(seq ${i_cuffdiff}); do
   wait ${cuffdiff_array[$(echo ${ii}-1 | bc)]}
   # wait until all tophat running over
done


################################################################################
##                               Cuffnorm                                     ##
################################################################################
for timeone in 0 6 24 48; do
    for timetwo in 6 24 48; do
        if [ ${timeone} -lt ${timetwo} ]; then
            createdirs ${cuffnorm_result_dir}/${cell_line}_${timeone}_${timetwo}
            cuffnorm \
                -o ${cuffnorm_result_dir}/${cell_line}_${timeone}_${timetwo} \
                -p ${cpu_num} -L H${timeone},H${timetwo} \
                ${Homo_sapiens_gtf} \
                ${cuffquant2_result_dir}/${cell_line}-${timeone}-1/abundances.cxb,${cuffquant2_result_dir}/${cell_line}-${timeone}-2/abundances.cxb,${cuffquant2_result_dir}/${cell_line}-${timeone}-3/abundances.cxb \
                ${cuffquant2_result_dir}/${cell_line}-${timetwo}-1/abundances.cxb,${cuffquant2_result_dir}/${cell_line}-${timetwo}-2/abundances.cxb,${cuffquant2_result_dir}/${cell_line}-${timetwo}-3/abundances.cxb &
        fi
    done
done

createdirs ${cuffnorm_result_dir}/${cell_line}

cuffnorm \
    -o ${cuffnorm_result_dir}/${cell_line} \
    -p ${cpu_big_num} -T -L H0,H6,H24,H48 \
    ${Homo_sapiens_gtf} \
    ${cuffquant2_result_dir}/${cell_line}-0-1/abundances.cxb,${cuffquant2_result_dir}/${cell_line}-0-2/abundances.cxb,${cuffquant2_result_dir}/${cell_line}-0-3/abundances.cxb \
    ${cuffquant2_result_dir}/${cell_line}-6-1/abundances.cxb,${cuffquant2_result_dir}/${cell_line}-6-2/abundances.cxb,${cuffquant2_result_dir}/${cell_line}-6-3/abundances.cxb \
    ${cuffquant2_result_dir}/${cell_line}-24-1/abundances.cxb,${cuffquant2_result_dir}/${cell_line}-24-2/abundances.cxb,${cuffquant2_result_dir}/${cell_line}-24-3/abundances.cxb \
    ${cuffquant2_result_dir}/${cell_line}-48-1/abundances.cxb,${cuffquant2_result_dir}/${cell_line}-48-2/abundances.cxb,${cuffquant2_result_dir}/${cell_line}-48-3/abundances.cxb &


################################################################################
##                                 EOF                                        ##
################################################################################



