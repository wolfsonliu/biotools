#!/usr/bin/env Rscript
#-*-coding:utf-8-*-

#####Main file to display chromosome information
#####Author: Wolfson
#####Date: 20151017

source("./chrdisplay.r")
setwd("/lustre/user/liclab/lisky/liuzh/scripts/nucleolus/display_chromosome")

#####file names

nemethnadfile <- "nemethnad.csv"
satellitefile <- "satellite.tsv"
rRNAfile <- "rRNA.tsv"
tRNAfile <- "tRNA.tsv"
genomefile <- "ge.csv"

#####read file

nemethnad_r <- read.csv(nemethnadfile, header = TRUE)
satellite_r <- read.csv(satellitefile, header = TRUE, sep = '\t')
rRNA_r <- read.csv(rRNAfile, header = TRUE, sep = '\t')
tRNA_r <- read.csv(tRNAfile, header = TRUE, sep = '\t')
genome_r <- read.csv(genomefile, header = TRUE)

bin <- list()
binsize <- c("1000", "10000", "100000", "250000", "50000", "500000")
for (m in binsize){
  bin[[m]] <- read.csv(paste(m, ".csv", sep = ''), header = FALSE)
  names(bin[[m]]) <- c("chromosome", "start", "end", "count")
}


#####format genome information

genome <- genome_infomation(data = genome_r,
                            cnchr = 1,
                            cncentrstart = 2,
                            cncentrend = 3,
                            cnlength= 4)

#####format data to be paint 

paintdata <- list(nemethnad_r, satellite_r, rRNA_r)
paintdataname <- c("nemethnad", "satellite", "rDNA")
#[name][chromosome][start][end]
paintdatacol <- rbind(c(5, 1, 2, 3),
                      c(12, 6, 7, 8),
                      c(12, 6, 7, 8))
data <- domain_data_list(data = paintdata,
                         dataname = paintdataname,
                         datacol = paintdatacol)

#####paint the chromosomes

chromosomes(data = data,
            reference = genome)


#####paint the reads

for (n in c(1:length(binsize))){
  binreads_chromosome(binread = bin[[n]],
                      data = data,
                      picname = paste("binchr", binsize[n], ".jpg", sep = ''),
                      reference = genome)
  
  binreads_chromosome(binread = bin[[n]],
                      data = data,
                      reference = genome,
                      picname = paste("binchr", binsize[n], ".jpg", sep = ''),
                      single = FALSE)
}

source("./chrdisplay.r")
binreads_chromosome(binread = bin[[4]],
                    data = data,
                    type = "l",
                    reference = genome,
                    picname = paste("binchr", binsize[4], ".jpg", sep = ''),
                    single = FALSE)

chrreads(binread = bin[[4]],
         data = data,
         type = "h",
         reference = genome)