#!/usr/bin/env Rscript
#-*-coding:utf-8-*-

##### Main file to display chromosome information
##### Author: Wolfson
##### Date: 20151029

source("./chromosome_display.R")
setwd("/lustre/user/liclab/lisky/liuzh/scripts/nucleolus/display_chromosome")

##### File names

nemethnadfile <- "nemethnad.csv"
satellitefile <- "satellite.tsv"
rRNAfile <- "rRNA.tsv"
tRNAfile <- "tRNA.tsv"
genomefile <- "ge.csv"

##### Read file

nemethnad_r <- read.csv(nemethnadfile, 
                        header = TRUE)
nemethaCGH <- rbind(nemethnad_r[nemethnad_r$detection_method == "aCGH",],
                    nemethnad_r[nemethnad_r$detection_method == "aCGH/sequencing",])
nemethseq <- rbind(nemethnad_r[nemethnad_r$detection_method == "sequencing",],
                    nemethnad_r[nemethnad_r$detection_method == "aCGH/sequencing",])
satellite_r <- read.csv(satellitefile, 
                        header = TRUE, 
                        sep    = '\t')
rRNA_r <- read.csv(rRNAfile, 
                   header = TRUE, 
                   sep    = '\t')
tRNA_r <- read.csv(tRNAfile, 
                   header = TRUE, 
                   sep    = '\t')
genome_r <- read.csv(genomefile, 
                     header = TRUE)


binsize <- c("1000", 
             "10000", 
             "100000", 
             "250000", 
             "50000", 
             "500000")
bin <- list()    # reference genome without mask.
for (m in binsize){
  bin[[m]] <- read.csv(paste(m, ".csv", sep = ''), header = FALSE)
  names(bin[[m]]) <- c("chromosome", 
                       "start", 
                       "end", 
                       "count")
  levels(bin[[m]]$chromosome) <- paste("chr",levels(bin[[m]]$chromosome), sep = "")
}
binm <- list()
for (m in binsize){
  binm[[m]] <- read.csv(paste("m", m, ".csv", sep = ''), header = FALSE)
  names(binm[[m]]) <- c("chromosome", 
                        "start", 
                        "end", 
                        "count")
  levels(binm[[m]]$chromosome) <- paste("chr",levels(binm[[m]]$chromosome), sep = "")
}

bincom <- list()    # reference genome with repeats mask.
for (m in binsize) {
  bincom[[m]]$toplevel <- read.csv(paste(m, ".csv", sep = ''), 
                                   header = FALSE)
  names(bincom[[m]]$toplevel) <- c("chromosome", 
                                   "start", 
                                   "end", 
                                   "count")
  bincom[[m]]$mask <- read.csv(paste("m",m, ".csv", sep = ''), 
                               header = FALSE)
  names(bincom[[m]]$mask) <- c("chromosome", 
                               "start", 
                               "end", 
                               "count")
}


##### Format genome information

genome <- genomeInformation(data         = genome_r,
                             cnchr        = 1,
                             cncentrstart = 2,
                             cncentrend   = 3,
                             cnlength     = 4)

##### Format data to be paint 

paintdata <- list(nemethnad_r, satellite_r)    #, rRNA_r)
paintdataname <- c("nemethnad", "satellite")    #, "rDNA")
# [name][chromosome][start][end]
paintdatacol <- rbind(c(5, 1, 2, 3),
                      c(12, 6, 7, 8))    #,c(12, 6, 7, 8))
data <- domainDataList(data     = paintdata,
                       dataname = paintdataname,
                       datacol  = paintdatacol)
### Generate compair data
paintdata1 <- list(nemethaCGH, satellite_r)    #, rRNA_r)
paintdataname1 <- c("nemethaCGH", "satellite")    #, "rDNA")
# [name][chromosome][start][end]
paintdatacol <- rbind(c(5, 1, 2, 3),
                      c(12, 6, 7, 8))    #, c(12, 6, 7, 8))
data1 <- domainDataList(data     = paintdata1,
                        dataname = paintdataname1,
                        datacol  = paintdatacol)

paintdata2 <- list(nemethseq, satellite_r)    #, rRNA_r)
paintdataname2 <- c("nemethseq", "satellite")    #, "rDNA")
# [name][chromosome][start][end]
paintdatacol <- rbind(c(5, 1, 2, 3),
                      c(12, 6, 7, 8))    #, c(12, 6, 7, 8))
data2 <- domainDataList(data     = paintdata2,
                        dataname = paintdataname2,
                        datacol  = paintdatacol)
#datac <- list()
#datac$aCGH <- data1
#datac$seq <- data2

##### Paint the chromosomes

chromosomePlotRect(data      = data,
                   picname   = "chrbar.jpg",
                   reference = genome,
                   reverse   = TRUE)


##### Paint the reads

for (n in binsize) {
  chromosomePlotBin(binread   = bin[[n]],
                 data      = data,
                 type      = "l",
                 lwd       = 1,
                 picname   = paste("binchrl", n, ".jpg", sep = ''),
                 reference = genome,
                 reverse   = TRUE)
  chromosomePlotBin(binread   = bin[[n]],
                 data      = data,
                 type      = "h",
                 lwd       = 1,
                 picname   = paste("binchrh", n, ".jpg", sep = ''),
                 reference = genome,
                 reverse   = TRUE)
  chromosomePlotBin(binread   = bin[[n]],
                 data      = data,
                 reference = genome,
                 picname   = paste("binchr", n, ".jpg", sep = ''),
                 single    = FALSE,
                 reverse   = TRUE)
  chromosomePlotLine(binread   = bin[[n]],
                       data      = data,
                       type      = "h",
                       reference = genome,
                       picname   = paste("chrsline", n, ".jpg", sep = ''))
}

for (m in binsize) {
  chromosomePlotBin(binread   = binm[[m]],
                 data      = data,
                 type      = "l",
                 lwd       = 1,
                 picname   = paste("mbinchrl", m, ".jpg", sep = ''),
                 reference = genome,
                 reverse   = TRUE)
  chromosomePlotBin(binread   = binm[[m]],
                 data      = data,
                 type      = "h",
                 lwd       = 1,
                 picname   = paste("mbinchrh", m, ".jpg", sep = ''),
                 reference = genome,
                 reverse   = TRUE)
  chromosomePlotBin(binread   = binm[[m]],
                 data      = data,
                 reference = genome,
                 picname   = paste("mbinchr", m, ".jpg", sep = ''),
                 single    = FALSE,
                 reverse   = TRUE)
  chromosomePlotLine(binread   = binm[[m]],
                       data      = data,
                       type      = "h",
                       reference = genome,
                       picname   = paste("mchrsline", m, ".jpg", sep = ''))
}

for (n in binsize) {
  compair2Mathods(methodname = c("aCGH", "sequencing"),
                  binread1   = bin[[n]],
                  data1      = data1,
                  binread2   = binm[[n]],
                  data2      = data2,
                  reference  = genome,
                  type       = "l",
                  lwd        = 1,
                  picname    = paste(n,"compair2.jpg", sep = ""))
}


##### Get the peak regions
source("./chromosome_display.R")
#sink("debug")
for (n in binsize) {
  cat("start", n, "\n")
  peaks <- chromosomeCallPeak(binm[[n]])
  write.table(x    = peaks,
              file = paste(n,"peakregion.csv", sep = ""),
              sep  = ",",
              row.names = FALSE)    # do not export row number.
  cat("end", n, "\n")
}
#sink()


#for (n in binsize) {
#  chrname <- chromosomesName(reverse = FALSE)
#  for (chri in chrname) {
#    peaks <- chromosomeDataPeak(binm[[n]][binm[[n]]$chromosome == chri,])
#    write.table(x    = peaks,
#                file = paste(n,"peak",chri,".csv", sep = ""),
#                sep  = ",",
#                row.names = FALSE)    # do not export row number.
#  }
#}

##### Try

source("./chromosome_display.R")

compair2Mathods(methodname = c("aCGH", "454"),
                binread1   = bin[[4]],
                data1      = data1,
                binread2   = binm[[4]],
                data2      = data2,
                reference  = genome,
                type       = "l",
                lwd        = 1,
                picname    = "compair2.jpg",
                single     = TRUE)

peaks <- chromosomeCallPeak(bin[[3]])


chromosomePlotBin(binread    = binm[[4]],
                   data      = data,
                   type      = "l",
                   lwd       = 1,
                   reference = genome,
                   picname   = paste("mbinchr", binsize[4], ".jpg", sep = ''),
                   single    = FALSE,
                   reverse   = TRUE)


compairMathod(binread   = bincom[[4]],
               data      = datac,
               reference = genome,
               picname   = "compair.jpg")



binReadDensity(binread = bin[[3]])

source("./chromosome_display.R")
nad250k <- read.csv("CopyOf250000peakregion.csv", 
                        header = TRUE)
paintdata <- list(nad250k, satellite_r)    #, rRNA_r)
paintdataname <- c("nad", "satellite")    #, "rDNA")
paintdatacol <- rbind(c(4, 1, 2, 3),
                      c(12, 6, 7, 8))    #,c(12, 6, 7, 8))
data <- domainDataList(data     = paintdata,
                       dataname = paintdataname,
                       datacol  = paintdatacol)
data <- list(nad250k)
names(data) <- "nad"
chromosomePlotRect(data      = data,
                   picname   = "rectfor250k.jpg",
                   reference = genome,
                   reverse   = TRUE)

paintdata <- list(nemethnad_r)    #, rRNA_r)
paintdataname <- "Nemeth"    #, "rDNA")
paintdatacol <- rbind(c(4, 1, 2, 3),
                      c(12, 6, 7, 8))    #,c(12, 6, 7, 8))
data <- domainDataList(data     = paintdata,
                       dataname = paintdataname,
                       datacol  = paintdatacol)
chromosomePlotRect(data      = data,
                   picname   = "Nemeth.jpg",
                   reference = genome,
                   reverse   = TRUE)