#!/usr/bin/env Rscript
#-*-coding:utf-8-*-

##### Main file to display chromosome information

source("./chromosome_display.R")

##### File names

hubfile <- "hub.csv"
hotfile <- 't.csv'
satellitefile <- "satellite.tsv"
rRNAfile <- "rRNA.tsv"
tRNAfile <- "tRNA.tsv"
genomefile <- "ge.csv"

##### Read file

hub <- read.csv(
    hubfile,
    header=FALSE,
    stringsAsFactors=FALSE
)
colnames(hub) <- c(
    'chromosome', 'start', 'end', 'V4', 'V5', 'count',
    'V8', 'V9', 'V10', 'V11', 'V12', 'V13'
)

hot <- read.csv(
    hotfile,
    header=TRUE,
    stringsAsFactors=FALSE
)

tophub <- hub[sample(seq(1139), 20), ]

tophub$end <- tophub$end + 1000000



nemethaCGH <- rbind(
    nemethnad_r[nemethnad_r$detection_method == "aCGH",],
    nemethnad_r[nemethnad_r$detection_method == "aCGH/sequencing",]
)
nemethseq <- rbind(
    nemethnad_r[nemethnad_r$detection_method == "sequencing",],
    nemethnad_r[nemethnad_r$detection_method == "aCGH/sequencing",]
)
satellite_r <- read.csv(
    satellitefile,
    header = TRUE,
    sep    = '\t'
)
rRNA_r <- read.csv(
    rRNAfile,
    header = TRUE,
    sep    = '\t'
)
tRNA_r <- read.csv(
    tRNAfile,
    header = TRUE,
    sep    = '\t'
)
genome_r <- read.csv(
    genomefile,
    header = TRUE
)


##### Format genome information

genome <- data.frame(
    chromosome=c(paste0('chr', seq(1, 22)), 'chrX', 'chrY'),
    length=c(
        249250621, 243199373, 198022430, 191154276, 180915260,
        171115067, 159138663, 146364022, 141213431, 135534747,
        135006516, 133851895, 115169878, 107349540, 102531392,
        90354753, 81195210, 78077248, 59128983, 63025520,
        48129895, 51304566, 155270560, 59373566
    ),
    centromerestart=c(
        121535434, 92326171, 90504854, 49660117, 46405641,
        58830166, 58054331, 43838887, 47367679, 39254935,
        51644205, 34856694, 16000000, 16000000, 17000000,
        35335801, 22263006, 15460898, 24681782, 26369569,
        11288129, 13000000, 58632012, 10104553
    ),
    centromereend=c(
        124535434, 95326171, 93504854, 52660117, 49405641,
        61830166, 61054331, 46838887, 50367679, 42254935,
        54644205, 37856694, 19000000, 19000000, 20000000,
        38335801, 25263006, 18460898, 27681782, 29369569,
        14288129, 16000000, 61632012, 13104553
    )
)

##### Format data to be paint

paintdata <- list(hub, tophub)    #, rRNA_r)
paintdataname <- c("Hub", 'selected Hub')    #, "rDNA")
# [name][chromosome][start][end]
paintdatacol <- rbind(
    c(12, 1, 2, 3),
    c(12, 1, 2, 3)
)    #,c(12, 6, 7, 8))
data <- domainDataList(
    data=paintdata,
    dataname=paintdataname,
    datacol=paintdatacol
)
##### Paint the chromosomes

chromosomePlotRect(
    data=data,
    picname="chrbar.jpg",
    reference=genome,
    reverse=TRUE,
    chrcol='grey65',
    cencol='grey80',
    datacol=c('blue', 'red')
)


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
