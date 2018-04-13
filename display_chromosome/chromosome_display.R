#!/usr/bin/env Rscript
#-*-coding:utf-8-*-

##### Functions to display chromosomes information
##### Author: Wolfson
##### Date: 20151029

####################
##### List
####################
##### Function-chromosomesName: return chromosome name vector: chr1.
##### Function-genomeInformation: change reference genome information to rigth format.
##### Function-domainDataList: format data to [name][chromosome][start][end].
##### Function-getFileType: get the type of input filename.
##### Function-chooseFunction: choose functions according to words.
##### Function-startPicture: change the picture output type.
##### Function-dataMerge: merge the start and end counts for paint.
##### Function-dataFrameValueNorm: normalize reads count number to c(bottom,top).   # need more features.
##### Function-chromosomeDataPeak: get the peak of reads count.    # threshold should be rethinked.
##### Function-chromosomeCallPeak: get the peak region of data.
##### Function-chromosomePlotRect: paint chromomes with different kinds of units with bar.
##### Function-chromosomePlotLine: paint chromomes with reads using lines.
##### Function-chromosomePlotBin: paint read numbers of bins and chromosomes.
##### Function-compairMathod: compair different mathods results.    # not been finished.
##### Function-compair2Mathods: compair different mathods results.    # need adjustment.
##### Function-binReadDensity: show the reads density.
####################


##### Function-chromosomesName: return chromosome name vector: chr1.

chromosomesName <- function(reverse=TRUE)
{
  chrname <- as.vector(c("chrY", "chrX", "chr22", "chr21", "chr20", "chr19", 
                         "chr18", "chr17", "chr16", "chr15", "chr14", "chr13", 
                         "chr12", "chr11", "chr10", "chr9", "chr8", "chr7", 
                         "chr6", "chr5", "chr4", "chr3", "chr2", "chr1"))
  if (reverse == FALSE) {
    chrname <- rev(chrname)
  }
  ### Remove useless data
  return(chrname)
}


##### Function-genomeInformation: change reference genome information to rigth format.

genomeInformation <- function(data,    # data frame of genome informations: columns(chr)(centrstart)(centrend)(length) 
                              cnchr,    # column number of chr
                              cncentrstart,    # column number of start of centromere
                              cncentrend,    # column number of end of centromere
                              cnlength)    # column number of length of chromosome
{
  tmp <- cbind(data[cnchr],
               data[cncentrstart],
               data[cncentrend],
               data[cnlength])
  colnames(tmp) <- c("chromosome",
                     "centromerestart",
                     "centromereend",
                     "length")
  return(tmp)
}

##### Function-domainDataList: format data to [name][chromosome][start][end].

domainDataList <- function(data,    # list of data file
                           dataname,    # vector of data names
                           datacol)    # table of column numbers of data list
{
  datalen <- length(data)
  if (length(dataname) != datalen) {
    print("Data and dataname are not the same length.",
          "Please check again!")
  }
  tmpdata <- list()
  for (n in 1:datalen) {
    tmp <- data[[n]][,datacol[n,]]
    names(tmp) <- c("name", 
                    "chromosome", 
                    "start", 
                    "end")
    tmpdata[[dataname[n]]] <- tmp
  }
  names(tmpdata) <- dataname
  ### Remove useless data
  rm(tmp)
  return(tmpdata)
}


##### Function-getFileType: get the type of input filename.

getFileType <- function(filename)
{
  filetype <- rev(strsplit(filename, split = '[.]')[[1]])[1]
  #or using grepl()
  return(filetype)
}


##### Function-chooseFunction: choose functions according to words.

chooseFunction <- function(filename,    # filename.
                           funcs,    # functions.
                           words)    # reference words vector.
{
  wordslen <- length(words)
  funcslen <- length(funcs)
  if (wordslen != funcslen) {
    print("Functions and words are not the same length.")
    return(0)
  }
  filetype <- getFileType(filename = filename)
  for (i in 1:funcslen) {
    if (filetype == words[i]) {
      func <- funcs[[i]]
      break
    }
  }
  return(func)
}


##### Function-startPicture: change the picture output type.

startPicture <- function(filename,
                         width=800,
                         height=600,
                         ...) {
                                        # this function is used to change different kinds of pictures.
                                        # remember to dev.off() after the picture.
    if (unlist(strsplit(filename, '[.]'))[-1] == 'pdf') {
        pdf(filename, width=width/72, height=height/72)
    } else {
        picfunc <- chooseFunction(
            filename=filename,
            funcs=c(jpeg, png, tiff, bmp),
            words=c("jpg", "png", "tiff", "bmp")
        )
        picfunc(
            filename=filename,    # file name, end with picture style.
            width=width,    # width.
            height=height,    # height.
            ...
        )
    }
}


##### Function-dataMerge: merge the start and end counts for paint.

dataMerge <- function(binread,
                      chr,    # chromosome
                      stc=c(1:2,4),    # order: chr start count.
                      edc=c(1,3:4),    # order: chr end count.
                      bstart=0,    # set the bin start from 0 or minit of bin.
                      normalize=TRUE)    # make the normalization of data.
{
  readsstart <- binread[binread$chromosome == chr, stc]    # get the bin start data of chromosome chri.
  readsend <- binread[binread$chromosome == chr, edc]    # get the bin end data of chromosome chri.
  names(readsstart) <- c("chromosome", 
                         "site", 
                         "count")
  names(readsend) <- c("chromosome", 
                       "site", 
                       "count")
  normthreshold <- summary(readsstart$count)[5] * 3   # 3 times of 3rd quartile to do the normalization.
  binsize <- min(readsend$site) - min(readsstart$site)    # bin size.
  binstart <- bstart    # the start of bin site
  binend <- max(readsend$site)    # the end of bin site.
  binnum <- (binend - binstart) / binsize    # number of bins.
  chrseq <- rep(chr, times = binnum)    # new chromosome.
  startseq <- seq(binstart, binend - binsize, by = binsize)    # full site of start.
  endseq <- seq(binstart + binsize, binend, by = binsize)    # full site of end.
  zeroseq <- rep(0, times = binnum)    # initiate count.
  bintablestart <- cbind(chrseq, as.data.frame(startseq), zeroseq)    # "chromosome", "site", "count". 
  # WARNING: using as.data.frame to keep the startseq as num seq
  bintableend <-  cbind(chrseq, as.data.frame(endseq), zeroseq)    # "chromosome", "site", "count".
  names(bintablestart) <- c("V1", "V2", "V3")
  names(bintableend) <- c("V1", "V2", "V3")
  for (n in readsstart$site) {
    # count number
    startcount <- readsstart[readsstart$site == n, 3]
    endcount <- readsend[readsstart$site == n, 3]
    if (readsstart[readsstart$site == n, 3] > normthreshold) {
      startcount <- normthreshold
    }
    if (readsend[readsstart$site == n, 3] > normthreshold) {
      endcount <- normthreshold
    }
    bintablestart[bintablestart$V2 == n, 3] <- startcount
    bintableend[bintableend$V2 ==  n, 3] <- endcount
  }
  tmp <- rbind(bintablestart, bintableend)
  names(tmp) <- c("chromosome",
                  "site",
                  "count")
  reads <- tmp[order(tmp$site),]    # make the data frame with a right order.
  ### Remove useless data, maybe no use.
  rm(chrseq)
  rm(startseq)
  rm(endseq)
  rm(tmp)
  rm(readsstart)
  rm(readsend)
  return(reads)
}

##### Function-dataFrameValueNorm: normalize reads count number to c(bottom,top).

dataFrameValueNorm <- function(data,    # data frame to be normalized.
                               column,    # column to be normalized.
                               range=c(0,1),    # bottom and top of min, default to be 0 and 1.
                               limit=FALSE,    # whether limit the max and min value.
                               bound=c(0,1))    # when limit is TRUE, the bound of data, bottom and top.
{
  if (limit == TRUE) {
    data[data[column] <= bound[1], column] <- bound[1]    # limit bottom values.
    data[data[column] >= bound[2], column] <- bound[2]    # limit top values.
  }
  top.value <- max(data[column])
  bottom.value <- min(data[column])
  data[column] <- data[column] / top.value
  return(data)
}

##### Function-chromosomeDataPeak: get the peak of reads count.

chromosomeDataPeak <- function(data,    # data should only include information of one chromosome.
                               nam=c("chromosome", "start", "end", "count"))    # order: chr start count.
{
  names(data) <- nam
  #data <- dataFrameValueNorm(data, 4)
  normthreshold <- summary(data$count)[5] * 0.7 + summary(data$count)[6] * 0.3    # 3rd quartile to do the normalization.
  tmp <- cbind(as.data.frame(character()),
               start=integer(),
               end=integer(),
               count=integer())
  # WARNING: using as.data.frame to keep the startseq as num seq.
  names(tmp) <- nam
  for (n in data$start) {
    # count number
    if (data[data$start == n, 4] > normthreshold) {
      tmp <- rbind(tmp, data[data$start == n,])
    }
  }
  reads <- tmp[order(tmp$start),]    # make the data frame with a right order.
  ### Remove useless data, maybe no use.
  return(reads[reads$count != 0,])
}


##### Function-chromosomeCallPeak: get the peak region of data.

chromosomeCallPeak <- function(data,    # data frame of bin reads, [chromosome, start, end, count].
                               nam=c("chromosome", "start", "end", "count"),    # order: chr start count.
                               chrname=0,    # the list of chromosomes to be paint with reversed order.
                               reverse=FALSE)
{
  if (chrname == 0) {
    chrname <- chromosomesName(reverse = reverse)    # chrname is the vector of chromosome names.
  }
  # initiate temporary data frames to store peaks information.
  # tmpdf1 is used to store peaks data.
  # tmpdf2 is used to store regions data.
  tmpdf2 <- cbind(as.data.frame(character()),
                  start=integer(),
                  end=integer(),
                  count=integer())
  names(tmpdf2) <- nam
  for (chri in chrname) {
    # generate the peaks.
    #cat("start", chri, "\n")
    tmpdf1 <- chromosomeDataPeak(data[data$chromosome == chri,])
    names(tmpdf1) <- nam
    if (nrow(tmpdf1[tmpdf1$chromosome == chri,]) == 0) {
      next
    } else if (nrow(tmpdf1[tmpdf1$chromosome == chri,]) == 1){
      tmpdf2 <- rbind(tmpdf2, tmpdf1[tmpdf1$chromosome == chri,][1,]) 
      next
    }
    tmpdf2 <- rbind(tmpdf2, tmpdf1[tmpdf1$chromosome == chri,][1,]) 
    #str(tmpdf1)
    binsize <- tmpdf1[1, "end"] - tmpdf1[1, "start"]
    sameregion <- 100000    # disdance to be considered in one NAD
    for (n in 2:nrow(tmpdf1)) {
      # generate the peak regions.
      tmpdf2row <- nrow(tmpdf2[tmpdf2$chromosome == chri,])
      tmpdf2rowend <- tmpdf2[tmpdf2$chromosome == chri,][tmpdf2row, "end"]
      ifbool <- FALSE
      #for (binsizenum in c(0:3)) {
      #  ifbool <- ifbool || tmpdf2rowend + binsize * binsizenum == tmpdf1[n, "start"]
      #}
      ifbool <- ifbool || (tmpdf2rowend <= tmpdf1[n, "start"] &&tmpdf2rowend + sameregion >= tmpdf1[n, "start"])
      if (ifbool) {
        tmpdf2[tmpdf2$chromosome == chri,][tmpdf2row, "count"] <- tmpdf2[tmpdf2row, "count"] + tmpdf1[n, "count"]
        tmpdf2[tmpdf2$chromosome == chri,][tmpdf2row, "end"] <- tmpdf1[n, "end"]
        next
      }
      tmpdf2 <- rbind(tmpdf2, tmpdf1[tmpdf1$chromosome == chri,][n,])
    }
    #cat("end", chri, "\n")
  }
  return(tmpdf2)
}
##### Function-chromosomePlotRect: paint chromomes with different kinds of units with bar.


chromosomePlotRect <- function(data,    # list of units to paint.
                               reference,    # reference genome informations, use genome_information() to load genome information first.
                               datacol=c("red", "royalblue4", "darkolivegreen", "blanchedalmond"),
                               chrname=NULL,    # the list of chromosomes to be paint with reversed order.
                               chrcol="deepskyblue",    # color of chromosomes.
                               cencol="darkorange",    # color of centromeres.
                               picname="chrbar.jpg",
                               w=1366,
                               h=768,
                               u="px",
                               reverse=FALSE)
{
  if (is.null(chrname)) {
    chrname <- chromosomesName(reverse = reverse)    # chrname is the vector of chromosome names.
  }
  startPicture(
      filename = picname,
      width    = w,
      height   = h,
      units    = u
  )    #need dev.off() later
  basepairlength <- max(reference[, "length"])    # length of chr1 the longest.
  rectheight <- basepairlength/96    # height of rectangules.
  #### Paint blank plot
  plot(
      0,
      0,
      type="n",
      xlim=c(0, basepairlength),
      ylim=c(0, rectheight * length(chrname) * 2),
      col="white",
      xlab="",
      ylab="",
      axes=FALSE
  )
  #par(new = TRUE)
  #### Paint chromosomes with information
  for (n in c(1:length(chrname))) {
    righty <- rectheight * (2 * n - 1)
    lefty <- rectheight * 2 * (n - 1)
    ### Paint chromosome
    ## Chromosome part before centromere
    rect(
        xleft=0,  #left bottom x
        ybottom=lefty + 0.25 * rectheight,    # left bottom y.
        xright=reference[reference$chromosome == chrname[n], "centromerestart"],    # right top x.
        ytop=righty - 0.25 * rectheight,    # right top y.
        col=chrcol,
        border=chrcol
    )
    ## Chromosome part after centromere
    rect(
        xleft=reference[reference$chromosome == chrname[n], "centromereend"],    # left bottom x.
        ybottom=lefty + 0.25 * rectheight,    # left bottom y
        xright=reference[reference$chromosome == chrname[n], "length"],    # right top x.
        ytop=righty - 0.25 * rectheight,    # right top y.
        col=chrcol,
        border=chrcol
    )

    ### Paint centromeres
    rect(
        xleft=reference[reference$chromosome == chrname[n], "centromerestart"],
        ybottom=lefty + 0.3 * rectheight,
        xright=reference[reference$chromosome == chrname[n], "centromereend"],
        ytop=righty - 0.3 * rectheight,
        col=cencol,
        border=cencol
    )
    ## chromosome round end
    points(
        x=c(
            0,
            reference[reference$chromosome == chrname[n], "centromerestart"],
            reference[reference$chromosome == chrname[n], "centromereend"],
            reference[reference$chromosome == chrname[n], "length"]
        ),
        y=c(
            lefty + 0.5 * rectheight,
            righty - 0.5 * rectheight,
            lefty + 0.5 * rectheight,
            righty - 0.5 * rectheight
        ),    # left bottom y
        col=chrcol,
        pch=19,
        cex=1
    )
    ### Paint data regions
    for (m in c(1:length(data))) {
        if (!chrname[n] %in% unique(data[[m]]$chromosome)) {
            next
        }
        rect(
            xleft=data[[m]][data[[m]]$chromosome == chrname[n], "start"],
            ybottom=lefty,
            xright=data[[m]][data[[m]]$chromosome == chrname[n], "end"],
            ytop=righty,
            col=datacol[m],
            border=datacol[m]
        )
    }
  }
  #### Add label for chromosomes
  text(
      x=-rectheight * 3,
      y=((1:length(chrname)) * 4 - 3) * 0.5 * rectheight,
      labels=substr(chrname, 4, 6),
      cex=1.1
  )
  #### Legend.
  legend(
      "right",
      legend=c(names(data)[1:length(data)], "Centromere"),
      pch=15,
      col=c(datacol[1:length(data)], cencol),
      border="white",
      cex=1.5
  )
  dev.off()
}



##### Function-chromosomePlotLine: paint chromomes with reads using lines.

chromosomePlotLine <- function(binread,    # bin reads data, [chromosome][start][end][count].
                               data,    # list of units to paint.
                               reference,    # reference genome informations, use genome_information() to load genome information first.
                               type="h",
                               lwd=1,
                               datacol=c("red", "royalblue4", "seagreen", "darkorchid4"),
                               chrname=0,    # the list of chromosomes to be paint with reversed order.
                               picname="chrline.jpg",    # a single name of picture(s).
                               w=1366,    # width of pictures.
                               h=768,    # height of pictures.
                               u="px",    # units of pictures.
                               reverse=FALSE)    # reverse the order of pictures.
{
  if (chrname == 0) {
    chrname <- chromosomesName(reverse = reverse)    # chrname is the vector of chromosome names.
  }
  startPicture(filename = picname,    #filenames.
               width    = w,    # picture width.
               height   = h,    # picture height.
               units    = u)    # need dev.off() later.
  basepairlength <- max(reference[, "length"])    # length of chr1 the longest.
  rectheight <- basepairlength/96    # height of rectangules.
  #### Paint blank plot
  plot(x    = 0,
       y    = 0,
       type = "n",
       xlim = c(-rectheight * 3, basepairlength),
       ylim = c(-rectheight, rectheight * length(chrname) * 2.5),
       col  = "white",
       xlab = "",
       ylab = "",
       axes = FALSE)
  #par(new = TRUE)
  #### Paint chromosomes with information
  for (n in c(1:length(chrname))) {
    righty <- rectheight * (2 * n - 1)
    lefty <- rectheight * 2 * (n - 1)
    ### Paint chromosome
    ## Chromosome part before centromere
    ### Paint data regions
    reads <- dataMerge(binread, chrname[n])
    readsmax <- max(reads$count)
    par(new = TRUE)
    lines(x = reads$site,
          y = reads$count/readsmax * rectheight +lefty)    # adjust the place of lines, multiply recheight and plus left y coordinate.
    lines(x   = c(-rectheight, basepairlength),
          y   = c(lefty, lefty),
          lty = 2)
    
  }
  #### Add label and axes for chromosomes
  text(x      = -rectheight,
       y      = ((1:length(chrname))*4 - 3) * 0.5 * rectheight,
       labels = substr(chrname,4,10), cex = 0.8, col = "red")
  axis(side   = 1,
       at     = c(0, 50000000, 100000000, 150000000, 200000000, 250000000),
       labels = c("0", "50", "100", "150", "200", "250"),
       pos    = 0)
  axis(side   = 2,
       at     = c(rectheight * 46, rectheight * 47),
       labels = c("0", "10"),
       pos    = -rectheight * 2)
  dev.off()
}


##### Function-chromosomePlotBin: paint read numbers of bins and chromosomes.

chromosomePlotBin <- function(binread,    # bin reads data, [chromosome][start][end][count].
                              data,    # list of units to paint.
                              reference,    # reference genome informations, use genome_information() to load genome information first.
                              type="h",
                              lwd=0.08,
                              datacol=c("red", "royalblue4", "seagreen", "darkorchid4"),
                              chrname=0,    # the list of chromosomes to be paint with reversed order.
                              chrcol="gray88",    # color of chromosomes.
                              cencol="darkorange",    # color of centromeres.
                              picname="binchr.jpg",    # a single name of picture(s).
                              w=1366,    # width of pictures.
                              h=768,    # height of pictures.
                              u="px",    # units of pictures.
                              single=TRUE,    # paint on one file or many files.
                              prow=6,    # if single is TRUE, the row of pictures in the file.
                              pcol=4,    # if single is TURE, the col of pictures in the file.
                              reverse=FALSE)    # reverse the order of pictures.
{
  if (chrname == 0) {
    chrname <- chromosomesName(reverse = reverse)    # chrname is the vector of chromosome names.
  }
  if (single == TRUE) {
    # if single is TRUE, chromosomes will be displayed in one file.
    picnames = picname
  } else {
    # if single is FALSE, chromosomes will be displayed in single files.
    pcol = 1
    prow = 1
    picnames = c(1:length(chrname))    # initiate picnames.
    for (chrnumi in c(1:length(chrname))) {
      # generate name sequences for a series of pictures.
      start = strsplit(picname, split = '[.]')[[1]][1]    # get the general name.
      end = rev(strsplit(picname, split = '[.]')[[1]])[1]    # get the file type.
      picnames[chrnumi] = paste(start,'_', chrname[chrnumi], '.', end, sep = '')    # make single names.
    }
  }
  for (picnumi in c(1:length(picnames))) {
    startPicture(filename = picnames[picnumi],    # filenames.
                 width    = w * pcol,    # picture width.
                 height   = h * prow,    # picture height.
                 units    = u)    # need dev.off() later.
    chrs = chrname[picnumi]    # get the picnumi's chrname value.
    if (single == TRUE) {
      layout(matrix(1:length(chrname), prow, pcol))    # setting one file to be separate into several place.
      chrs = chrname    # reset chrs to chrname.
    }
    for (chri in chrs) {
      reads <- dataMerge(binread,chri)
      limcount <- summary(reads$count)[5]
      #reads[reads$count > limcount, 3] <- limcount    # limit the reads count.
      xlimit <- reference[reference$chromosome == chri,"length"]   # max x.
      ylimit <- max(reads$count)    # max y.
      if (single == TRUE) {
        xlimit <- max(reference[,"length"])   # max x.
      }
      ybottom <- -0.1 * ylimit    # rectangle bottom y.
      ytop <- 0    # rectangle top y.
      ### Paint blank plot
      plot(x    = c(0, max(reference$length)),
           y    = c(0, 0),
           type = "s",
           lwd  = 1.5,
           xlim = c(-2500, 1.1 * xlimit),
           xlab = "",
           ylim = c(-0.15 * ylimit, 1.1 * ylimit),
           ylab = "",
           col  = "white",    # rgb arg must less than 1 rgb(222/255, 45/255, 38/255, 1).
           axes = FALSE)
      
      ### Paint chromosome
      ## Chromosome part before centromere
      rect(xleft    = 0,    #left bottom x
           ybottom  = ybottom,    # left bottom y
           xright   = reference[reference$chromosome == chri, "centromerestart"],    # right top x.
           ytop     = ytop,    # right top y.
           col      = chrcol,
           border   = chrcol)
      ## Chromosome part after centromere
      rect(xleft    = reference[reference$chromosome == chri, "centromereend"],    # left bottom x.
           ybottom  = ybottom,    # left bottom y.
           xright   = reference[reference$chromosome == chri, "length"],    # right top x.
           ytop     = ytop,    # right top y.
           col      = chrcol,
           border   = chrcol)
      ### Paint centromeres
      #rect(reference[reference$chromosome == chrname[n], "centromerestart"],    # left bottom x.
      #     ybottom,  # left bottom y . 
      #     reference[reference$chromosome == chrname[n], "centromereend"],    # right top x.
      #     ytop,  # right top y.
      #     col = cencol,
      #     border = cencol)
      ### Paint data regions
      for (m in c(1:length(data))) {
        rect(xleft    = data[[m]][data[[m]]$chromosome == chri, "start"],
             ybottom  = ybottom,
             xright   = data[[m]][data[[m]]$chromosome == chri,"end"],
             ytop     = ytop,
             col      = datacol[m],
             border   = datacol[m])
      }
      ### Paint Horizontal lines    # need revise
      lines(reads$count ~ reads$site,
            col  = rgb(8/255, 69/255, 148/255, 1),
            lwd  = lwd)
      axis(side = 2,
           at   = c(-0.3 * ylimit, 1.1 * ylimit),
           line = FALSE,
           cex  = 4)
      mtext(text = chri,
            side = 3,
            at   = 0,
            cex  = 2)
    }
    legend("right",
           legend = c(names(data)),
           pch    = 15,
           col    = c(datacol),
           border = "white",
           cex    = 1.5)
    dev.off()
  }
}


##### Function-compairMathod: compair different mathods results.
##### not finished.
compairMathod <- function(binread,    # bin reads data, list, $ [chromosome][start][end][count].
                          data,    # list of list of units to paint.
                          reference,    # reference genome informations, use genome_information() to load genome information first.
                          type="h",
                          lwd=0.08,
                          datacol=c("red", "royalblue4", "seagreen", "darkorchid4"),
                          chrname=0,    # the list of chromosomes to be paint with reversed order.
                          chrcol="gray88",    # color of chromosomes.
                          cencol="darkorange",    # color of centromeres.
                          picname="compair_mathods.jpg",    # a single name of picture(s).
                          w=1366,    # width of pictures.
                          h=768,    # height of pictures.
                          u="px",    # units of pictures.
                          single=TRUE,    # paint on one file or many files.
                          prow=6,    # if single is TRUE, the row of pictures in the file.
                          pcol=4,    # if single is TURE, the col of pictures in the file.
                          reverse=FALSE)    # reverse the order of pictures.
{
  if (chrname == 0) {
    chrname <- chromosomesName(reverse = reverse)    # chrname is the vector of chromosome names.
  }
  if (single == TRUE) {
    # if single is TRUE, chromosomes will be displayed in one file.
    picnames = picname
  } else {
    # if single is FALSE, chromosomes will be displayed in single files.
    pcol = 1
    prow = 1
    picnames = c(1:length(chrname))    # initiate picnames.
    for (chrnumi in c(1:length(chrname))) {
      # generate name sequences for a series of pictures.
      start = strsplit(picname, split = '[.]')[[1]][1]    # get the general name.
      end = rev(strsplit(picname, split = '[.]')[[1]])[1]    # get the file type.
      picnames[chrnumi] = paste(start,'_', chrname[chrnumi], '.', end, sep = '')    # make single names.
    }
  }
  for (picnumi in c(1:length(picnames))) {
    startPicture(filename = picnames[picnumi],    # filenames.
                 width    = w * pcol,    # picture width.
                 height   = h * prow,    # picture height.
                 units    = u)    # need dev.off() later.
    chrs = chrname[picnumi]    # get the pi's chrname value.
    if (single == TRUE) {
      layout(matrix(1:length(chrname), prow, pcol))    # setting one file to be separate into several place.
      chrs = chrname    # reset chrs to chrname.
    }
    for ( chri in chrs) {
      reads <- list()
      ymax <- array()
      for (binreadleni in c(1:length(binread))) {
        reads[[binreadleni]] <- dataMerge(binread[[binreadleni]], chri) 
        ymax[binreadleni] <- max(reads[[binreadleni]]$count)
      }
      #reads <- rbind(readsstart, readsend)
      #limcount <- summary(reads$count)[5]
      #reads[reads$count > limcount, 3] <- limcount    # limit the reads count.
      xlimit <- reference[reference$chromosome == chri,"length"]   # max x.
      ylimit <- max(ymax)     # max y.
      if (single == TRUE) {
        xlimit <- max(reference[,"length"])   # max x.
      }
      
      ### Paint blank plot
      plot(x    = c(0, max(reference$length)),
           y    = c(0, 0),
           type = "s",
           lwd  = 1.5,
           xlim = c(-2500, 1.1 * xlimit),
           xlab = "",
           ylim = c(-0.15 * ylimit, 1.1 * ylimit),
           ylab = "",
           col  = "white",    # rgb arg must less than 1 rgb(222/255, 45/255, 38/255, 1).
           axes = FALSE)
      for (binreadleni in c(1:length(binread))) {
        ybottom <- ylimit * (10 - binreadleni) / 10   # rectangle bottom y.
        ytop <- ylimit * (11 - binreadleni) / 10 # rectangle top y.
        ### Paint chromosome
        ## Chromosome part before centromere
        rect(xleft   = 0,    # left bottom x.
             ybottom = ybottom,    # left bottom y.
             xright  = reference[reference$chromosome == chri, "centromerestart"],    # right top x.
             ytop    = ytop,    # right top y.
             col     = chrcol,
             border  = chrcol)
        ## Chromosome part after centromere
        rect(xleft   = reference[reference$chromosome == chri, "centromereend"],    # left bottom x.
             ybottom = ybottom,    # left bottom y.
             xright  = reference[reference$chromosome == chri, "length"],    # right top x.
             ytop    = ytop,    # right top y.
             col     = chrcol,
             border  = chrcol)
        ###paint centromeres
        #rect(reference[reference$chromosome == chrname[n], "centromerestart"],    #left bottom x
        #     ybottom,  #left bottom y
        #     reference[reference$chromosome == chrname[n], "centromereend"],    #right top x
        #     ytop,  #right top y
        #     col = cencol,
        #     border = cencol)
        ### Paint data regions
        for (dataleni in c(1:length(data[[binreadleni]]))) {
          # check whether the data2 has data in such situation.
          if (length(data[[binreadleni]][[dataleni]][data[[binreadleni]][[dataleni]]$chromosome == chri, "start"]) == 0) {
            next
          }
          if (length(data[[binreadleni]][[dataleni]][data[[binreadleni]][[dataleni]]$chromosome == chri, "end"]) == 0) {
            next
          }
          rect(xleft   = data[[binreadleni]][[dataleni]][data[[binreadleni]][[dataleni]]$chromosome == chri, "start"],
               ybottom = ybottom,
               xright  = data[[binreadleni]][[dataleni]][data[[binreadleni]][[dataleni]]$chromosome == chri, "end"],
               ytop    = ytop,
               col     = datacol[dataleni],
               border  = datacol[dataleni])
        }
        ### Paint Horizontal lines
        par(new = TRUE)
        plot(reads[[binreadleni]]$count ~ reads[[binreadleni]]$site,
             type = type,
             col  = rgb(8/255, 69/255, 148/255, 1),
             lwd  = lwd,
             xlim = c(-2500, 1.1 * xlimit),
             ylim = c(-0.3 * ylimit, 1.1 * ylimit),
             ylab = "",
             xlab = "",
             axes = FALSE)
        
      }
      axis(side = 2,
           at   = c(-0.3 * ylimit, 1.1 * ylimit),
           line = FALSE)
      mtext(chri,
            side = 3,
            at   = 0,
            cex  = 1.5)
    }
    legend("right",
           legend = c(names(data)),
           pch    = 15,
           col    = c(datacol),
           border = "white",
           cex    = 1.5)
    dev.off()
  }
}


##### Function-compair2Mathods: compair different mathods results.

compair2Mathods <- function(methodname,
                            binread1,    # bin reads data, list, $ [chromosome][start][end][count].
                            data1,    # list of list of units to paint.
                            binread2,    # bin reads data, list, $ [chromosome][start][end][count].
                            data2,    # list of list of units to paint.
                            reference,    # reference genome informations, use genome_information() to load genome information first.
                            col1=c("royalblue4", "seagreen", "darkorchid4"),
                            col2=c("red", "seagreen", "darkorchid4"),
                            type="h",    # type, choose h for histogram lines, l for lines
                            lwd=0.08,    # line width.
                            #datacol=c("royalblue4", "seagreen", "darkorchid4"),
                            chrname=0,    # the list of chromosomes to be paint with reversed order.
                            chrcol="gray88",    # color of chromosomes.
                            cencol="darkorange",    # color of centromeres.
                            picname="compair_mathods.jpg",    # a single name of picture(s).
                            w=1366,    # width of pictures.
                            h=768,    # height of pictures.
                            u="px",    # units of pictures.
                            single=TRUE,    # paint on one file or many files.
                            prow=6,    # if single is TRUE, the row of pictures in the file.
                            pcol=4,    # if single is TURE, the col of pictures in the file.
                            reverse=FALSE)    # reverse the order of pictures.
{
  if (chrname == 0) {
    chrname <- chromosomesName(reverse = reverse)    # chrname is the vector of chromosome names.
  }
  if (single == TRUE) {
    # if single is TRUE, chromosomes will be displayed in one file.
    picnames = picname
  } else {
    # if single is FALSE, chromosomes will be displayed in single files.
    pcol = 1
    prow = 1
    picnames = c(1:length(chrname))    # initiate picnames.
    for (chrnumi in c(1:length(chrname))) {
      # generate name sequences for a series of pictures.
      start = strsplit(picname, split = '[.]')[[1]][1]    # get the general name.
      end = rev(strsplit(picname, split = '[.]')[[1]])[1]    # get the file type.
      picnames[chrnumi] = paste(start,'_', chrname[chrnumi], '.', end, sep = '')    # make single names.
    }
  }
  for (picnumi in c(1:length(picnames))) {
    startPicture(filename = picnames[picnumi],    # filenames.
                 width    = w * pcol,    # picture width.
                 height   = h * prow,    # picture height.
                 units    = u)    # need dev.off() later.
    chrs = chrname[picnumi]    # get the pi's chrname value.
    if (single == TRUE) {
      layout(matrix(1:length(chrname), prow, pcol))    # setting one file to be separate into several place.
      chrs = chrname    # reset chrs to chrname.
    }
    for ( chri in chrs) {
      #### Manage data
      reads1 <- dataMerge(binread1, chri)    # merge data.
      reads1 <- dataFrameValueNorm(reads1, 3)    # normalize data.
      ymax1 <- max(reads1$count)
      reads2 <- dataMerge(binread2, chri) 
      reads2 <- dataFrameValueNorm(reads2, 3)
      ymax2 <- max(reads2$count)
      
      xlimit <- reference[reference$chromosome == chri,"length"]   # max x.
      ylimit <- max(ymax1, ymax2)     # max y.
      if (single == TRUE) {
        xlimit <- max(reference[,"length"])   # max x.
      }
      ### Paint blank plot
      plot(x    = c(0, max(reference$length)),
           y    = c(0, 0),
           type = "s",
           lwd  = 1.5,
           xlim = c(-2500, 1.1 * xlimit),
           xlab = "",
           ylim = c(-0.5 * ylimit, 1.1 * ylimit),
           ylab = "",
           col  = "white",    # rgb arg must less than 1 rgb(222/255, 45/255, 38/255, 1).
           axes = FALSE)
      ### Calculate ybottoms and ytops
      ybottom1 <- -1 * ylimit / 10   # rectangle bottom y.
      ytop1 <- -1 * 0.5 * ylimit / 10  # rectangle top y.
      ybottom2 <- -1 * 2 * ylimit / 10   # rectangle bottom y.
      ytop2 <- -1 * 1.5 * ylimit / 10  # rectangle top y.
      ### Paint chromosome for method 1
      ## Chromosome part before centromere
      rect(xleft   = 0,    # left bottom x.
           ybottom = ybottom1,    # left bottom y.
           xright  = reference[reference$chromosome == chri, "centromerestart"],    # right top x.
           ytop    = ytop1,    # right top y.
           col     = chrcol,
           border  = chrcol)
      ## Chromosome part after centromere
      rect(xleft   = reference[reference$chromosome == chri, "centromereend"],    # left bottom x.
           ybottom = ybottom1,    # left bottom y.
           xright  = reference[reference$chromosome == chri, "length"],    # right top x.
           ytop    = ytop1,    # right top y.
           col     = chrcol,
           border  = chrcol)
      ## Label of method1
      #text(x      = -5000,
      #     y      = ybottom1,
      #     labels = methodname[1],
      #     cex    = 1.5)
      ### Paint chromosome for method 2
      ## Chromosome part before centromere
      rect(xleft   = 0,    # left bottom x.
           ybottom = ybottom2,    # left bottom y.
           xright  = reference[reference$chromosome == chri, "centromerestart"],    # right top x.
           ytop    = ytop2,    # right top y.
           col     = chrcol,
           border  = chrcol)
      ## Chromosome part after centromere
      rect(xleft   = reference[reference$chromosome == chri, "centromereend"],    # left bottom x.
           ybottom = ybottom2,    # left bottom y.
           xright  = reference[reference$chromosome == chri, "length"],    # right top x.
           ytop    = ytop2,    # right top y.
           col     = chrcol,
           border  = chrcol)
      ## Label of method1
      #text(x      = -5000,
      #     y      = ybottom2,
      #     labels = methodname[2],
      #     cex    = 1.5)
      ###paint centromeres
      #rect(reference[reference$chromosome == chrname[n], "centromerestart"],    #left bottom x
      #     ybottom,  #left bottom y
      #     reference[reference$chromosome == chrname[n], "centromereend"],    #right top x
      #     ytop,  #right top y
      #     col = cencol,
      #     border = cencol)
      ### Paint data 1 regions 
      for (dataleni in c(1:length(data1))) {
        # check whether the data2 has data in such situation.
        if (length(data1[[dataleni]][data1[[dataleni]]$chromosome == chri, "start"]) == 0) {
          next
        }
        if (length(data1[[dataleni]][data1[[dataleni]]$chromosome == chri, "end"]) == 0) {
          next
        }
        rect(xleft   = data1[[dataleni]][data1[[dataleni]]$chromosome == chri, "start"],
             ybottom = ybottom1,
             xright  = data1[[dataleni]][data1[[dataleni]]$chromosome == chri, "end"],
             ytop    = ytop1,
             col     = col1[dataleni],
             border  = col1[dataleni])
      }
      ### Paint data 2 regions
      for (dataleni in c(1:length(data2))) {
        # check whether the data2 has data in such situation.
        if (length(data2[[dataleni]][data2[[dataleni]]$chromosome == chri, "start"]) == 0) {
          next
        }
        if (length(data2[[dataleni]][data2[[dataleni]]$chromosome == chri, "end"]) == 0) {
          next
        }
        rect(xleft   = data2[[dataleni]][data2[[dataleni]]$chromosome == chri, "start"],
             ybottom = ybottom2,
             xright  = data2[[dataleni]][data2[[dataleni]]$chromosome == chri, "end"],
             ytop    = ytop2,
             col     = col2[dataleni],
             border  = col2[dataleni])
      }
      ### Paint bin reads 
      ## 1
      lines(reads1$count ~ reads1$site, 
            col  = col1[1],
            lwd  = lwd)
      ## 2
      lines(reads2$count ~ reads2$site,
            type = type,
            col  = col2[1],
            lwd  = lwd) 
      
      axis(side = 2,
           at   = c(-0.3 * ylimit, 1.1 * ylimit),
           line = FALSE)
      mtext(chri,
            side = 3,
            at   = 0,
            cex  = 1.5)
    }
    if (single == FALSE) {
      legend("right",
             legend = c(names(data1), names(data2)),
             pch    = 15,
             col    = c(col1[1:length(names(data1))], col2[1:length(names(data2))]),
             border = "white",
             cex    = 1.5)
    }
    
    dev.off()
  }
}


##### Function-binReadDensity: show the reads density.

binReadDensity <- function(binread,    # input the the binread data.frame with count at the 4th column.
                           chrname=0,    # the list of chromosomes to be paint with reversed order.
                           picname="density.jpg",    # a single name of picture(s).
                           w=1366,    # width of pictures.
                           h=768,    # height of pictures.
                           u="px",    # units of pictures.
                           reverse=TRUE)
{
  if (chrname == 0) {
    chrname <- chromosomesName(reverse = reverse)    # chrname is the vector of chromosome names.
  }
  
  picnames = c(1:length(chrname))    # initiate picnames.
  for (chrnumi in c(1:length(chrname))) {
    # generate name sequences for a series of pictures.
    start = strsplit(picname, split = '[.]')[[1]][1]    # get the general name.
    end = rev(strsplit(picname, split = '[.]')[[1]])[1]    # get the file type.
    picnames[chrnumi] = paste(start,'_', chrname[chrnumi], '.', end, sep = '')    # make single names.
  }
  for (picnumi in c(1:length(picnames))) {
    startPicture(filename = picnames[picnumi],    # filenames.
                 width    = w,    # picture width.
                 height   = h,    # picture height.
                 units    = u)    # need dev.off() later.
    rd.ds <- density(binread[binread$chromosome == chrname[picnumi], "count"])
    plot(rd.ds)
    dev.off()
  }
}
#####EOF
