#!/usr/bin/env Rscript
#-*-coding:utf-8-*-

#####Functions to display chromosomes information
#####Author: Wolfson


#####function-readgenomeinfo:change reference genome information to rigth format

genome_infomation <- function(data,  #data frame of genome informations: columns(chr)(centrstart)(centrend)(length) 
                              cnchr,  #column number of chr
                              cncentrstart,  #column number of start of centromere
                              cncentrend,  #column number of end of centromere
                              cnlength)  #column number of length of chromosome
{
  tmp <- cbind(data[cnchr], data[cncentrstart], data[cncentrend], data[cnlength])
  colnames(tmp) <- c("chromosome", "centromerestart", "centromereend", "length")
  return(tmp)
}

#####function-domain_data_list:format data to [name][chromosome][start][end]

domain_data_list <- function(data,  #list of data file
                             dataname,  #vector of data names
                             datacol)  #table of column numbers of data list
{
  datalen <- length(data)
  if (length(dataname) != datalen){
    print("Data and dataname are not the same length.","Please check again!")
  }
  tmpdata <- list()
  for (n in 1:datalen){
    tmp <- data[[n]][,datacol[n,]]
    names(tmp) <- c("name", "chromosome", "start", "end")
    tmpdata[[dataname[n]]] <- tmp
  }
  names(tmpdata) <- dataname
  return(tmpdata)
}

#####function-get_file_type:get the type of input filename

get_file_type <- function(filename)
{
  filetype <- rev(strsplit(filename, split = '[.]')[[1]])[1]
  return(filetype)
}

#####function-choose_function:choose functions according to words

choose_function <- function(filename,  #filename
                            funcs,  #functions
                            words)  #reference words vector
{
  wordslen <- length(words)
  funcslen <- length(funcs)
  if (wordslen != funcslen){
    print("Functions and words are not the same length.")
    return(0)
  }
  filetype <- get_file_type(filename)
  for (i in 1:funcslen){
    if(filetype == words[i]) {
      funcs[i]
      func <- funcs[i]
      break
    }
  }
  if (is.na(func)){
    print("Function was not found.")
  }
  return(func)
}


#####function-merge_startend: merge the start and end counts for paint

merge_startend <- function(binread,
                           chr,  #chromosome
                           stc = c(1:2,4),  #order: chr start count
                           edc = c(1,3:4),  #order: chr end count
                           bstart = 0)  #set the bin start from 0 or minit of bin
{
  readsstart <- binread[binread$chromosome == substr(chr, 4, 10), stc]  #get the bin start data of chromosome ci
  readsend <- binread[binread$chromosome == substr(chr, 4, 10), edc]  #get the bin end data of chromosome ci
  names(readsstart) <- c("chromosome", "site", "count")
  names(readsend) <- c("chromosome", "site", "count")
  binsize <- min(readsend$site) - min(readsstart$site)  #bin size
  binstart <- bstart  #the start of bin site
  binend <- max(readsend$site)  #the end of bin site
  binnum <- (binend - binstart) / binsize  #number of bins
  bintablestart <- as.data.frame(matrix(0,nrow = binnum, ncol = 3))  #"chromosome", "site", "count"
  bintableend <- as.data.frame(matrix(0,nrow = binnum, ncol = 3))  #"chromosome", "site", "count"
  bintablestart[,1] <- ci
  bintableend[,1] <- ci
  bintablestart[,2] <- seq(binstart, binend - binsize, by = binsize)
  bintablestart[,2] <- seq(binstart + binsize, binend, by = binsize)
  for (n in readsstart$site){
    bintablestart[bintablestart$V2 == n, 3] <- readsstart[readsstart$site == n, 3]
    bintableend[bintableend$V2 ==  n, 3] <- readsend[readsstart$site == n, 3]
  }
  tmp <- rbind(bintablestart, bintableend)
  names(tmp) <- c("chromosome", "site", "count")
  reads <- tmp[order(tmp$site),]  #make the data frame with a right order
  return(reads)
}
#####function-start_picture:change the picture output type

#start_picture <- function(file, width = 800, height = 600, units = "px")
#{  #this function is used to change different kinds of pictures
#                                        #remember to dev.off() after the picture
#  picfunc <- choose_function(file, c(jpeg, png, tiff, bmp), c("jpg", "png", "tiff", "bmp"))
#  return(picfunc)
#    #picfunc(filename = file,  #file name, end with picture style
     #     width = width,  #width
      #    height = height,  #height
       #   units = units))  #unit for width and height, usually "px"
#}




#####function-chromosomes:paint chromomes with different kinds of units

chromosomes <- function(data,  #list of units to paint
                        reference,  #reference genome informations, use genome_information() to load genome information first
                        datacol = c("red", "royalblue4", "darkolivegreen", "blanchedalmond"),
                        chri = 0,  #the list of chromosomes to be paint with reversed order
                        chrcol = "deepskyblue",  #color of chromosomes
                        cencol = "darkorange",  #color of centromeres
                        picname = "chromosomes.jpg",
                        w = 1366,
                        h = 768,
                        u = "px",
                        reverse = FALSE)
{
  if (chri == 0){
    chri <- as.vector(c("chrY", "chrX", "chr22", "chr21", "chr20", "chr19", 
                        "chr18", "chr17", "chr16", "chr15", "chr14", "chr13", 
                        "chr12", "chr11", "chr10", "chr9", "chr8", "chr7", 
                        "chr6", "chr5", "chr4", "chr3", "chr2", "chr1"))  #chromosome order
  }
  if (reverse == TRUE){
      chri <- rev(chri)
  }
  #start_picture <- choose_function(picname, c(jpeg, png, tiff, bmp), c("jpg", "png", "tiff", "bmp"))
  jpeg(filename = picname,
                width = w,
                height = h,
                units = u)  #need dev.off() later
  basepairlength <- max(reference[, "length"])  #length of chr1 the longest
  rectheight <- basepairlength/96 #height of rectangules
  ####paint blank plot
  plot(0,
       0,
       type = "n",
       xlim = c(0, basepairlength),
       ylim = c(0, rectheight * length(chri) * 2),
       col = "white",
       xlab = "",
       ylab = "",
       axes = FALSE)
  #par(new = TRUE)
  ####paint chromosomes with information
  for (n in c(1:length(chri))){
    righty <- rectheight * (2 * n - 1)
    lefty <- rectheight * 2 * (n - 1)
    ###paint chromosome
    ##chromosome part before centromere
    rect(0,  #left bottom x
         lefty + 0.25 * rectheight,  #left bottom y
         reference[reference$chromosome == chri[n], "centromerestart"],  #right top x
         righty - 0.25 * rectheight,  #right top y
         col = chrcol,
         border = chrcol)
    ##chromosome part after centromere
    rect(reference[reference$chromosome == chri[n], "centromereend"],  #left bottom x
         lefty + 0.25 * rectheight,  #left bottom y
         reference[reference$chromosome == chri[n], "length"],  #right top x
         righty - 0.25 * rectheight,  #right top y
         col = chrcol,
         border = chrcol)
    ###paint centromeres
    rect(reference[reference$chromosome == chri[n], "centromerestart"],
         lefty + 0.25 * rectheight,
         reference[reference$chromosome == chri[n], "centromereend"],
         righty - 0.25 * rectheight,
         col = cencol,
         border = cencol)
    ###paint data regions
    for (m in c(1:length(data))){
      rect(data[[m]][data[[m]]$chromosome == chri[n], "start"],
           lefty,
           data[[m]][data[[m]]$chromosome == chri[n], "end"],
           righty,
           col = datacol[m],
           border = datacol[m])
    }
  }
    ###paint satellites
    #rect(satellite$genoStart[satellite$genoName == chri[n]],
    #     lefty,
    #     satellite$genoEnd[satellite$genoName == chri[n]],
    #     righty,
    #     col = satcol,
    #     border = satcol)
    #}
  ####add label for chromosomes
  text(x = -rectheight * 3,
       y=((1:length(chri))*4 - 3) * 0.5 * rectheight,
       labels = chri, cex = 1.5)
  ####legend
  legend("right",
         c(names(data)[c(1,length(data))], "centromere"),
         pch = 15,
         col = c(c(1,datacol[length(data)]), cencol),
         border = "white",
         cex = 1.5)
  dev.off()
}


#####function-chrreads:paint chromomes with reads

chrreads <- function(binread,  #bin reads data, [chromosome][start][end][count]
                     data,  #list of units to paint
                     reference,  #reference genome informations, use genome_information() to load genome information first
                     type = "h",
                     lwd = 1,
                     datacol = c("red", "royalblue4", "seagreen", "darkorchid4"),
                     chri = 0,  #the list of chromosomes to be paint with reversed order
                     picname = "chr.jpg",  #a single name of picture(s)
                     w = 1366,  #width of pictures
                     h = 768,  #height of pictures
                     u = "px",  #units of pictures
                     reverse = FALSE)  #reverse the order of pictures
{
  if (chri == 0){
    chri <- as.vector(c("chrY", "chrX", "chr22", "chr21", "chr20", "chr19", 
                        "chr18", "chr17", "chr16", "chr15", "chr14", "chr13", 
                        "chr12", "chr11", "chr10", "chr9", "chr8", "chr7", 
                        "chr6", "chr5", "chr4", "chr3", "chr2", "chr1"))
    
  }
  if (reverse == TRUE){
    chri <- rev(chri)  #reverse chromosome names vector
  }
  jpeg(picname,  #filenames
       width = w,  #picture width
       height = h,  #picture height
       units = u)  #need dev.off() later
  basepairlength <- max(reference[, "length"])  #length of chr1 the longest
  rectheight <- basepairlength/96 #height of rectangules
  ####paint blank plot
  plot(0,
       0,
       type = "n",
       xlim = c(-rectheight * 3, basepairlength),
       ylim = c(-rectheight, rectheight * length(chri) * 2.5),
       col = "white",
       xlab = "",
       ylab = "",
       axes = FALSE)
  #par(new = TRUE)
  ####paint chromosomes with information
  for (n in c(1:length(chri))){
    righty <- rectheight * (2 * n - 1)
    lefty <- rectheight * 2 * (n - 1)
    ###paint chromosome
    ##chromosome part before centromere
    ###paint data regions
    reads <- merge_startend(binread, chri[n])
    readsmax <- max(reads$count)
    par(new = TRUE)
    lines(x = reads$site,
          y = reads$count/readsmax * rectheight +lefty)  #adjust the place of lines, multiply recheight and plus left y coordinate
    lines(x = c(-rectheight, basepairlength),
          y = c(lefty, lefty),
          lty = 2)

  }
  ####add label and axes for chromosomes
  text(x = -rectheight,
       y=((1:length(chri))*4 - 3) * 0.5 * rectheight,
       labels = substr(chri,4,10), cex = 0.8, col = "red")
  axis(side = 1,
       at = c(0, 50000000, 100000000, 150000000, 200000000, 250000000),
       labels = c("0", "50", "100", "150", "200", "250"),
       pos = 0)
  axis(side = 2,
       at = c(rectheight * 46, rectheight * 47),
       labels = c("0", "10"),
       pos = -rectheight * 2)
  dev.off()
}


#####function-binreads_chromosome: paint read numbers of bins and chromosomes

binreads_chromosome <- function(binread,  #bin reads data, [chromosome][start][end][count]
                                data,  #list of units to paint
                                reference,  #reference genome informations, use genome_information() to load genome information first
                                type = "h",
                                lwd = 0.08,
                                datacol = c("red", "royalblue4", "seagreen", "darkorchid4"),
                                chri = 0,  #the list of chromosomes to be paint with reversed order
                                chrcol = "gray88",  #color of chromosomes
                                cencol = "darkorange",  #color of centromeres
                                picname = "binchr.jpg",  #a single name of picture(s)
                                w = 1366,  #width of pictures
                                h = 768,  #height of pictures
                                u = "px",  #units of pictures
                                single = TRUE,  #paint on one file or many files
                                prow = 6,  #if single is TRUE, the row of pictures in the file
                                pcol = 4,  #if single is TURE, the col of pictures in the file
                                reverse = FALSE)  #reverse the order of pictures
{
  if (chri == 0){
    chri <- as.vector(c("chrY", "chrX", "chr22", "chr21", "chr20", "chr19", 
                        "chr18", "chr17", "chr16", "chr15", "chr14", "chr13", 
                        "chr12", "chr11", "chr10", "chr9", "chr8", "chr7", 
                        "chr6", "chr5", "chr4", "chr3", "chr2", "chr1"))
    
  }
  if (reverse == TRUE){
    chri <- rev(chri)  #reverse chromosome names vector
  }
  if (single == TRUE){
    #if single is TRUE, chromosomes will be displayed in one file
    picnames = picname
  }
  else{
    #if single is FALSE, chromosomes will be displayed in single files
    pcol = 1
    prow = 1
    picnames = c(1:length(chri))  #initiate picnames
    for (i in c(1:length(chri))){
      #generate name sequences for a series of pictures
      start = strsplit(picname, split = '[.]')[[1]][1]  #get the general name
      end = rev(strsplit(picname, split = '[.]')[[1]])[1]  #get the file type
      picnames[i] = paste(start,'_', chri[i], '.', end, sep = '')  #make single names
    }
  }
  for (pi in c(1:length(picnames))){
    #start_picture <- choose_function(picname, c(jpeg, png, tiff, bmp), c("jpg", "png", "tiff", "bmp"))
    jpeg(picnames[pi],  #filenames
         width = w * pcol,  #picture width
         height = h * prow,  #picture height
         units = u)  #need dev.off() later
    chrs = chri[pi]  #get the pi's chri value
    if (single == TRUE){
      layout(matrix(1:length(chri), prow, pcol))  #setting one file to be separate into several place
      chrs = chri  #reset chrs to chri
    }
    for ( ci in chrs){
      reads <- merge_startend(binread,ci)
      #reads <- rbind(readsstart, readsend)
      limcount <- summary(reads$count)[5]
      #reads[reads$count > limcount, 3] <- limcount  #limit the reads count
      xlimit <- reference[reference$chromosome == ci,"length"] #max x
      ylimit <- max(reads$count)  #max y
      if (single == TRUE){
        xlimit <- max(reference[,"length"]) #max x
      }
      ###check
      #str(bintablestart)
      #str(reads)
      #max(reads$count)
      #levels(reads$count)
      #print(binstart)
      #print(binend)
      #print(max(reads$count))
      #print(ylimit)
      #print(xlimit)
      ybottom <- -0.1 * ylimit  #rectangle bottom y
      ytop <- 0  #rectangle top y
      ###paint blank plot
      plot(x = c(0, max(reference$length)),
           y = c(0, 0),
           type = "s",
           lwd = 1.5,
           xlim = c(-2500, 1.1 * xlimit),
           xlab = "",
           ylim = c(-0.15 * ylimit, 1.1 * ylimit),
           ylab = "",
           col = "white",  #rgb arg must less than 1 rgb(222/255, 45/255, 38/255, 1)
           axes = FALSE)
      
      ###paint chromosome
      ##chromosome part before centromere
      rect(0,  #left bottom x
           ybottom,  #left bottom y
           reference[reference$chromosome == ci, "centromerestart"],  #right top x
           ytop,  #right top y
           col = chrcol,
           border = chrcol)
      ##chromosome part after centromere
      rect(reference[reference$chromosome == ci, "centromereend"],  #left bottom x
           ybottom,  #left bottom y
           reference[reference$chromosome == ci, "length"],  #right top x
           ytop,  #right top y
           col = chrcol,
           border = chrcol)
      ###paint centromeres
      #rect(reference[reference$chromosome == chri[n], "centromerestart"],    #left bottom x
      #     ybottom,  #left bottom y
      #     reference[reference$chromosome == chri[n], "centromereend"],  #right top x
      #     ytop,  #right top y
      #     col = cencol,
      #     border = cencol)
      ###paint data regions
      for (m in c(1:length(data))){
        rect(data[[m]][data[[m]]$chromosome == ci, "start"],
             ybottom,
             data[[m]][data[[m]]$chromosome == ci,"end"],
             ytop,
             col = datacol[m],
             border = datacol[m])
      }
      ###paint Horizontal lines
      par(new = TRUE)
      plot(reads$count ~ reads$site,
           type = type,
           col = rgb(8/255, 69/255, 148/255, 1),
           lwd = lwd,
           xlim = c(-2500, 1.1 * xlimit),
           ylim = c(-0.3 * ylimit, 1.1 * ylimit),
           ylab = "",
           xlab = "",
           axes = FALSE)
      axis(side = 2,
           at = c(-0.3 * ylimit, 1.1 * ylimit),line = FALSE)
      mtext(ci,
            side = 3,
            at = 0,
            cex = 1.5)
    }
    legend("right",
           c(names(data)),
           pch = 15,
           col = c(datacol),
           border = "white",
           cex = 1.5)
    dev.off()
  }
}

#####EOF