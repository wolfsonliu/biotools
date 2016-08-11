# -*- coding: utf-8 -*-

################################################################################
##                                                                            ##
##                             Circos plot                                    ##
##                                                                            ##
## Author     : Wolfson                                                       ##
## Date       : 2016-07-05                                                    ##
## Modified   : 2016-07-05                                                    ##
## Version    : 0.01                                                          ##
## Description:                                                               ##
##                                                                            ##
##----------------------------------------------------------------------------##
##                                                                            ##
##                                                                            ##
################################################################################

suppressPackageStartupMessages(require(methods,      quietly = TRUE))
suppressPackageStartupMessages(require(grDevices,    quietly = TRUE))
suppressPackageStartupMessages(require(RColorBrewer, quietly = TRUE))
suppressPackageStartupMessages(require(colorspace,   quietly = TRUE))
suppressPackageStartupMessages(require(reshape2,   quietly = TRUE))

changeMatrix <- function(change) {
###{{{ Function: changeMatrix

    function(data) {
       storecall <- match.call()
       func.name <- as.character(storecall[1]) # get the function name.
       ## Check input data
       if (!is.matrix(data)) {
           stop(paste0(func.name, ": Please enter a matrix as input data."))
       } else {}
       if (median(data) > 0 && median(data) < 1) {
           change(data * 10)
       } else {
           change(data)
       }          
    }

###}}}
}

ceilMatrix <- changeMatrix(ceiling)

table2KaryotypeLength <- function(data) {
###{{{ Function: table2KaryotypeLength

    ## 
    ## Function: table2KaryotypeLength
    ## This function is used to calculate the length of item should be plotted
    ## in circos.
    ## 
    storecall <- match.call()
    func.name <- as.character(storecall[1]) # get the function name.
    ## Check input data
    if (!is.matrix(data)) {
        stop(paste0(func.name, ": Please enter a matrix as input data."))
    } else {}
    ## The input matrix should have rownames and colnames.
    lapply(
        c(rownames, colnames),
        function(x) {
            if (is.null(x(data))) {
                stop(
                    paste0(
                        func.name,
                        ": Please give data ",
                        deparse(substitute(x)),
                        "."
                    )                    
                )
            } else {}
        }
    )
    ## Calculate the length of item
    data        <- ceilMatrix(data)     # scale data and ceiling float to int.
    l.collength <- colSums(data)        # End of links.
    l.rowlength <- colSums(t(data))     # Start of links.
    l.item      <- unique(c(rownames(data), colnames(data)))
    l.length    <- unlist(
        ## Sum the start and end length of the same item
        lapply(
            l.item,
            function(x) {
                len        <- sum(l.collength[x], l.rowlength[x], na.rm = TRUE)
                names(len) <- x
                len
            }
        )
    )
    l.length

###}}}
}


makeKaryotypeColor <- function(name,
                               alpha  = 0.4) {
###{{{ Function: makeKaryotypeColor

    ##
    ## Function: makeKaryotypeColor
    ## This function is used to assign colors to item names in circos.
    ##
    storecall <- match.call()
    func.name <- as.character(storecall[1]) # get the function name.
    ## Check input data
    if (!is.character(name)) {
        stop(paste0(func.name, ": Please enter a character vector as input."))
    } else {}
    ## if col1 or col2 not assigned, make colors.
    col1 <- colorspace::hex2RGB(
        substr(rainbow(length(name)),start = 1, stop = 7)
    )
    col2 <- colorspace::RGB(1, 1, 1)
    l.color <- colorspace::hex(
        colorspace::mixcolor(alpha, color1 = col1, color2 = col2)
    )
    names(l.color) <- name
    if (length(l.color) != length(name)) {
        stop(paste0(func.name, ": Color length and name length not equal."))
    } else {}

    l.color

###}}}
}


outputColorConf <- function(color,
                            colornameprefix = "color",
                            file) {
###{{{ Function: outputColorConf

    ##
    ## Function: outputColorConf
    ## output the color configure file
    ##
    storecall <- match.call()
    func.name <- as.character(storecall[1]) # get the function name.
    ## compare data
    color.conf <- data.frame(                    
        name = paste0(colornameprefix, tolower(names(color))),
        sign = rep("=", length(color)),
        rgb  = apply(
            t(col2rgb(color)),
            1,
            function(x) paste0(x, collapse = ",")
        )
    )
    rownames(color.conf) <- names(color)
    if (!is.null(as.list(storecall[-1])[["file"]])) {
        write.table(
            color.conf,
            file      = file,
            quote     = FALSE,
            sep       = " ",
            row.names = FALSE,
            col.names = FALSE
        )        
    } else {}
    color.conf

###}}}
}


outputKaryotype <- function(data,
                            alpha = 0.4,
                            file) {
###{{{ Function: outputKaryotypeConf

    ## 
    ## Function: outputKaryotypeConf
    ## output karyotype data file of circos
    ## 
    storecall <- match.call()
    func.name <- as.character(storecall[1]) # get the function name.
    ## Check input data
    if (!is.matrix(data)) {
        stop(paste0(func.name, ": Please enter a matrix as input data."))
    } else {}
    ## The input matrix should have rownames and colnames.
    lapply(
        c(rownames, colnames),
        function(x) {
            if (is.null(x(data))) {
                stop(
                    paste0(
                        func.name,
                        ": Please give data ",
                        deparse(substitute(x)),
                        "."
                    )                    
                )
            } else {}
        }
    )
    ## Calculate length and color
    l.length <- table2KaryotypeLength(data)
    l.color  <- makeKaryotypeColor(names(l.length), alpha)
    l.color[colnames(data)] <- gray(
        seq_along(colnames(data)) / (length(colnames(data)) + 1)
    )
    l.karyotype <- data.frame(
        Type  = rep("chr", length(l.length)),
        Place = rep("-", length(l.length)),
        ID    = names(l.length),
        Lable = names(l.length),
        Start = rep(0, length(l.length)),
        End   = l.length,
        Color = outputColorConf(l.color)[names(l.length), "name"]
    )
    ## Write to file if file path is gaven
    if (!is.null(as.list(storecall[-1])[["file"]])) {
        write.table(
            l.karyotype,
            file      = file,
            quote     = FALSE,
            sep       = " ",
            row.names = FALSE,
            col.names = FALSE
        )
    } else {}    
    l.karyotype

###}}}
}


outputLink <- function(data,
                       file,
                       setcolor = TRUE,
                       ...) {
###{{{

    storecall <- match.call()
    func.name <- as.character(storecall[1]) # get the function name.
    ## Check input data
    if (!is.matrix(data)) {
        stop(paste0(func.name, ": Please enter a matrix as input data."))
    } else {}
    ## The input matrix should have rownames and colnames.
    lapply(
        c(rownames, colnames),
        function(x) {
            if (is.null(x(data))) {
                stop(
                    paste0(
                        func.name,
                        ": Please give data ",
                        deparse(substitute(x)),
                        "."
                    )                    
                )
            } else {}
        }
    )
    ## Make linkage file
    l.length                <- table2KaryotypeLength(data)
    l.color                 <- makeKaryotypeColor(names(l.length), ...)
    l.color[colnames(data)] <- gray(
        seq_along(colnames(data)) / (length(colnames(data)) + 1)
    )
    l.record                <- data.frame(
        mark   = rep(0, length(l.length)),
        length = l.length
    )
    rownames(l.record)      <- names(l.length)
    l.interaction           <- melt(ceilMatrix(data))
    l.interaction           <- l.interaction[l.interaction$value != 0,]
    l.link                  <- matrix(nrow = dim(l.interaction)[1], ncol = 6)
    for (i in seq(1, dim(l.interaction)[1])) {
        one.start  <- l.record[as.character(l.interaction[i, "Var1"]),"mark"]
        one.end    <- as.numeric(one.start) + as.numeric(l.interaction[i, "value"])
        l.record[as.character(l.interaction[i, "Var1"]),"mark"] <- one.end
        two.start  <- l.record[as.character(l.interaction[i, "Var2"]), "mark"]
        two.end    <- as.numeric(two.start) + as.numeric(l.interaction[i,"value"])
        l.record[as.character(l.interaction[i,"Var2"]), "mark"] <- two.end
        l.link[i,] <- c(
            as.character(l.interaction[i,"Var1"]),
            one.start,
            one.end,
            as.character(l.interaction[i,"Var2"]),
            two.start,
            two.end
        )
    }
    l.linkdata <- data.frame(
        One       = l.link[,1],
        One.start = l.link[,2],
        One.end   = l.link[,3],
        Two       = l.link[,4],
        Two.start = l.link[,5],
        Two.end   = l.link[,6]
    )
    if (setcolor) {
        l.linkdata$color <- paste0(
            "color=",
            as.character(outputColorConf(l.color, ...)[l.linkdata$One, "name"])
        )
    } else {}
    ## Write to file if file path is gaven
    if (!is.null(as.list(storecall[-1])[["file"]])) {
        write.table(
            l.linkdata,
            file      = file,
            quote     = FALSE,
            sep       = " ",
            row.names = FALSE,
            col.names = FALSE
        )
    } else {}    
    l.linkdata
    
###}}}
}



