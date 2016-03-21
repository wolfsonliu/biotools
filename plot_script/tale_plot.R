#! /bin/env R

################################################################################
## Author:      Wolfson                                                       ##
## Date:        2016-02-26                                                    ##
## Modified:    2016-02-29                                                    ##
## Description: Plot the heatmap for TALE                                     ##
################################################################################

tale.heatmap <- function(data,
                         xlab         = "",
                         ylab         = "",
                         legend.label = c("", "", "", ""),
                         label.size   = 1,
                         title.size   = 1,
                         bias         = 0.8,
                         col          = c("green", "red", "yellow", "blue"),
                         col.bg       = "black",
                         col.line     = "gray47",
                         has.legend   = TRUE)
    {
        dl                <- list()
        dl[["data"]]      <- data
        dl[["row"]]       <- levels(dl[["data"]][, 1])
        dl[["column"]]    <- levels(dl[["data"]][, 2])
        dl[["xchange"]]   <- c(-1, 0, -1,  0)
        dl[["ychange"]]   <- c( 0, 0, -1, -1)
        dl[["maxvalue"]]  <- max(dl[["data"]][, c(3:6)])
        dl[["colorname"]] <- col
        dl[[col[1]]]      <- colorRampPalette(c(col.bg, col[1]),
                                              bias = bias)(as.integer(dl[["maxvalue"]] * 10))
        dl[[col[2]]]      <- colorRampPalette(c(col.bg, col[2]),
                                              bias = bias)(as.integer(dl[["maxvalue"]] * 10))
        dl[[col[3]]]      <- colorRampPalette(c(col.bg, col[3]),
                                              bias = bias)(as.integer(dl[["maxvalue"]] * 10))
        dl[[col[4]]]      <- colorRampPalette(c(col.bg, col[4]),
                                              bias = bias)(as.integer(dl[["maxvalue"]] * 10))
        ##----draw a blank plot
        if (has.legend) {
            layout(mat = matrix(c(rep(1, 54), 2, rep(3:6, each = 2)), nrow = 9, ncol = 7))
        }
        par(mar    = c(0.1, 4.5, 4.1, 0.5))
        plot(x     = 0,
             y     = 0,
             type  = "n",
             xlab  = "",
             ylab  = "",
             xlim  = c(0, 2 * length(dl[["column"]])),
             ylim  = c(0, 2 * length(dl[["row"]])),
             axes  = FALSE)
        ##----draw heatmap
        column.order <- as.numeric(dl[["data"]][, 2])    # used for x 
        row.order    <- length(dl[["row"]]) + 1 - as.numeric(dl[["data"]][, 1])    # used for y
        for (i in c(1:4)) {
            ##----draw rectangle
            
            color.list   <- dl[[dl[["colorname"]][i]]][as.integer(dl[["data"]][, 2 + i] * 10)]
            rect(xleft   = 2 * column.order + dl[["xchange"]][i] - 1,
                 ybottom = 2 * row.order    + dl[["ychange"]][i] - 1,
                 xright  = 2 * column.order + dl[["xchange"]][i],
                 ytop    = 2 * row.order    + dl[["ychange"]][i],
                 border  = NA,
                 col     = color.list)
        }
        par(xpd = TRUE)
        for (i in c(0:length(dl[["row"]]))) {
            ##----draw horizontal lines
            lines(x   = c(-0.25, 2 * length(dl[["column"]])),
                  y   = 2 * c(i, i),
                  col = col.line,
                  lwd = 2)
        }
        for (i in c(0:length(dl[["column"]]))) {
            ##----draw vertical lines
            lines(x   = 2 * c(i, i),
                  y   = c(0, 2 * length(dl[["row"]]) + 0.25),
                  col = col.line,
                  lwd = 2)
        }
        ##----add labels
        text(y      = 2 * c(length(dl[["row"]]) : 1) - 1,
             x      = par("usr")[3] - 0.25,
             srt    = 0,
             labels = dl[["row"]],
             cex    = label.size)    # y
        mtext(text  = ylab,
              side  = 2,
              line  = 2,
              cex   = title.size)
        text(x      = 2 * c(1:length(dl[["column"]])) - 1,
             y      = 2 * length(dl[["row"]]) + 1.25,
             labels = dl[["column"]],
             cex    = label.size) # x
        mtext(text  = xlab,
              side  = 3,
              line  = 1,
              cex   = title.size)
        ##----add legends
        if (has.legend) {
            par(mar = c(0, 0, 0, 0))
            plot(x     = 0,
                 y     = 0,
                 type  = "n",
                 xlim  = c(0, 1),
                 ylim  = c(0, 1),
                 xlab  = "",
                 ylab  = "",
                 axes  = FALSE)
            par(mar = c(1.1, 0.3, 1.1, 5.1))
            for (i in c(1:length(dl[["colorname"]]))) {
                plot(x     = 0,
                     y     = 0,
                     type  = "n",
                     xlim  = c(0, length(dl[[dl[["colorname"]][i]]]) / 5),
                     ylim  = c(0, length(dl[[dl[["colorname"]][i]]])),
                     xlab  = "",
                     ylab  = "",
                     axes  = FALSE)
                rect(xleft   = rep(par("usr")[1], length(dl[[dl[["colorname"]][i]]])),
                     ybottom = c(0: (length(dl[[dl[["colorname"]][i]]]) - 1)),
                     xright  = rep(par("usr")[2], length(dl[[dl[["colorname"]][i]]])),
                     ytop    = c(1: length(dl[[dl[["colorname"]][i]]])),
                     border  = NA,
                     col     = dl[[dl[["colorname"]][i]]])
                axis(side    = 4,
                     at      = c(0,
                         (length(dl[[dl[["colorname"]][i]]])) / 2,
                         length(dl[[dl[["colorname"]][i]]])),
                     labels  = c(0, dl[["maxvalue"]] / 2, dl[["maxvalue"]]),
                     cex     = label.size * 2,
                     lwd     = label.size)
                mtext(text   = legend.label[i],
                      side   = 3)
            }
        }
    }





