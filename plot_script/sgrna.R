#! /bin/env R

################################################################################
## Author:      Wolfson                                                       ##
## Date:        2016-02-29                                                    ##
## Modified:    2016-02-29                                                    ##
## Description: Plot sgRNA screen data                                        ##
################################################################################


symbols(x       = mydata1$id,
        y       = -log2(mydata1$pos.score),
        circles = mydata1$pos.goodsgrna,
        inches  = 1/5,
        add     = FALSE,
        fg      = "gray70",
        bg      = "gray70",
        xaxs    = "i",
        yaxs    = "i",
        xlab    = "gene id",
        ylab    = "significance (-log2(pos score))")
symbols(x       = mydata1$id[mydata1$pos.fdr < 0.25],
        y       = -log2(mydata1$pos.score[mydata1$pos.fdr < 0.25]),
        circles = mydata1$pos.goodsgrna[mydata1$pos.fdr < 0.25],
        inches  = 1/5,
        add     = TRUE,
        fg      = "gray20",
        bg      = 2,
        xaxs    = "i",
        yaxs    = "i",
        xaxt    = "n",
        yaxt    = "n")
text(x          = mydata1$id[mydata1$pos.fdr < 0.25],
     y          = -log2(mydata1$pos.score[mydata1$pos.fdr < 0.25]),
     labels     = paste0(mydata1$id[mydata1$pos.fdr < 0.25], "(", mydata1$pos.goodsgrna[mydata1$pos.fdr < 0.25], ")"),
     pos        = 4,
     cex        = 0.6)
