################################################################################


################################################################################
##                           Required Library                                 ##
################################################################################
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(reshape2))

################################################################################
##                              Functions                                     ##
################################################################################

pairsPanelAblineRug <- function(x,
                                y,
                                col.point  = 2,
                                col.abline = 1,
                                col.rug    = 1,
                                ...) {
###{{{
    points(
        x, y, 
        pch = 20, 
        col = col.point,
        ...)
    abline(
        a   = 0,
        b   = 1,
        lwd = 2,
        col = col.abline,
        ...)
    rug(x, side = 1, col = col.rug)
    rug(y, side = 2, col = col.rug)
###}}}
}

pairsPanelCor <- function(x,
                          y,
                          method    = "pearson",
                          words     = "",
                          words.cex = 2,
                          ...) {
###{{{
    text(
        x      = mean(par()$usr[1:2]),
        y      = mean(par()$usr[3:4]),
        labels = paste(
            words, "\n",
            format(cor(x, y, method = method), digits = 2),
            sep = ""),
        cex    = words.cex,
        ...)
###}}}
}


termPlot <- function(data,
                     term,
                     value,
                     fill,
                     title   = "",
                     xlab    = "",
                     ylab    = "",
                     palette = c("Reds", "Blues")) {
###{{{

    formal <- as.list(match.call()[-1])
    p <- ggplot(
        get(x = as.character(formal$data)),
        aes_string(
            x    = as.character(formal$term),
            y    = as.character(formal$value),
            fill = as.character(formal$fill))) +
        geom_bar(stat = "identity") + theme_classic() +
        theme(
            legend.position  = "none",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.ticks       = element_blank()) +
        xlab(xlab) + ylab(ylab) + ggtitle(title) + 
        coord_flip() +
        scale_fill_distiller(palette = I(palette[1]))
    p

###}}}
}

clustChangeFactorOrder <- function(data, 
                                   factor) {
###{{{

    ## data is used to clust
    ## factor is usde to change
    clust.object <- hclust(dist(data))
    factor(
        factor,
        clust.object$labels[clust.object$order],
        ordered = TRUE)        

###}}}
}

heatMapPlot <- function(data,
                        xlab.size  = 15,
                        xlab.angle = 45,
                        xlab.vjust = 0.4,
                        xlab.color = "black",
                        ylab.size  = 15,
                        ylab.angle = 0,
                        ylab.hjust = 0,
                        ylab.vjust = 0,
                        ylab.color = "black",
                        xcluster   = FALSE,
                        ycluster   = FALSE,
                        xorder     = "",
                        yorder     = "",
                        color      = "") {
###{{{

    if (!is.matrix(data)) {
        ## Input data should be a matrix with rownames.
        stop("Please enter a numeric matrix with rownames as input data.")
    }
    melt.data <- reshape2::melt(data)
    if (sum(nchar(xorder) != 0) != 0 & xcluster) {
        stop("Please choose either using cluster or using order, not both.")
    } else {}
    if (sum(nchar(yorder) != 0) != 0 & ycluster) {
        stop("Please choose either using cluster or using order, not both.")
    } else {}
    if (xcluster) {
        melt.data$Var2 <- clustChangeFactorOrder(data, melt.data$Var2)
    } else {}
    if (ycluster) {
        melt.data$Var1 <- clustChangeFactorOrder(data, melt.data$Var1)
    } else {}
    if (sum(nchar(xorder) != 0) != 0) {
        if (length(xorder) != length(unique(as.character((melt.data$Var2))))) {
            stop("Order of x not the same lenght with level of x.")
        } else {}    
        melt.data$Var2 <- factor(
            melt.data$Var2,
            xorder,
            ordered = TRUE)
    } else {}
    if (sum(nchar(yorder) != 0) != 0) {
        if (length(yorder) != length(unique(as.character((melt.data$Var1))))) {
            stop("Order of y not the same lenght with level of y.")
        } else {}    
        melt.data$Var1 <- factor(
            melt.data$Var1,
            yorder,
            ordered = TRUE)
    } else {}
    if (length(color) == 1 & sum(nchar(color)) == 0) {
        ## Setting color palette.
        palette.setting <- scale_fill_gradient(
            low  = "white",
            high = "red")
    } else if (length(color) == 1 & sum(nchar(color)) != 0) {
        palette.setting <- scale_fill_gradient(
            low  = "white",
            high = color[1])
    } else if (length(color) == 2) {
        palette.setting <- scale_fill_gradient(
            low  = color[1],
            high = color[2])
    } else if (length(color) == 3) {
        palette.setting <- scale_fill_gradient2(
            low  = color[1],
            mid  = color[2],
            high = color[3])
    }
    p <- ggplot(
        melt.data,
        aes(
            x    = Var2,
            y    = Var1,
            fill = value)) +
        geom_tile(colour = "white") + palette.setting +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.ticks       = element_blank(),
              axis.text.x      = element_text(
                  size   = xlab.size,
                  angle  = xlab.angle,
                  vjust  = xlab.vjust,
                  colour = xlab.color),
              axis.text.y      = element_text(
                  hjust  = ylab.hjust,
                  size   = ylab.size,
                  angle  = ylab.angle,
                  vjust  = ylab.vjust,
                  colour = ylab.color),
              axis.title       = element_blank(),
              plot.title       = element_blank())
    p
    
###}}}
}


plotDodgeBarplot <- function(data,
                             xlab.size  = 10,
                             xlab.angle = 45,
                             xlab.vjust = 0.4) {
###{{{

    p.data <- melt(data)
    p <- ggplot(
        p.data,
        aes(x     = Sample,
            y     = value,
            group = variable,
            fill  = variable
            )
    ) +
        geom_bar(stat = "identity", position = "dodge") +
        theme_bw() +
        theme(
            axis.text.x      = element_text(
                size   = xlab.size,
                angle  = xlab.angle,
                vjust  = xlab.vjust
            )
        )

###}}}
}




                             
################################################################################
##                                 EOF                                        ##
################################################################################
