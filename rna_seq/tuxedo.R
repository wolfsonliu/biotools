rm(list=ls())

################################################################################
##                      Define Functions and Class                            ##
################################################################################

## setClass(
##     Class = "RNASeq",
##     slots = c(
##         set       = "CuffSet",
##         condition = "character",
##         rep.num   = "numeric",
.sigrow <- function(data,
                   log2.threshold,
                   trend) {
    which(
        data$significant == "yes" &
        eval(
            call(
                ifelse(trend == "up", ">", "<"),
                data$log2_fold_change,
                log2.threshold
            )
        )
    )
}


################################################################################
##                Preparation: Library, Directory, File                       ##
################################################################################


###{{{ Library

suppressPackageStartupMessages(require(methods,    quietly = TRUE))
suppressPackageStartupMessages(require(reshape2,   quietly = TRUE))
suppressPackageStartupMessages(require(ggplot2,    quietly = TRUE))
suppressPackageStartupMessages(require(cummeRbund, quietly = TRUE))
suppressPackageStartupMessages(require(biomaRt,    quietly = TRUE))
suppressPackageStartupMessages(require(goseq,      quietly = TRUE))
source("keggenrichment.R")
source("rnafunction.R")

###}}}

###{{{ Working Directory and Variables

rna.dir                 <- list()         # store dirs.
rna.dir[["basic"]]      <- "/home/wolf/Documents/Work/Project/Yangyl/New" 
rna.dir[["cuffdiff"]]   <- file.path(rna.dir["basic"], "cuffdiff")
rna.dir[["info"]]       <- file.path(rna.dir["basic"], "info")
rna.dir[["cummeRbund"]] <- file.path(rna.dir["basic"], "cummeRbund")
setwd(rna.dir[["basic"]])                 # set working directory to basic.

## create dir if not exists.
invisible(
    lapply(
        rna.dir,
        function(x)
        {
            if (!dir.exists(x)) {
                dir.create(x)
            }
        }
    )
)

rna.file.path <- list()                 # store file paths.
rna.file.path[["tophat.result"]] <- file.path(rna.dir["info"], "tophat.csv")


###}}}

###{{{ Loading Data

rna.data <- list()

for (nam.data.dir in dir(rna.dir[["cuffdiff"]])) {
    suppressMessages(
        rna.data[[nam.data.dir]] <- readCufflinks(
            dir = file.path(
                rna.dir[["cuffdiff"]],
                nam.data.dir
            )
        )
    )
}

rna.info <- list()

rna.info[["tophat.result"]] <- read.csv(
    rna.file.path[["tophat.result"]],
    stringsAsFactors = FALSE
)

biomart <- list()
biomart[["ensembl_human"]] <- useMart(
    "ensembl",
    dataset="hsapiens_gene_ensembl"
)

###}}}

################################################################################
##                        Analysis and Graphics                               ##
################################################################################


###{{{ Using cummeRbund

rna.sig.geneid  <- list()
rna.test.result <- list()
rna.sig.gene    <- list()
rna.ggplot      <- list()
rna.de.geneid   <- list()


for (nam.data in names(rna.data)) {
    ## tmp.dir.path <- file.path(rna.dir[["cummeRbund"]], nam.data) 
    ## if (!dir.exists(tmp.dir.path)) {
    ##     ## Create directory to store files.
    ##     dir.create(tmp.dir.path)
    ## }
    ## ## Differential gene analysis
    ## rna.sig.geneid[[nam.data]]  <- getSig(
    ##     rna.data[[nam.data]],
    ##     alpha = 0.05,
    ##     level = "genes"
    ## )
    rna.de.geneid[[nam.data]]          <- list()
    rna.de.geneid[[nam.data]][["all"]] <- getBM(
        attributes = c("hgnc_symbol", "entrezgene"),
        filters    = "hgnc_symbol",
        values     = rna.sig.geneid[[nam.data]],
        mart       = biomart[["ensembl_human"]]
    )
    ## write(
    ##     ## write significant gene id to file.
    ##     rna.sig.geneid[[nam.data]],
    ##     file = file.path(
    ##         tmp.dir.path,
    ##         paste0(nam.data, "_siggeneid.txt")
    ##     )
    ## )    
    ## rna.test.result[[nam.data]] <- diffData(genes(rna.data[[nam.data]]))
    rna.de.geneid[[nam.data]][["up.0"]] <- getBM(
        attributes = c("hgnc_symbol", "entrezgene"),
        filters    = "hgnc_symbol",
        values     = rna.test.result[[nam.data]][
            .sigrow(rna.test.result[[nam.data]], 0, "up"), "gene_id"],
        mart       = biomart[["ensembl_human"]]
    )
    ## write(
    ##     ## write significant up gene id to file.
    ##     rna.test.result[[nam.data]][
    ##         .sigrow(rna.test.result[[nam.data]], 0, "up"),
    ##         "gene_id"],
    ##     file = file.path(
    ##         tmp.dir.path,
    ##         paste0(nam.data, "_siggeneid_up0.txt")
    ##     )
    ## )
    rna.de.geneid[[nam.data]][["up.1"]] <- getBM(
        attributes = c("hgnc_symbol", "entrezgene"),
        filters    = "hgnc_symbol",
        values     = rna.test.result[[nam.data]][
            .sigrow(rna.test.result[[nam.data]], 1, "up"), "gene_id"],
        mart       = biomart[["ensembl_human"]]
    )
    ## write(
    ##     ## write significant up with log2_fold_change > 1 gene id to file.
    ##     rna.test.result[[nam.data]][
    ##         .sigrow(rna.test.result[[nam.data]], 1, "up"),
    ##         "gene_id"],
    ##     file = file.path(
    ##         tmp.dir.path,
    ##         paste0(nam.data, "_siggeneid_up1.txt")
    ##     )
    ## )
    rna.de.geneid[[nam.data]][["down.0"]] <- getBM(
        attributes = c("hgnc_symbol", "entrezgene"),
        filters    = "hgnc_symbol",
        values     = rna.test.result[[nam.data]][
            .sigrow(rna.test.result[[nam.data]], 0, "down"), "gene_id"],
        mart       = biomart[["ensembl_human"]]
    )
    ## write(
    ##     ## write significant down gene id to file.
    ##     rna.test.result[[nam.data]][
    ##         .sigrow(rna.test.result[[nam.data]], 0, "down"),
    ##         "gene_id"],
    ##     file = file.path(
    ##         tmp.dir.path,
    ##         paste0(nam.data, "_siggeneid_down0.txt")
    ##     )
    ## )
    rna.de.geneid[[nam.data]][["down.1"]] <- getBM(
        attributes = c("hgnc_symbol", "entrezgene"),
        filters    = "hgnc_symbol",
        values     = rna.test.result[[nam.data]][
            .sigrow(rna.test.result[[nam.data]], -1, "down"), "gene_id"],
        mart       = biomart[["ensembl_human"]]
    )
    ## write(
    ##     ## write significant down with log2_fold_change < -1 gene id to file.
    ##     rna.test.result[[nam.data]][
    ##         .sigrow(rna.test.result[[nam.data]], -1, "down"),
    ##         "gene_id"],
    ##     file = file.path(
    ##         tmp.dir.path,
    ##         paste0(nam.data, "_siggeneid_down1.txt")
    ##     )
    ## )       
    ## rna.sig.gene                <- getGenes(
    ##     rna.data[[nam.data]],
    ##     rna.sig.geneid[[nam.data]]
    ## )
    ## ## Plot
    ## rna.ggplot[[nam.data]] <- list()
    ## pdf(
    ##     file.path(
    ##         tmp.dir.path,
    ##         paste0(nam.data, "_dendro.pdf")
    ##     )
    ## )
    ## assign(nam.data, rna.data[[nam.data]])
    ## csDendro(
    ##     genes(get(nam.data)),
    ##     replicates = TRUE
    ## )
    ## ## rm(list = nam.data)
    ## dev.off()
    ## rna.ggplot[[nam.data]][["boxplot"]] <- csBoxplot(
    ##     genes(rna.data[[nam.data]]),
    ##     replicates = TRUE
    ## ) + theme_bw()
    ## rna.ggplot[[nam.data]][["density"]] <- csDensity(
    ##     genes(rna.data[[nam.data]]),
    ##     replicates = TRUE
    ## ) + theme_bw()
    ## rna.ggplot[[nam.data]][["pairs"]]   <- csScatterMatrix(
    ##     genes(rna.data[[nam.data]]),
    ##     replicates = TRUE
    ## )
    ## rna.ggplot[[nam.data]][["scatter"]] <- csScatter(
    ##     genes(rna.data[[nam.data]]),
    ##     samples(rna.data[[nam.data]])$sample_name[1],
    ##     samples(rna.data[[nam.data]])$sample_name[2],
    ##     smooth = TRUE
    ## ) + theme_bw()
    ## rna.ggplot[[nam.data]][["maplot"]]  <- MAplot(
    ##     genes(rna.data[[nam.data]]),
    ##     samples(rna.data[[nam.data]])$sample_name[1],
    ##     samples(rna.data[[nam.data]])$sample_name[2]
    ## ) + theme_bw()
    ## rna.ggplot[[nam.data]][["volcano"]] <- csVolcano(
    ##     genes(rna.data[[nam.data]]),
    ##     samples(rna.data[[nam.data]])$sample_name[1],
    ##     samples(rna.data[[nam.data]])$sample_name[2],
    ##     alpha           = 0.05,
    ##     showSignificant = TRUE
    ## )
    ## for (nam.fig in names(rna.ggplot[[nam.data]])) {
    ##     ggsave(
    ##         file.path(
    ##             tmp.dir.path,
    ##             paste0(nam.data, "_", nam.fig, ".pdf")
    ##         ),
    ##         plot = rna.ggplot[[nam.data]][[nam.fig]]
    ##     )
    ## }
}

###}}}

###{{{ KEGG

rna.kegg <- list()

for (nam.de in names(rna.de.geneid)) {
    rna.kegg[[nam.de]] <- list()
    for (nam.situ in names(rna.de.geneid[[nam.de]])) {
        rna.kegg[[nam.de]][[nam.situ]] <- pathwayEnrichment(
            na.omit(rna.de.geneid[[nam.de]][[nam.situ]]$entrezgene),
            pathway.list      = grep(
                "hsa05",
                pathwayList("hsa")$ID,
                invert = TRUE,
                value  = TRUE),
            whole.gene.number = dim(rna.test.result[[nam.de]])[1]
        )
        tmp.order <- order(rna.kegg[[nam.de]][[nam.situ]]$adjusted.P.value)
        tmp.data  <- data.frame(
            pathway     = factor(
                as.character(rna.kegg[[nam.de]][[nam.situ]]$Term),
                levels = as.character(
                    rna.kegg[[nam.de]][[nam.situ]]$Term)[tmp.order],
            ),          
            log10.adj.p = -log10(rna.kegg[[nam.de]][[nam.situ]]$adjusted.P.value),
            order       = tmp.order
        )        
        p         <- termPlot(
            data    = tmp.data,
            term    = pathway,
            value   = log10.adj.p,
            fill    = order,
            title   = paste(
                nam.de,
                nam.situ,
                "KEGG pathway",
                sep = " "
            ),
            xlab    = "Pathway",
            ylab    = "",
            palette = ifelse(grepl("up", nam.situ), "Reds", "Blues")
        )
        ggsave(
            filename = file.path(
                rna.dir[["cummeRbund"]],
                nam.de,
                paste(nam.de, nam.situ, "kegg.pdf", sep = "_")
            ),
            plot = p
        )
        rm(tmp.order, tmp.data)
    }
}




###}}}

###{{{ GO: goseq

rna.go <- list()

for (nam.data in names(rna.de.geneid)) {
    for (nam.situ in names(rna.de.geneid[[nam.data]])) {
        tmp.genes <- as.integer(
            rna.test.result[[nam.data]]$gene_id %in%
            rna.de.geneid[[nam.data]][[nam.situ]]$hgnc_symbol
        )
        names(tmp.genes) <- rna.test.result[[nam.data]]$gene_id
        tmp.pwf <- nullp(tmp.genes,"hg19","geneSymbol")
        rna.go[[nam.data]][[nam.situ]] <- goseq(
            tmp.pwf,
            "hg19",
            "geneSymbol",
            test.cats=c("GO:BP")
        )
        rna.go[[nam.data]][[nam.situ]]$adjusted.P.value <- p.adjust(
            rna.go[[nam.data]][[nam.situ]]$over_represented_pvalue
        )
        tmp.go <- rna.go[[nam.data]][[nam.situ]][
            which(
                rna.go[[nam.data]][[nam.situ]]$adjusted.P.value < 0.05
            ),]
        if (dim(tmp.go)[1] > 20) {
            tmp.go <- tmp.go[order(tmp.go$adjusted.P.value)[1:20],]
        }
        tmp.order <- order(tmp.go$adjusted.P.value)
        tmp.data  <- data.frame(
            go          = factor(
                tmp.go$term,
                levels = as.character(tmp.go$term)[tmp.order],
            ),          
            log10.adj.p = -log10(tmp.go$adjusted.P.value),
            order       = tmp.order
        )                    
        p         <- termPlot(
            data    = tmp.data,
            term    = go,
            value   = log10.adj.p,
            fill    = order,
            title   = paste(
                nam.data,
                nam.situ,
                "GO BP",
                sep = " "
            ),
            xlab    = "GO Term",
            ylab    = "",
            palette = ifelse(grepl("up", nam.situ), "Reds", "Blues")
        )
        ggsave(
            filename = file.path(
                rna.dir[["cummeRbund"]],
                nam.data,
                paste(nam.data, nam.situ, "go.pdf", sep = "_")
            ),
            plot = p
        )
        write.table(
            tmp.go,
            file = file.path(
                rna.dir[["cummeRbund"]],
                nam.data,
                paste(nam.data, nam.situ, "go.csv", sep = "_")
            ),
            sep = ",",
            row.names = FALSE
        )
        rm(tmp.order, tmp.data, tmp.go, p)
    }
}


###}}}


###{{{ Infomation


## Mapping Result
for (nam.col in c("Left.map",
                  "Left.multimap",
                  "Right.map",
                  "Right.multimap",
                  "Pairs.map",
                  "Pairs.multimap")) {
    ## Calculate the ratio of mapped reads or pairs.
    rna.info[["tophat.result"]][[paste0(
                                    nam.col,
                                    "ratio"
                                )]] <- rna.info[["tophat.result"]][[nam.col]] /
        rna.info[["tophat.result"]][["Input"]]
}

p.map.data <- melt(
    ## Reshape data for plot.
    rna.info[["tophat.result"]][,c(1,9:14)]
)

p.map.data <- rbind(
    ## Calculate the mean of ratios.
    p.map.data,
    data.frame(
        Sample   = rep("mean", 6),
        variable = c(
            "Left.mapratio",
            "Left.multimapratio",
            "Right.mapratio",
            "Right.multimapratio",
            "Pairs.mapratio",
            "Pairs.multimapratio"),
        value    = as.numeric(
            addmargins(
                as.matrix(rna.info[["tophat.result"]][,9:14]),
                1,
                FUN = mean)[25,]
        )
    )
)

p.map <- plotDodgeBarplot(
    p.map.data
)

ggsave(
    filename = file.path(
        rna.dir[["info"]],
        paste0("tophat.pdf")
    ),
    p.map,
    width  = 20,
    height = 10,
    units  = "cm"
)

rm(
    p.map.data,
    p.map
)

## Differential Expression Result

p.diff.data <- melt(
    data.frame(
        Sample        = names(rna.sig.geneid),
        Sig.gene      = as.numeric(lapply(rna.sig.geneid, length)),
        Sig.up.gene   = as.numeric(
            lapply(
                rna.test.result,
                function(x) {
                    length(as.character(x[.sigrow(x, 0, "up"), "gene_id"]))
                }
            )
        ),
        Sig.up.gene.1 = as.numeric(
            lapply(
                rna.test.result,
                function(x) {
                    length(as.character(x[.sigrow(x, 1, "up"), "gene_id"]))
                }
            )
        ),
        Sig.down.gene   = as.numeric(
            lapply(
                rna.test.result,
                function(x) {
                    length(as.character(x[.sigrow(x, 0, "down"), "gene_id"]))
                }
            )
        ),
        Sig.down.gene.1 = as.numeric(
            lapply(
                rna.test.result,
                function(x) {
                    length(as.character(x[.sigrow(x, -1, "down"), "gene_id"]))
                }
            )
        )
    )
)

p.diff <- plotDodgeBarplot(p.diff.data)

ggsave(
    filename = file.path(
        rna.dir[["info"]],
        paste0("de_gene.pdf")
    ),
    p.diff,
    width  = 20,
    height = 10,
    units  = "cm"
)

rm(
    p.diff.data,
    p.diff
)
            
    
###}}}

################################################################################
##                               SCRATCH                                      ##
################################################################################

###{{{ Scratch

tmp <- getBM(
    attributes = c("hgnc_symbol", "entrezgene"),
    filters    = "hgnc_symbol",
    values     = rna.sig.geneid[["jeko_0_6"]],
    mart       = biomart[["ensembl_human"]]
)

tmp.pathway <- pathwayEnrichment(na.omit(tmp$entrezgene), whole.gene.number = 20362)


tmp               <- readCufflinks(
    dir = file.path(
        rna.dir[["cuffdiff"]],
        "jeko_0_6"
    )
)

tmp.sigGeneIds    <- getSig(tmp, alpha=0.05,level="genes")

tmp.gene.diff     <- diffData(genes(tmp))

tmp.myGenes       <- getGenes(tmp, tmp.sigGeneIds)

tmp.dispersion    <- dispersionPlot(genes(tmp))

tmp.pBoxRep       <- csBoxplot(genes(tmp), replicates=T)

tmp.pDendro       <- csDendro(genes(tmp), replicates=T)

tmp.genes.scv     <- fpkmSCVPlot(genes(tmp))

tmp.densRep       <- csDensity(genes(tmp),replicates=T)

tmp.scattermatrix <- csScatterMatrix(genes(tmp),replicates = T)

tmp.scatter       <- csScatter(genes(tmp), "H0", "H6", smooth=T)

tmp.ma            <- MAplot(genes(tmp),"H0","H6") +
    theme_bw()

tmp.volcanomatrix <- csVolcanoMatrix(genes(tmp))

tmp.volcano       <- csVolcano(
    genes(tmp),
    samples(tmp)$sample_name[1],
    samples(tmp)$sample_name[2],
    alpha           = 0.05,
    showSignificant = TRUE
)

tmp.heatmap       <-csHeatmap(tmp.myGenes,cluster="both",replicates = TRUE)

rm(list = grep("tmp", ls(), value = TRUE))



myfun <- function(){
    for (h in 1:50) {
        Sys.sleep(1)
        cat("Fetching pathway information", sep = "")
        cat(" * ", sep = "")
        cat(switch(h%%3 + 1, "/", "\\", "-"), sep = "")
        cat(" * \r", sep = "")
    }
    cat("\n")
}
###}}}
