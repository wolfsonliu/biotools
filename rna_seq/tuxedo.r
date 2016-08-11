rm(list=ls())

################################################################################
##                      Define Functions and Class                            ##
################################################################################

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

###{{{ Working Directory and Variables

## Directories
rna.dir                 <- list()         # store dirs.
rna.dir[["basic"]]      <- "/home/wolf/Documents/Work/Project/Yangyl/New" 
rna.dir[["cuffdiff"]]   <- file.path(rna.dir[["basic"]], "cuffdiff")
rna.dir[["info"]]       <- file.path(rna.dir[["basic"]], "info")
rna.dir[["cummeRbund"]] <- file.path(rna.dir[["basic"]], "cummeRbund")
rna.dir[["circos"]]     <- file.path(rna.dir[["basic"]], "circos")
rna.dir[["htseq"]]      <- file.path(rna.dir[["basic"]], "htseq_count")
rna.dir[["genes"]]      <- file.path(rna.dir[["basic"]], "genes")
setwd(rna.dir[["basic"]])                 # set working directory to basic.

## Create directories if not exists.
invisible(
    lapply(
        rna.dir,
        function(x)
        {
            if (!dir.exists(x)) {
                dir.create(x)
            } else {}
        }
    )
)

## File paths
rna.file.path <- list()                 # store file paths.
rna.file.path[["tophat.result"]] <- file.path(
    rna.dir[["info"]],
    "tophat.csv"
)
rna.file.path[["gsea.tft"]]       <- file.path(
    rna.dir[["info"]],
    "msigdb/c3.tft.v5.1.symbols.gmt"
)
rna.file.path[["gsea"]] <- file.path(
    rna.dir[["info"]],
    "msigdb/msigdb_v5.1.xml"
)
rna.file.path[["aa.gene"]] <- file.path(
    rna.dir[["info"]],
    "aa.gene.csv"
)

###}}}

###{{{ Library

suppressPackageStartupMessages(require(methods,        quietly = TRUE))
suppressPackageStartupMessages(require(reshape2,       quietly = TRUE))
suppressPackageStartupMessages(require(ggplot2,        quietly = TRUE))
suppressPackageStartupMessages(require(gplots,         quietly = TRUE))
suppressPackageStartupMessages(require(cummeRbund,     quietly = TRUE))
suppressPackageStartupMessages(require(biomaRt,        quietly = TRUE))
suppressPackageStartupMessages(require(goseq,          quietly = TRUE))
suppressPackageStartupMessages(require(preprocessCore, quietly = TRUE))
suppressPackageStartupMessages(require(XML,            quietly = TRUE))
suppressPackageStartupMessages(require(RColorBrewer,   quietly = TRUE))
suppressPackageStartupMessages(require(VennDiagram,    quietly = TRUE))
source(file.path(dirname(rna.dir[["basic"]]), "keggenrichment.r"))
source(file.path(dirname(rna.dir[["basic"]]), "rnafunction.r"))
source(file.path(dirname(rna.dir[["basic"]]), "mycircos.r"))

###}}}

###{{{ Loading Data


## Loading Main Data
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

rna.count <- lapply(
    dir(rna.dir[["htseq"]]),
    function(x) {
        tmp <- read.table(
            file             = file.path(rna.dir[["htseq"]], x, "count.txt"),
            header           = FALSE,
            stringsAsFactors = FALSE
        )
        colnames(tmp) <- c("gene", "count")
        rownames(tmp) <- tmp$gene
        tmp[grep("__", tmp$gene, invert = TRUE),]
    }
)

names(rna.count) <- dir(rna.dir[["htseq"]])
                   
    


            
## Loading Information Data
rna.info <- list()                      # Store information data and public data.

rna.info[["tophat.result"]] <- read.csv(
    rna.file.path[["tophat.result"]],
    stringsAsFactors = FALSE
)

rna.info[["gsea.tft"]] <- lapply(
    strsplit(
        readLines(
            rna.file.path[["gsea.tft"]],
            n = -1
        ),
        split = "\t"
    ),
    function(x) {
        list(
            tft  = x[1],
            web  = x[2],
            gene = x[-c(1,2)]
        )
    }
)

names(rna.info[["gsea.tft"]]) <- unlist(
    lapply(
        rna.info[["gsea.tft"]],
        "[[",
        "tft"
    )
)

msigdb <- xmlToList(xmlParse(rna.file.path[["gsea"]]))


rna.info[["aa.gene"]] <- read.csv(
    rna.file.path[["aa.gene"]],
    header           = TRUE,
    sep              = ",",
    stringsAsFactors = FALSE
)

## Loading Package Data
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
    tmp.dir.path <- file.path(rna.dir[["cummeRbund"]], nam.data) 
    if (!dir.exists(tmp.dir.path)) {
        ## Create directory to store files.
        dir.create(tmp.dir.path)
    } else {}
    ## Differential gene analysis
    rna.sig.geneid[[nam.data]]  <- getSig(
        rna.data[[nam.data]],
        alpha = 0.05,
        level = "genes"
    )
    rna.de.geneid[[nam.data]]          <- list()
    rna.de.geneid[[nam.data]][["all"]] <- getBM(
        attributes = c("hgnc_symbol", "entrezgene"),
        filters    = "hgnc_symbol",
        values     = rna.sig.geneid[[nam.data]],
        mart       = biomart[["ensembl_human"]]
    )
    write(
        ## write significant gene id to file.
        rna.sig.geneid[[nam.data]],
        file = file.path(
            tmp.dir.path,
            paste0(nam.data, "_siggeneid.txt")
        )
    )    
    rna.test.result[[nam.data]] <- diffData(genes(rna.data[[nam.data]]))
    rna.de.geneid[[nam.data]][["up.0"]] <- getBM(
        attributes = c("hgnc_symbol", "entrezgene"),
        filters    = "hgnc_symbol",
        values     = rna.test.result[[nam.data]][
            .sigrow(rna.test.result[[nam.data]], 0, "up"), "gene_id"],
        mart       = biomart[["ensembl_human"]]
    )
    write(
        ## write significant up gene id to file.
        rna.test.result[[nam.data]][
            .sigrow(rna.test.result[[nam.data]], 0, "up"),
            "gene_id"],
        file = file.path(
            tmp.dir.path,
            paste0(nam.data, "_siggeneid_up0.txt")
        )
    )
    rna.de.geneid[[nam.data]][["up.1"]] <- getBM(
        attributes = c("hgnc_symbol", "entrezgene"),
        filters    = "hgnc_symbol",
        values     = rna.test.result[[nam.data]][
            .sigrow(rna.test.result[[nam.data]], 1, "up"), "gene_id"],
        mart       = biomart[["ensembl_human"]]
    )
    write(
        ## write significant up with log2_fold_change > 1 gene id to file.
        rna.test.result[[nam.data]][
            .sigrow(rna.test.result[[nam.data]], 1, "up"),
            "gene_id"],
        file = file.path(
            tmp.dir.path,
            paste0(nam.data, "_siggeneid_up1.txt")
        )
    )
    rna.de.geneid[[nam.data]][["down.0"]] <- getBM(
        attributes = c("hgnc_symbol", "entrezgene"),
        filters    = "hgnc_symbol",
        values     = rna.test.result[[nam.data]][
            .sigrow(rna.test.result[[nam.data]], 0, "down"), "gene_id"],
        mart       = biomart[["ensembl_human"]]
    )
    write(
        ## write significant down gene id to file.
        rna.test.result[[nam.data]][
            .sigrow(rna.test.result[[nam.data]], 0, "down"),
            "gene_id"],
        file = file.path(
            tmp.dir.path,
            paste0(nam.data, "_siggeneid_down0.txt")
        )
    )
    rna.de.geneid[[nam.data]][["down.1"]] <- getBM(
        attributes = c("hgnc_symbol", "entrezgene"),
        filters    = "hgnc_symbol",
        values     = rna.test.result[[nam.data]][
            .sigrow(rna.test.result[[nam.data]], -1, "down"), "gene_id"],
        mart       = biomart[["ensembl_human"]]
    )
    write(
        ## write significant down with log2_fold_change < -1 gene id to file.
        rna.test.result[[nam.data]][
            .sigrow(rna.test.result[[nam.data]], -1, "down"),
            "gene_id"],
        file = file.path(
            tmp.dir.path,
            paste0(nam.data, "_siggeneid_down1.txt")
        )
    )       
    rna.sig.gene                <- getGenes(
        rna.data[[nam.data]],
        rna.sig.geneid[[nam.data]]
    )
    ## Plot
    rna.ggplot[[nam.data]] <- list()
    pdf(
        file.path(
            tmp.dir.path,
            paste0(nam.data, "_dendro.pdf")
        )
    )
    assign(nam.data, rna.data[[nam.data]])
    csDendro(
        genes(get(nam.data)),
        replicates = TRUE
    )
    ## rm(list = nam.data)
    dev.off()
    rna.ggplot[[nam.data]][["boxplot"]] <- csBoxplot(
        genes(rna.data[[nam.data]]),
        replicates = TRUE
    ) + theme_bw()
    rna.ggplot[[nam.data]][["density"]] <- csDensity(
        genes(rna.data[[nam.data]]),
        replicates = TRUE
    ) + theme_bw()
    rna.ggplot[[nam.data]][["pairs"]]   <- csScatterMatrix(
        genes(rna.data[[nam.data]]),
        replicates = TRUE
    )
    rna.ggplot[[nam.data]][["scatter"]] <- csScatter(
        genes(rna.data[[nam.data]]),
        samples(rna.data[[nam.data]])$sample_name[1],
        samples(rna.data[[nam.data]])$sample_name[2],
        smooth = TRUE
    ) + theme_bw()
    rna.ggplot[[nam.data]][["maplot"]]  <- MAplot(
        genes(rna.data[[nam.data]]),
        samples(rna.data[[nam.data]])$sample_name[1],
        samples(rna.data[[nam.data]])$sample_name[2]
    ) + theme_bw()
    rna.ggplot[[nam.data]][["volcano"]] <- csVolcano(
        genes(rna.data[[nam.data]]),
        samples(rna.data[[nam.data]])$sample_name[1],
        samples(rna.data[[nam.data]])$sample_name[2],
        alpha           = 0.05,
        showSignificant = TRUE
    )
    tmp.data <- rna.test.result[[nam.data]]
    tmp.data$DE <- c("gray","red")[(tmp.data$significant == "yes") + 1]
    tmp.p <- ggplot(
        tmp.data,
        aes(
            x      = value_1 + 1,
            y      = value_2 + 1,
            colour = I(DE)
        )
    ) + geom_point() + geom_rug() +
        geom_abline(slope = 1, intercept = 0, colour = "blue", size = 1.5) +
        scale_x_log10() + scale_y_log10() +
        ggtitle(paste(nam.data, "Differencial Expressed Genes")) +
        xlab(samples(rna.data[[nam.data]])$sample_name[1]) +
        ylab(samples(rna.data[[nam.data]])$sample_name[2]) +
        theme_bw()
    ggsave(
        file.path(
            tmp.dir.path,
            paste0(nam.data, "_descatter.pdf")
        ),
        plot = tmp.p
    )
    for (nam.fig in names(rna.ggplot[[nam.data]])) {
        ggsave(
            file.path(
                tmp.dir.path,
                paste0(nam.data, "_", nam.fig, ".pdf")
            ),
            plot = rna.ggplot[[nam.data]][[nam.fig]]
        )
    }
    rm(tmp.data, tmp.p)
}

###}}}

###{{{ KEGG

rna.kegg <- list()




for (nam.de in names(rna.de.geneid)[1]) {
    rna.kegg[[nam.de]] <- list()
    for (nam.situ in names(rna.de.geneid[[nam.de]])[1]) {
        rna.kegg[[nam.de]][[nam.situ]] <- pathwayEnrichment(
            na.omit(rna.de.geneid[[nam.de]][[nam.situ]]$entrezgene),
            pathway.list      = grep(
                "hsa05",
                pathwayList("hsa")$ID,
                invert = TRUE,
                value  = TRUE),
            whole.gene.number = dim(rna.test.result[[nam.de]])[1],
            threshold         = 1,
            organism          = "hsa"
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
                paste(nam.de, nam.situ, "keggord.pdf", sep = "_")
            ),
            plot = p
        )
        rm(tmp.order, tmp.data)
    }
}


tmp.kegg <- lapply(
    rna.kegg,
    function(x) {
        lapply(
            x,
            function(y) {
                na.omit(y[order(y$adjusted.P.value)[1:5],])
            }
        )
    }
)

tmp.pathwayorder <- list(
    jeko = c(
        "Amino sugar and nucleotide sugar metabolism",
        "Apoptosis",
        "B cell receptor signaling pathway",
        "Cell adhesion molecules (CAMs)",
        "Endocytosis",
        "ErbB signaling pathway",
        "Estrogen signaling pathway",
        "FoxO signaling pathway",
        "Glycosylphosphatidylinositol(GPI)-anchor biosynthesis",
        "HIF-1 signaling pathway",
        "Lysosome",
        "Inositol phosphate metabolism",
        "Insulin signaling pathway",
        "MAPK signaling pathway",
        "mTOR signaling pathway",
        "Neurotrophin signaling pathway",
        "N-Glycan biosynthesis",
        "Osteoclast differentiation",
        "p53 signaling pathway",
        "Phagosome",
        "Phosphatidylinositol signaling system",
        "Protein export",
        "Protein processing in endoplasmic reticulum",
        "SNARE interactions in vesicular transport",
        "Sphingolipid signaling pathway",
        "Toll-like receptor signaling pathway",
        "Ubiquitin mediated proteolysis",
        "VEGF signaling pathway",
        "Wnt signaling pathway",
        "Aminoacyl-tRNA biosynthesis",
        "Antigen processing and presentation",
        "Base excision repair",
        "Cell cycle",
        "DNA replication",
        "Fanconi anemia pathway",
        "Homologous recombination",
        "Intestinal immune network for IgA production",
        "Mismatch repair",
        "Non-alcoholic fatty liver disease (NAFLD)",
        "Nucleotide excision repair",
        "One carbon pool by folate",
        "Oxidative phosphorylation",
        "Proteasome",
        "Purine metabolism",
        "Pyrimidine metabolism",
        "Pyruvate metabolism",
        "Ribosome",
        "RNA transport",
        "Spliceosome",
        "Type I diabetes mellitus"
    ),
    raji = c(
        "Apoptosis",
        "B cell receptor signaling pathway",
        "Base excision repair",
        "Cell cycle",
        "Citrate cycle (TCA cycle)",
        "DNA replication",
        "Endocytosis",
        "Insulin signaling pathway",
        "Mismatch repair",
        "Non-alcoholic fatty liver disease (NAFLD)",
        "Non-homologous end-joining",
        "Nucleotide excision repair",
        "Oxidative phosphorylation",
        "Oocyte meiosis",
        "Proteasome",
        "Protein export",
        "Protein processing in endoplasmic reticulum",
        "Purine metabolism",
        "Pyrimidine metabolism",
        "Ribosome biogenesis in eukaryotes",
        "RNA degradation",
        "RNA transport",
        "Spliceosome",
        "Steroid biosynthesis",
        "Toll-like receptor signaling pathway",
        "Ubiquitin mediated proteolysis",
        "Valine, leucine and isoleucine degradation",
        "Antigen processing and presentation",
        "cGMP-PKG signaling pathway",
        "Cytokine-cytokine receptor interaction",
        "Fanconi anemia pathway",
        "Glycosaminoglycan biosynthesis - chondroitin sulfate / dermatan sulfate",
        "Glycosylphosphatidylinositol(GPI)-anchor biosynthesis",
        "Hedgehog signaling pathway",
        "Intestinal immune network for IgA production",
        "MAPK signaling pathway",
        "Natural killer cell mediated cytotoxicity",
        "NF-kappa B signaling pathway",
        "NOD-like receptor signaling pathway",
        "One carbon pool by folate",
        "Osteoclast differentiation",
        "Other glycan degradation",
        "p53 signaling pathway",
        "Pantothenate and CoA biosynthesis",
        "Peroxisome",
        "Phagosome",
        "Ribosome",
        "RIG-I-like receptor signaling pathway",
        "Selenocompound metabolism",
        "TNF signaling pathway",
        "Type I diabetes mellitus"
    )
)

tmp.pathwayorder2 <- list(
    jeko = c(
        "Apoptosis",
        "B cell receptor signaling pathway",
        "Endocytosis",
        "ErbB signaling pathway",
        "FoxO signaling pathway",
        "Lysosome",
        "Insulin signaling pathway",
        "Neurotrophin signaling pathway",
        "N-Glycan biosynthesis",
        "Phagosome",
        "Phosphatidylinositol signaling system",
        "Protein export",
        "Protein processing in endoplasmic reticulum",
        "SNARE interactions in vesicular transport",
        "Wnt signaling pathway",
        "Aminoacyl-tRNA biosynthesis",
        "Antigen processing and presentation",
        "Base excision repair",
        "DNA replication",
        "Intestinal immune network for IgA production",
        "Mismatch repair",
        "Oxidative phosphorylation",
        "Proteasome",
        "Pyrimidine metabolism",
        "Ribosome",
        "Spliceosome",
        "Type I diabetes mellitus"
    ),
    raji = c(
        "Apoptosis",
        "B cell receptor signaling pathway",
        "Base excision repair",
        "Cell cycle",
        "DNA replication",
        "Endocytosis",
        "Insulin signaling pathway",
        "Non-homologous end-joining",
        "Proteasome",
        "Protein export",
        "Protein processing in endoplasmic reticulum",
        "Pyrimidine metabolism",
        "RNA degradation",
        "RNA transport",
        "Spliceosome",
        "Steroid biosynthesis",
        "Antigen processing and presentation",
        "Fanconi anemia pathway",
        "Intestinal immune network for IgA production",
        "MAPK signaling pathway",
        "NF-kappa B signaling pathway",
        "Osteoclast differentiation",
        "p53 signaling pathway",
        "Ribosome",
        "RIG-I-like receptor signaling pathway",
        "TNF signaling pathway",
        "Type I diabetes mellitus"
    )
)




tmp.pathwaylist <- pathwayList("hsa")

for (cell.line in c("jeko", "raji")) {
    if (cell.line == "jeko") {
        tmp.order <- c(1, 2, 3, 4,  5)
    } else {
        tmp.order <- c(7, 8, 9, 10, 11)
    }
    tmp.pathwayname <- lapply(
        tmp.kegg[tmp.order],
        function(x) {
            lapply(
                x[c(2,4)],
                function(y) {
                    as.character(rownames(y))
                }
            )
        }
    )
    tmp.pathwayall <- unique(unlist(tmp.pathwayname))
    rownames(tmp.pathwaylist) <- tmp.pathwaylist$ID
    #    pathway = tmp.pathwaylist[tmp.pathwayall, "Name"]
    .adjPvalue <- function(x, ..., nam.data, nam.situ) {
        ifelse(
            is.na(rna.kegg[[nam.data]][[nam.situ]][x, "adjusted.P.value"]),
            1,
            rna.kegg[[nam.data]][[nam.situ]][x, "adjusted.P.value"]
        )
    }
    tmp.heatmatrix <- data.frame(
        H0_H6.down = log10(
            unlist(
                lapply(
                    tmp.pathwayall,
                    .adjPvalue,
                    nam.data = paste0(cell.line, "_0_6"),
                    nam.situ = "down.0"
                )
            )
        ),
        H0_H6.up = -log10(
            unlist(
                lapply(
                    tmp.pathwayall,
                    .adjPvalue,
                    nam.data = paste0(cell.line, "_0_6"),
                    nam.situ = "up.0"
                )
            )
        ),
        H6_H24.down = log10(
            unlist(
                lapply(
                    tmp.pathwayall,
                    .adjPvalue,
                    nam.data = paste0(cell.line, "_6_24"),
                    nam.situ = "down.0"
                )
            )
        ),
        H6_H24.up = -log10(
            unlist(
                lapply(
                    tmp.pathwayall,
                    .adjPvalue,
                    nam.data = paste0(cell.line, "_6_24"),
                    nam.situ = "up.0"
                )
            )
        ),
        H0_H24.down = log10(
            unlist(
                lapply(
                    tmp.pathwayall,
                    .adjPvalue,
                    nam.data = paste0(cell.line, "_0_24"),
                    nam.situ = "down.0"
                )
            )
        ),
        H0_H24.up = -log10(
            unlist(
                lapply(
                    tmp.pathwayall,
                    .adjPvalue,
                    nam.data = paste0(cell.line, "_0_24"),
                    nam.situ = "up.0"
                )
            )
        ),
        H24_H48.down = log10(
            unlist(
                lapply(
                    tmp.pathwayall,
                    .adjPvalue,
                    nam.data = paste0(cell.line, "_24_48"),
                    nam.situ = "down.0"
                )
            )
        ),
        H24_H48.up = -log10(
            unlist(
                lapply(
                    tmp.pathwayall,
                    .adjPvalue,
                    nam.data = paste0(cell.line, "_24_48"),
                    nam.situ = "up.0"
                )
            )
        ),
        H0_H48.down = log10(
            unlist(
                lapply(
                    tmp.pathwayall,
                    .adjPvalue,
                    nam.data = paste0(cell.line, "_0_48"),
                    nam.situ = "down.0"
                )
            )
        ),
        H0_H48.up = -log10(
            unlist(
                lapply(
                    tmp.pathwayall,
                    .adjPvalue,
                    nam.data = paste0(cell.line, "_0_48"),
                    nam.situ = "up.0"
                )
            )
        )
    )
    rownames(tmp.heatmatrix) <- tmp.pathwaylist[tmp.pathwayall, "Name"]
    ## p <- heatMapPlot(
    ##     as.matrix(tmp.heatmatrix),
    ##     ycluster  = TRUE,
    ##     color     = c("blue", "white", "red"),
    ##     ylab.size = 11
    ## )
    ## ggsave(
    ##     filename = file.path(
    ##         rna.dir[["cummeRbund"]],
    ##         paste0(cell.line, "_kegg.pdf")
    ##     ),
    ##     plot     = p
    ## )
    tmp.heatmatrix2 <- data.frame(
        H0_H6 = ifelse(
            abs(tmp.heatmatrix$H0_H6.down) > abs(tmp.heatmatrix$H0_H6.up),
            tmp.heatmatrix$H0_H6.down,
            tmp.heatmatrix$H0_H6.up
        ),
        H6_H24 = ifelse(
            abs(tmp.heatmatrix$H6_H24.down) > abs(tmp.heatmatrix$H6_H24.up),
            tmp.heatmatrix$H6_H24.down,
            tmp.heatmatrix$H6_H24.up
        ),        
        H0_H24 = ifelse(
            abs(tmp.heatmatrix$H0_H24.down) > abs(tmp.heatmatrix$H0_H24.up),
            tmp.heatmatrix$H0_H24.down,
            tmp.heatmatrix$H0_H24.up
        ),
        H24_H48 = ifelse(
            abs(tmp.heatmatrix$H24_H48.down) > abs(tmp.heatmatrix$H24_H48.up),
            tmp.heatmatrix$H24_H48.down,
            tmp.heatmatrix$H24_H48.up
        ),        
        H0_H48 = ifelse(
            abs(tmp.heatmatrix$H0_H48.down) > abs(tmp.heatmatrix$H0_H48.up),
            tmp.heatmatrix$H0_H48.down,
            tmp.heatmatrix$H0_H48.up
        )
    )
    rownames(tmp.heatmatrix2) <- tmp.pathwaylist[tmp.pathwayall, "Name"]
    p2 <- heatMapPlot(
        as.matrix(tmp.heatmatrix2),
        ycluster  = FALSE,
        color     = c("blue", "white", "red"),
        xlab.size = 10,
        ylab.size = 10,
        yorder    = tmp.pathwayorder2[[cell.line]][length(tmp.pathwayorder2[[cell.line]]):1]
    )
    ggsave(
        filename = file.path(
            rna.dir[["cummeRbund"]],
            paste0(cell.line, "_kegg_merge.pdf")
        ),
        plot     = p2,
        width    = 30,
        height   = 15,
        units    = "cm"
    )
}

### Plot heatmap

tmp


###}}}

###{{{ GO: goseq and kobas

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
        } else {}
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

rna.kobasgo <- list()

for (nam.data in names(rna.go)) {
    rna.kobasgo[[nam.data]] <- list()
    for (nam.situ in c("all", "down.0", "down.1", "up.0", "up.1")) {
        rna.kobasgo[[nam.data]][[nam.situ]] <- read.table(
            file = file.path(
                rna.dir[["cummeRbund"]],
                nam.data,
                paste(nam.data, nam.situ,"kobasgo.csv", sep = "_")
            ),
            header = TRUE,
            sep    = "\t",
            stringsAsFactors = FALSE
        )
        tmp.go <- rna.kobasgo[[nam.data]][[nam.situ]][
            which(
                rna.kobasgo[[nam.data]][[nam.situ]]$Corrected.P.Value < 0.05
            ),]
        if (dim(tmp.go)[1] > 20) {
            tmp.go <- tmp.go[order(tmp.go$Corrected.P.Value)[1:20],]
        } else {}
        tmp.order <- order(tmp.go$Corrected.P.Value)
        tmp.data  <- data.frame(
            go          = factor(
                tmp.go$Term,
                levels = as.character(tmp.go$Term)[tmp.order],
            ),          
            log10.adj.p = -log10(tmp.go$Corrected.P.Value),
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
                paste(nam.data, nam.situ, "kobasgo.pdf", sep = "_")
            ),
            plot = p
        )
        write.table(
            tmp.go,
            file = file.path(
                rna.dir[["cummeRbund"]],
                nam.data,
                paste(nam.data, nam.situ, "kobasgo_sig.csv", sep = "_")
            ),
            sep = ",",
            row.names = FALSE
        )
        rm(tmp.order, tmp.data, tmp.go, p)
    }              
}


###}}}

###{{{ TFT


rna.tft <- list()

for (nam.de in names(rna.de.geneid)) {
    rna.tft[[nam.de]] <- list()
    for (nam.situ in names(rna.de.geneid[[nam.de]])) {
        rna.tft[[nam.de]][[nam.situ]] <- tftEnrichment(
            gene.symbol.list = na.omit(
                rna.de.geneid[[nam.de]][[nam.situ]]$hgnc_symbol
            ),
            tft.list      = lapply(
                subset(
                    rna.info[["gsea.tft"]],
                    !grepl("UNKNOWN", names(rna.info[["gsea.tft"]]))
                ),
                function(x) {
                    return(x$gene)
                }
            ),
            whole.gene.number = dim(rna.test.result[[nam.de]])[1],
            threshold         = 0.05,
            adjust.p.method   = "BH"
        )
        tmp.order <- order(rna.tft[[nam.de]][[nam.situ]]$adjusted.P.value)
        tmp.data  <- data.frame(
            tft     = factor(
                as.character(rna.tft[[nam.de]][[nam.situ]]$TFT),
                levels = as.character(
                    rna.tft[[nam.de]][[nam.situ]]$TFT)[tmp.order],
            ),          
            log10.adj.p = -log10(rna.tft[[nam.de]][[nam.situ]]$adjusted.P.value),
            order       = tmp.order
        )
        if (dim(tmp.data)[1] > 20) {
            tmp.data <- tmp.data[1:20,]
        }
        p         <- termPlot(
            data    = tmp.data,
            term    = tft,
            value   = log10.adj.p,
            fill    = order,
            title   = paste(
                nam.de,
                nam.situ,
                "TFT",
                sep = " "
            ),
            xlab    = "Transcript Factor Targets",
            ylab    = "",
            palette = ifelse(grepl("up", nam.situ), "Reds", "Blues")
        )
        ggsave(
            filename = file.path(
                rna.dir[["cummeRbund"]],
                nam.de,
                paste(nam.de, nam.situ, "tft.pdf", sep = "_")
            ),
            plot = p
        )
        rm(tmp.order, tmp.data)
    }
}

###}}}

###{{{ For circos


rna.de.aa <- lapply(
    rna.de.geneid,
    function(x) {
        lapply(
            x,
            function(y) {
               y[y$hgnc_symbol %in% rna.info[["aa.gene"]]$Gene,]
            }
        )
    }
)

compare.time <- with(
    expand.grid(c(0, 6, 24, 48), c(0, 6, 24, 48)),
    {
        paste(Var1[Var1 < Var2], Var2[Var1 < Var2], sep = "_")
    }
)


situ <- names(rna.sig.aa[[1]])

rna.circos <- list()

for (i in seq_along(compare.time)) {
    rna.circos[[compare.time[i]]] <- lapply(
        situ[c(1,2,4)],
        function(x) {
            tmp.jeko      <- paste("jeko", compare.time[i], sep = "_")
            tmp.raji      <- paste("raji", compare.time[i], sep = "_")
            tmp.gene.list <- unique(
                c(
                    rna.de.aa[[tmp.jeko]][[x]]$hgnc_symbol,
                    rna.de.aa[[tmp.raji]][[x]]$hgnc_symbol
                )
            )
            tmp <- data.frame(
                gene = tmp.gene.list,
                jeko = unlist(
                    lapply(
                        tmp.gene.list,
                        function(x) {
                            with(
                                rna.test.result[[tmp.jeko]],
                                ifelse(
                                    sum(gene_id == x) == 1 &
                                    q_value[gene_id == x] < 0.05,
                                    -log10(q_value[gene_id == x][1]),
                                    0
                                )
                            )
                        }
                    )
                ),
                raji = unlist(
                    lapply(
                        tmp.gene.list,
                        function(x) {
                            with(
                                rna.test.result[[tmp.raji]],
                                ifelse(
                                    sum(gene_id == x) == 1 &
                                    q_value[gene_id == x] < 0.05,
                                    -log10(q_value[gene_id == x][1]),
                                    0
                                )
                            )
                        }
                    )
                ),
                stringsAsFactors = FALSE
            )
            rownames(tmp) <- tmp.gene.list
            tmp
        }        
    )
    names(rna.circos[[compare.time[i]]]) <- situ[c(1,2,4)]
}


for (nam.comp in names(rna.circos)) {
    for (nam.situ in names(rna.circos[[nam.comp]])) {
        tmp.path <- file.path(
            rna.dir[["circos"]],
            paste(nam.comp, nam.situ, sep = "_")
        )
        if (!dir.exists(tmp.path )) {
            dir.create(tmp.path)
            dir.create(file.path(tmp.path, "etc"))
            dir.create(file.path(tmp.path, "data"))
            dir.create(file.path(tmp.path, "results"))
            file.copy(
                file.path(
                    rna.dir[["circos"]],
                    "etc",
                    dir(file.path(rna.dir[["circos"]], "etc"))
                ),
                file.path(tmp.path, "etc")
            )
            file.copy(
                file.path(
                    rna.dir[["circos"]],
                    "data",
                    dir(file.path(rna.dir[["circos"]], "data"))
                ),
                file.path(tmp.path, "data")
            )
        } else {}
        data       <- as.matrix(rna.circos[[nam.comp]][[nam.situ]][,-1])
        colnames(data) <- c("JeKo-1", "Raji")
        tmp.length <- table2KaryotypeLength(data)
        tmp.color  <- makeKaryotypeColor(names(tmp.length), alpha = 0.3)
        tmp.color[colnames(data)] <- gray(
            seq_along(colnames(data)) / (length(colnames(data)) + 1)
        )
        outputColorConf(
            color = tmp.color,
            file  = file.path(tmp.path, "data/colors.conf")
        )
        outputKaryotype(
            data = data,
            file = file.path(tmp.path, "data/karyotype.txt")
        )
        outputLink(
            data = data,
            file = file.path(tmp.path, "data/links.txt")
        )
        write.table(
            rna.circos[[nam.comp]][[nam.situ]],
            file = file.path(
                tmp.path,
                paste0(nam.comp, "_", nam.situ, ".csv")
            ),
            row.names = FALSE,
            quote     = FALSE,
            sep       = "\t"
        )
    }
}

###}}}

###{{{ Genes

rna.count.mat <- list()

rna.count.mat[["raw"]] <- as.matrix(as.data.frame(
    lapply(
        rna.count,
        function(x) {
            x[rownames(rna.count[[1]]), "count"]
        }
    )
))

rownames(rna.count.mat[["raw"]]) <- rownames(rna.count[[1]])
colnames(rna.count.mat[["raw"]]) <- gsub(
    "[.]", "_",
    colnames(rna.count.mat[["raw"]])
)

rna.count.mat[["norm"]] <- normalize.quantiles(
    rna.count.mat[["raw"]],
    copy = TRUE
)

rownames(rna.count.mat[["norm"]]) <- rownames(rna.count.mat[["raw"]])
colnames(rna.count.mat[["norm"]]) <- colnames(rna.count.mat[["raw"]])

rna.fpkm <- list()
rna.fpkm[["jeko"]] <- as.matrix(data.frame(
    H0  = rna.test.result[[1]]$value_1,
    H6  = rna.test.result[[3]]$value_2 * 1.06,
    H24 = rna.test.result[[1]]$value_2,
    H48 = rna.test.result[[2]]$value_2 * 0.96
))
rownames(rna.fpkm[["jeko"]]) <- rna.test.result[[1]]$gene_id
rna.fpkm[["raji"]] <- as.matrix(data.frame(
    H0  = rna.test.result[[7]]$value_1,
    H6  = rna.test.result[[9]]$value_2 * 1.01,
    H24 = rna.test.result[[7]]$value_2,
    H48 = rna.test.result[[8]]$value_2 * 1.06
))
rownames(rna.fpkm[["raji"]]) <- rna.test.result[[7]]$gene_id

mygene <- c(
    "MTOR", "TP53",
    grep("[a-zA-z]", unique(unlist(rna.de.aa)), value = TRUE)
)

mygene <- "APAF1"

for (nam.gene in mygene) {
    for (nam.cell in c("jeko", "raji")) {
        tmp.expression <- data.frame(
            Time = c(0, 6, 24, 48),
            FPKM = rna.fpkm[[nam.cell]][nam.gene,]
        )
        p <- ggplot(tmp.expression, aes(x = Time, y = FPKM)) +
            geom_line(
                size     = 2,
                linetype = 2,
                lineend  = "round",
                colour   = I("lightskyblue")
            ) + geom_point(
                    size   = 6,
                    colour = I("coral")
                ) + theme_bw() + xlab("Time") + ylab("FPKM") +
            ggtitle(paste(nam.cell, nam.gene, sep = " "))
        ggsave(
            filename = file.path(
                rna.dir[["genes"]],
                paste(nam.gene, "_", nam.cell, ".pdf", sep = "")
            ),
            plot = p
        )
    }
    tmp.exp <- data.frame(
        Time = rep(c(0, 6, 24, 48), 2),
        Cell = rep(c("JeKo-1", "Raji"), each = 4),
        FPKM = c(rna.fpkm[["jeko"]][nam.gene,], rna.fpkm[["raji"]][nam.gene,])
    )
    p <- ggplot(
        tmp.exp,
        aes(x = Time, y = FPKM, group = Cell, colour = Cell)
    ) +
        geom_line(
            size     = 2,
            linetype = 2,
            lineend  = "round"
        ) + geom_point(
                size   = 6
            ) + theme_bw() + xlab("Time") + ylab("FPKM") + ggtitle(nam.gene)
    ggsave(
        filename = file.path(
            rna.dir[["genes"]],
            paste0(nam.gene, ".pdf")
        ),
        plot     = p
    )
}

###}}}

###{{{ Venn Plot

rna.venn <- list()

rna.venn[["jeko"]] <- venn(
    list(
        H0_H6   = as.character(rna.de.geneid[["jeko_0_6"]][["all"]]$hgnc_symbol),
        H0_H24  = as.character(rna.de.geneid[["jeko_0_24"]][["all"]]$hgnc_symbol),
        H0_H48  = as.character(rna.de.geneid[["jeko_0_48"]][["all"]]$hgnc_symbol)
    ),
    show.plot = FALSE
)

venn.diagram(
    list(
        H0_H6   = as.character(rna.de.geneid[["jeko_0_6"]][["all"]]$hgnc_symbol),
        H0_H24  = as.character(rna.de.geneid[["jeko_0_24"]][["all"]]$hgnc_symbol),
        H0_H48  = as.character(rna.de.geneid[["jeko_0_48"]][["all"]]$hgnc_symbol)
    ),
    filename  = file.path(rna.dir[["basic"]], "jeko.svg"),
    imagetype = "svg",
    width     = 6,
    height    = 6,
    units     = "cm",
    fill      = brewer.pal(3, "Set1"),
    cat.col   = brewer.pal(3, "Set1")
)

rna.venn[["raji"]] <- venn(
    list(
        H0_H6   = as.character(rna.de.geneid[["raji_0_6"]][["all"]]$hgnc_symbol),
        H0_H24  = as.character(rna.de.geneid[["raji_0_24"]][["all"]]$hgnc_symbol),
        H0_H48  = as.character(rna.de.geneid[["raji_0_48"]][["all"]]$hgnc_symbol)
    ),
    show.plot = FALSE
)

venn.diagram(
    list(
        H0_H6   = as.character(rna.de.geneid[["raji_0_6"]][["all"]]$hgnc_symbol),
        H0_H24  = as.character(rna.de.geneid[["raji_0_24"]][["all"]]$hgnc_symbol),
        H0_H48  = as.character(rna.de.geneid[["raji_0_48"]][["all"]]$hgnc_symbol)
    ),
    filename  = file.path(rna.dir[["basic"]], "raji.svg"),
    imagetype = "svg",
    width     = 6,
    height    = 6,
    units     = "cm",
    fill      = brewer.pal(3, "Set1"),
    cat.col   = brewer.pal(3, "Set1")
)

write.table(
    attr(rna.venn[["jeko"]], "intersections")[["H0_H6:H0_H24:H0_H48"]],
    file      = "jeko.csv",
    quote     = FALSE,
    row.names = FALSE,
    col.names = FALSE
)

write.table(
    attr(rna.venn[["raji"]], "intersections")[["H0_H6:H0_H24:H0_H48"]],
    file      = "raji.csv",
    quote     = FALSE,
    row.names = FALSE,
    col.names = FALSE
)

write.table(
    intersect(
        attr(rna.venn[["raji"]], "intersections")[["H0_H6:H0_H24:H0_H48"]],
        attr(rna.venn[["jeko"]], "intersections")[["H0_H6:H0_H24:H0_H48"]]
    ),
    file      = "both.csv",
    quote     = FALSE,
    row.names = FALSE,
    col.names = FALSE
)

rna.venn.go <- list()

for (nam in c("jeko", "raji", "both")) {
    rna.venn.go[[nam]] <- read.table(
        file   = paste0(nam, "_go.txt"),
        sep    = "\t",
        header = TRUE
    )
    tmp.go      <- rna.venn.go[[nam]][rna.venn.go[[nam]]$Benjamini < 0.05,]
    if (dim(tmp.go)[1] > 20) {
        tmp.go  <- tmp.go[order(tmp.go$Benjamini)[1:20],]
    } else {}
    tmp.go$Term <- unlist(lapply(strsplit(as.character(tmp.go$Term), "~"), "[[", 2))
    tmp.order   <- order(tmp.go$Benjamini)
    tmp.data    <- data.frame(
        GO          = factor(
            tmp.go$Term,
            levels = as.character(tmp.go$Term)[tmp.order]
        ),
        log10.benjamini = -log10(tmp.go$Benjamini),
        order       = tmp.order
    )
    p         <- termPlot(
        data    = tmp.data,
        term    = GO,
        value   = log10.benjamini,
        fill    = order,
        title   = paste(
            ifelse(nam == "jeko", "JeKo-1", "Raji"),
            "GO BP",
            sep = " "
        ),
        xlab    = "GO Term",
        ylab    = "",
        palette = "Blues"
    )
    ggsave(
        filename = paste(nam, "go.pdf", sep = "_"),
        plot = p
    )
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


rna.de.aaratio <- data.frame(
    comparison = rep(
        x    = names(rna.de.geneid),
        each = length(names(rna.de.geneid[[1]]))
    ),
    situation  = rep(
        x     = names(rna.de.geneid[[1]]),
        times = length(names(rna.de.geneid))
    ),
    ratio      = unlist(
        lapply(
            names(rna.de.geneid),
            function(x) {
                lapply(
                    names(rna.de.geneid[[x]]),
                    function(y) {
                        dim(rna.de.aa[[x]][[y]])[1] /
                            dim(rna.de.geneid[[x]][[y]])[1]
                    }
                )
            }
        )
    )
)

tmp.aaratio <- rna.de.aaratio[
(rna.de.aaratio$comparison %in% names(rna.de.geneid)[-c(6, 12)]) &
!grepl("[.]1", rna.de.aaratio$situation),]

tmp.aaratio$comparison <- factor(
    as.character(tmp.aaratio$comparison),
    levels = c("jeko_0_6", "jeko_6_24", "jeko_0_24", "jeko_24_48", "jeko_0_48",
               "raji_0_6", "raji_6_24", "raji_0_24", "raji_24_48", "raji_0_48")
)

tmp.aaratio$situation <- factor(
    unlist(lapply(
        strsplit(as.character(tmp.aaratio$situation), split = "[.]"),
        "[",
        1
    )),
    levels = c("all", "up", "down")
)

p <- ggplot(
    tmp.aaratio,
    aes(x = comparison, y = ratio, fill = situation)
) + geom_bar(
        stat     = "identity",
        position = "dodge"
    ) + theme_classic() +
    theme(
        axis.text.x      = element_text(
            angle  = 45,
            vjust  = 0.8,
            hjust  = 0.9
        )
    )

ggsave(
    filename = file.path(
        rna.dir[["info"]],
        paste0("aaratio.pdf")
    ),
    plot     = p,
    width    = 20,
    height   = 10,
    units    = "cm"
)


###}}}

################################################################################
##                               SCRATCH                                      ##
################################################################################

###{{{ Scratch


###}}}


        
