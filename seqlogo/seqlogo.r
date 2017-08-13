# unfinished
seqs <- t(as.data.frame(
    lapply(
        seq(50),
        function(x) {
            return(sample(c('A', 'T', 'C', 'G'), 6, replace=TRUE))
        }
    )
))


rownames(seqs) <- NULL

seqs2 <- t(matrix(
    c('A', 'T', 'C', 'C', 'A', 'C', 'A', 'G', 'A', 'T', 'G', 'G'),
    4, 3
))

allchar <- unique(as.vector(seqs2))

ntable <- apply(
    seqs2, MARGIN=2,
    function(x) {
        re <- table(x)
        for (ch in allchar) {
            if (!(ch %in% names(re))) {
                tmp <- c(re, 0)
                names(tmp) <- c(names(re), ch)
                re <- tmp
            } else {}
        }
        re <- re[allchar]
        return(re)
    }
)

f <- apply(
    ntable,
    MARGIN=2,
    function(x) {
        return(x / sum(x))
    }
)

h <- apply(
    f,
    MARGIN=2,
    function(x) {
        return(-sum(x * log2(x), na.rm=TRUE))
    }
)

r = log2(length(allchar)) - h

height <- f * r

xlim <- c(0, dim(seqs2)[2] + 0.5)

ylim <- c(0, max(colSums(height)) + 0.5)

default.par <- par()
par(bg='transparent')
plot(
    0, 0,
    col=NULL,
    xlim=c(0, 4), ylim=c(0, 4),
    axes=FALSE,
    xlab='', ylab=''
)
par(new=FALSE)
for (coli in seq(dim(height)[2])) {
    coldata <- sort(height[,coli])
}
