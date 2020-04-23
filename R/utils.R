##' construct RGChannelSet in parallel using BiocParallel
##'
##'
##' @title construct RGChannelSet in parallel
##' @param targets data.frame representing valid targets file as is
##' used within the minfi package
##' @param verbose logical default TRUE
##' @param ... optional arguments to read.metharray.exp
##' @return RGset
##' @author mvaniterson
##' @export
##' @import minfi
##' @importFrom BiocParallel bplapply bpworkers bpparam
##' @importFrom Biobase combine
##' @importFrom utils str
read.metharray.exp.par <- function(targets, verbose = TRUE, ...) {
    nworkers <- bpworkers(bpparam())
    if (nworkers <= 1)
        stop("Did you registered a biocparallel back-end?")
    y <- rep(1, ceiling(nrow(targets)/nworkers))
    for (i in 2:nworkers) y <- c(y, rep(i, ceiling(nrow(targets)/nworkers)))
    y <- y[1:nrow(targets)]
    jobs <- split(targets, y)

    fun <- function(x, ...) {
        ##these need to be loaded on the worker nodes explicitly for BatchJobs!
        requireNamespace("minfi")
        requireNamespace("Biobase")
        read.metharray.exp(targets = x, ...)
    }

    message("Reading multiple idat-files in parallel")
    res <- bplapply(jobs, FUN = fun, ...)
    if(verbose)
        message(str(res))
    message("Combining the RGsets to one big RGset")
    rgSet <- res[[1]]
    for (i in 2:length(res)) rgSet <- combine(rgSet, res[[i]])
    rgSet
}

##' detectionP that accepts NA's
##'
##' since probe-filtering can introduce NA's
##' detectionP is modified to handle NA's properly
##' beware the P-value matrix thus can contain NA's
##' @title detectionP that accepts NA's
##' @param rgSet RGset
##' @param type "m+u"
##' @param na.rm FALSE/TRUE
##' @return matrix with detection P-values
##' @author mvaniterson
##' @importFrom matrixStats colMedians colMads
##' @importFrom stats pnorm
detectionP <- function(rgSet, type = "m+u", na.rm=FALSE) {
    locusNames <- getManifestInfo(rgSet, "locusNames")
    detP <- matrix(NA_real_, ncol = ncol(rgSet), nrow = length(locusNames),
                   dimnames = list(locusNames, sampleNames(rgSet)))

    controlIdx <- getControlAddress(rgSet, controlType = "NEGATIVE")
    r <- getRed(rgSet)
    rBg <- r[controlIdx,]
    rMu <- matrixStats::colMedians(rBg, na.rm=na.rm)
    rSd <- matrixStats::colMads(rBg, na.rm=na.rm)

    g <- getGreen(rgSet)
    gBg <- g[controlIdx,]
    gMu <- matrixStats::colMedians(gBg, na.rm=na.rm)
    gSd <- matrixStats::colMads(gBg, na.rm=na.rm)

    TypeII <- getProbeInfo(rgSet, type = "II")
    TypeI.Red <- getProbeInfo(rgSet, type = "I-Red")
    TypeI.Green <- getProbeInfo(rgSet, type = "I-Green")
    for (i in 1:ncol(rgSet)) {
        ## Type I Red
        intensity <- r[TypeI.Red$AddressA, i] + r[TypeI.Red$AddressB, i]
        detP[TypeI.Red$Name, i] <- 1-pnorm(intensity, mean=rMu[i]*2, sd=rSd[i]*2)
        ## Type I Green
        intensity <- g[TypeI.Green$AddressA, i] + g[TypeI.Green$AddressB, i]
        detP[TypeI.Green$Name, i] <- 1-pnorm(intensity, mean=gMu[i]*2, sd=gSd[i]*2)
        ## Type II
        intensity <- r[TypeII$AddressA, i] + g[TypeII$AddressA, i]
        detP[TypeII$Name, i] <- 1-pnorm(intensity, mean=rMu[i]+gMu[i], sd=rSd[i]+gSd[i])
    }
    detP
}

##' Determine sex based on X chr methylation status
##'
##' @title Determine sex based on X chr methylation status
##' @param beta beta matrix
##' @param cutbeta beta values threshold use values between [0.2-0.6]
##' @param nx number of X chr probes to determine sex is Male default 3000
##' @return sex prediction
##' @author rcslieker
##' @import FDb.InfiniumMethylation.hg19
##' @importFrom GenomicFeatures features
##' @export 
getSex.DNAmArray <- function(beta, cutbeta=c(0.2, 0.6), nx = 3000){
    InfMet <- features(FDb.InfiniumMethylation.hg19)    
    chrX <- names(InfMet[as.vector(seqnames(InfMet)) %in% "chrX"])    
    chrX <- chrX[grep("cg", chrX)]
    BetaX <- beta[match(chrX, rownames(beta)),]
    nopoX <- colSums(BetaX >= cutbeta[1] & BetaX <= cutbeta[2], na.rm=TRUE)
    ifelse(nopoX <= 3000, "Male", "Female")
}

##' Printer friendly qq-plot
##'
##' Use interpolation on the 'uninteresting' P-values to reduce the number
##' of points that are drawn.
##' @title Printer friendly qq-plot
##' @param pval vector of p-values
##' @param p proportion that should be interpolated
##' @param k number of points used to interpolate
##' @param add add qq-plot to existing qq-plot
##' @param show.fit show the interpolation fit
##' @param pch default plot pars
##' @param col default plot pars
##' @param nhighlight number of highlighted points
##' @param bty default plot pars
##' @param main default plot pars
##' @param ... additional graphical parameters
##' @return plot
##' @author mvaniterson
##' @export
##' @importFrom stats quantile rexp splinefun pnorm prcomp
##' @importFrom graphics abline barplot curve plot points text
qqpf <- function(pval, p=0.90, k=7, add=FALSE, show.fit=FALSE, pch=16, col=1, nhighlight=6, bty="n", main="", ...) {

    lbs <- names(pval)
    pval <- pval[!is.na(pval) & !is.nan(pval) & !is.null(pval) & is.finite(pval) & pval < 1 & pval >= 0]

    pval[pval == 0] <- min(pval[pval>0])/10

    y <- sort(-log10(pval))
    lbs <- lbs[order(-log10(pval))]

    set.seed(12345) ##
    x <- sort(rexp(length(pval), log(10)))

    inputArgs <- list(...)

    if(any(names(inputArgs) == "ylim"))
        ylim <- inputArgs[["ylim"]]
    else
        ylim <- range(y)

    if(any(names(inputArgs) == "xlim"))
        xlim <- inputArgs[["xlim"]]
    else
        xlim <- range(x)

    id <- y < quantile(y, prob=p)
    xid <- x[id][seq(1, sum(id), (sum(id) - 1)/k)]
    yid <- y[id][seq(1, sum(id), (sum(id) - 1)/k)]
    f <- splinefun(xid, yid, method="natural")

    highlight <- function(x, y, lbs, nhighlight, col, cex=0.5){
        id <- order(y, decreasing=TRUE)[1:nhighlight]
        text(x[id]-1.2, y[id], lbs[id], col=col, cex=cex)
        ##arrows(x[id]-.2, y[id], x[id], y[id], length=0, col=col, cex=0.5)
    }

    if(add) {
        points(x[!id], y[!id], pch=pch, col=col)
        if(nhighlight>0)
            highlight(x[!id], y[!id], lbs[!id], nhighlight, col=col)
        curve(expr=f, add=TRUE, xlim=range(x[id]), lwd=7, col=col)
    }
    else {
        if(!show.fit) {
            plot(x[!id], y[!id], xlim=xlim, ylim=ylim,
                 xlab = expression(Expected ~ ~-log[10](italic(p))),
                 ylab = expression(Observed ~ ~-log[10](italic(p))),
                 pch=pch, col=col, bty=bty, main=main)
            if(nhighlight>0)
                highlight(x[!id], y[!id], lbs[!id], nhighlight, col=col)
            curve(expr=f, add=TRUE, xlim=range(x[id]), lwd=7, col=col)
            abline(0, 1, col=1, lty=1)
        }
        else {
            plot(x, y, xlim=range(x), ylim=range(y),
                 xlab="Theoretical (exp(ln10))", ylab="Observed (-log10(P-values))",
                 pch=pch, col=col, bty=bty, main=main)
            if(nhighlight>0)
                highlight(x, y, lbs, nhighlight)
            curve(expr=f, add=TRUE, xlim=range(x[id]), lwd=1, col=2)
            abline(0, 1, col=1, lty=1)
        }
    }
}

##' find gene nearest to list fo cpgs
##'
##' find gene nearest to list fo cpgs
##' @title nearestGenes
##' @param cpgs list of illumina cpg identifiers
##' @param TxDb name of Transcript database, i.e. TxDb.Hsapiens.UCSC.hg19.knownGene
##' @return data.frame with cpg annotations
##' @author mvaniterson
##' @import FDb.InfiniumMethylation.hg19 org.Hs.eg.db GenomicRanges
##' @importFrom AnnotationDbi select
##' @export
cpgInfo <- function(cpgs, TxDb) {
    
    if(!requireNamespace(TxDb, character.only = TRUE))
        stop("TxDb:", TxDb, " not available!")

    ##extract annotation
    genes <- genes(eval(parse(text=TxDb)))
    rowRanges <- getPlatform(platform = "HM450", genome = "hg19")
    ##find nearest genes
    gr <- rowRanges[names(rowRanges) %in% cpgs]
    ##genes[precede(gr, genes)]    
    genes <- genes[nearest(gr, genes)]

    ##add gene symbol
    if(grepl("ensGene", TxDb)){
        map <- select(org.Hs.eg.db, names(genes), columns="SYMBOL", keytype="ENSEMBL")
        map$SYMBOL[is.na(map$SYMBOL)] <- map$ENSEMBL[is.na(map$SYMBOL)]
        id <- match(names(genes), map$ENSEMBL)
    } else {
        map <- select(org.Hs.eg.db, names(genes), columns="SYMBOL", keytype="ENTREZID")
        map$SYMBOL[is.na(map$SYMBOL)] <- map$ENTREZID[is.na(map$SYMBOL)]
        id <- match(names(genes), map$ENTREZID)
    }
    mcols(genes)$SYMBOL <- map$SYMBOL[id]
    
    ##reorder
    genes <- as.data.frame(genes, row.names=names(gr))
    gr <- as.data.frame(mcols(gr), row.names=names(gr))
    genes <- merge(genes, gr, by="row.names")
    rownames(genes) <- genes$Row.names
    genes <- genes[,-1]    
    genes[match(cpgs, rownames(genes)),]   
}

