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
    for (i in 2:nworkers) rgSet <- combine(rgSet, res[[i]])
    rgSet
}

##' Determine sex based on X chr methylation status
##'
##' @title Determine sex based on X chr methylation status
##' @param beta beta matrix
##' @param cutbeta beta values threshold use values between [0.2-0.6]
##' @param nx number of X chr probes to determine sex is Male default 3000
##' @param array EPIC or 450K
##' @param genome hg19 or hg38
##' @return sex prediction
##' @author ljsinke
##' @importFrom readr read_tsv
##' @export 
getSex.DNAmArray <- function(beta, cutbeta=c(0.2, 0.6), nx = 3000, array = '450K', genome = 'hg19'){
  if(array=="EPIC" & genome=="hg19") {
    maskProbes <- read_tsv(paste0(path.package("DNAmArray"), "/extdata/EPIC.hg19.manifest.txt.gz"))
  }
  if(array=="EPIC" & genome=="hg38") {
    maskProbes <- read_tsv(paste0(path.package("DNAmArray"), "/extdata/EPIC.hg38.manifest.txt.gz"))
  }
  if(array=="450K" & genome=="hg19") {
    maskProbes <- read_tsv(paste0(path.package("DNAmArray"), "/extdata/HM450.hg19.manifest.txt.gz"))
  }
  if(array=="450K" & genome=="hg38") {
    maskProbes <- read_tsv(paste0(path.package("DNAmArray"), "/extdata/HM450.hg38.manifest.txt.gz"))
  }
  
  chrX <- names(maskProbes[seqnames(maskProbes) %in% 'chrX'])
  chrX <- chrX[grep("cg", chrX)]
  
  betaX <- beta[rownames(beta) %in% chrX,]
  betaX <- assay(betaX)
  
  nopoX <- colSums(betaX >= cutbeta[1] & betaX <= cutbeta[2], na.rm=TRUE)
  ifelse(nopoX <= nx, "Male", "Female")
}

##' find gene nearest to list of cpgs
##'
##' find gene nearest to list of cpgs
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
##' @importFrom stats pnorm median mad
detectionP <- function(rgSet, type = "m+u", na.rm=FALSE) {
  locusNames <- getManifestInfo(rgSet, "locusNames")
  detP <- matrix(NA_real_, ncol = ncol(rgSet), nrow = length(locusNames),
                 dimnames = list(locusNames, sampleNames(rgSet)))
  
  controlIdx <- getControlAddress(rgSet, controlType = "NEGATIVE")
  r <- getRed(rgSet)
  rBg <- r[controlIdx,]
  rMu <- apply(rBg, 2, median, na.rm=na.rm)
  rSd <- apply(rBg, 2, mad, na.rm=na.rm)
  
  g <- getGreen(rgSet)
  gBg <- g[controlIdx,]
  gMu <- apply(gBg, 2, median, na.rm=na.rm)
  gSd <- apply(gBg, 2, mad, na.rm=na.rm)
  
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
