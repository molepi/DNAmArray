##' filter on probe-level: number of beads and zero intensities
##'
##'
##' @title filter on probe-level: number of beads and zero intensities
##' @param RGset a RGChannelSetExtended
##' @param cutbead threshold number of beads
##' @param zeroint filter out zero intenseties default is TRUE
##' @param verbose Default is TRUE
##' @return RGChannelSet
##' @author mvaniterson
##' @export
##' @import minfi
##' @importFrom SummarizedExperiment colData
##' @importFrom Biobase varMetadata AnnotatedDataFrame phenoData featureData experimentData annotation protocolData assayDataElement
probeFiltering <- function(RGset, cutbead=3, zeroint=TRUE, verbose=TRUE){

    ##Filter on number of beads
    if(verbose)
        cat("Filtering on number of beads... \n")

    beadmat <- getNBeads(RGset)
    
    idBeadmat <- beadmat < cutbead
    ##beadmat[idBeadmat] <- NA

    if(verbose)
        cat("On average", round(100*sum(idBeadmat)/prod(dim(idBeadmat)), 2),"% of the probes (",nrow(idBeadmat),") have number of beads below", cutbead, "\n")

    ##Filter on Red and Green intensity <1
    if(zeroint) {
        if(verbose)
            cat("Filtering on zero intensities... \n")

        Grn <- getGreen(RGset)
        Red <- getRed(RGset)

        ##determine if Grn and/or Red intensities of type II probes are <1
        idT2 <- Grn[getProbeInfo(RGset, type = "II")$AddressA,] < 1 | Red[getProbeInfo(RGset, type = "II")$AddressA,] < 1

        ##determine if either Grn or Red intensities of Type I probes are <1
        idT1Grn <- Grn[c(getProbeInfo(RGset, type = "I-Green")$AddressA,
                         getProbeInfo(RGset, type = "I-Green")$AddressB),] < 1

        idT1Red <- Red[c(getProbeInfo(RGset, type = "I-Red")$AddressA,
                         getProbeInfo(RGset, type = "I-Red")$AddressB),] < 1

        if(verbose) {
            cat("On average", round(100*sum(idT2)/prod(dim(idT2)), 3),"% of the Type II probes (",nrow(idT2),") have Red and/or Green intensity below 1 \n")
            cat("On average", round(100*sum(idT1Grn)/prod(dim(idT1Grn)), 3),"% of the Type I probes (",nrow(idT1Grn),"), measured in Green channel, have intensity below 1 \n")
            cat("On average", round(100*sum(idT1Red)/prod(dim(idT1Red)), 3),"% of the Type I probes (",nrow(idT1Red),"), measured in Red channel, have intensity below 1 \n")
        }
    }

    ##combine all filtered results and set NA in Red and/or Green channels
    Red[idBeadmat] <- Grn[idBeadmat] <- NA

    if(zeroint) {
        if(verbose){
            cat("Set filtered probes in Red and/or Green channels to NA... \n")
        }

        for(i in 1:ncol(RGset)) {
            if(verbose & i%%100 == 0)
                cat("... done ",i," out of ",ncol(RGset)," ... \n")
            idRed <- c(names(which(idT2[,i])), names(which(idT1Red[,i])))
            midRed <- match(idRed, rownames(Red))
            Red[midRed, i] <- NA
            idGrn <- c(names(which(idT2[,i])), names(which(idT1Grn[,i])))
            midGrn <- match(idGrn, rownames(Grn))
            Grn[midGrn, i] <- NA
        }
    }

    RGChannelSet(Green = Grn, Red = Red,
                 colData = colData(RGset),
                 annotation = annotation(RGset))
}


##' Extract functional normalized data according to filter data
##'
##' Since functional normalization doesn't accept NA e.g. after probe filter the
##' normalized 'GenomicRatioSet' and filtered 'RGChannelSet' need to be merged.
##' Additionally if M-values are calculated we set -/+Inf values to +/-16
##' more details
##' @title Extract functional normalized data according to filter data
##' @param GRset 'GenomicRatioSet' output from functional normalization
##' @param RGset 'RGset' output from filtering
##' @param what return type either M- or beta-values
##' @param cutp detection p-value threshold
##' @param cutsamples threshold removing samples
##' @param cutcpgs threshold removing probes
##' @param verbose Default is TRUE
##' @param ... optional arguments getBeta
##' @return matrix with  M- or beta-values
##' @author mvaniterson
##' @export
##' @importFrom utils data
reduce <- function(GRset, RGset, what=c("beta", "M"), cutp=0.01, cutsamples=0.95, cutcpgs=0.95, verbose=TRUE, ...) {

    what <- match.arg(what)
    
    if (!inherits(GRset, "GenomicRatioSet"))
        stop("First argument should be an object of class 'GenomicRatioSet' from preprocessFunnorm!")
    if (!inherits(RGset, "RGChannelSet") | !inherits(RGset, "RGChannelSetExtended"))
        stop("Second argument should be an object of class RGChannelSet or RGChannelSetExtended from probeFiltering!")

    ##Filter on detection P-value using minfi's detectionP
    if(verbose)
        cat("Calculate and filter on detection P-value... \n")

    pvalmat <- detectionP(RGset, na.rm=TRUE)
    idPvalmat <- pvalmat > cutp
    idPvalmat[is.na(idPvalmat)] <- TRUE ##set those for which a detection P-value could not be calculate TRUE to be filtered out
    ##pvalmat[idPvalmat] <- NA

    if(verbose)
        cat("On average", round(100*sum(idPvalmat)/prod(dim(idPvalmat)), 2),"% of the CpGs (",nrow(idPvalmat),") have detection P-value above the threshold ",cutp, "\n")

    if(verbose)
        cat("Transform to ",what,"-values... \n")

    if(what == "M"){
        matfilt <- getM(preprocessRaw(RGset), ...)
        matnorm <- getM(GRset, ...)
    }
    else {
        matfilt <- getBeta(RGset, ...)
        matnorm <- getBeta(GRset, ...)
    }

    ##set max/min M-values to +/- 16
    if(verbose & what=="M")
        cat("Set +/-Inf to +/-16... \n")
    
    if(what == "M")
        matnorm[!is.finite(matnorm)] <-  sign(matnorm[!is.finite(matnorm)])*16
    
    ##set NA from probeFiltering
    ##!!!NOTE orders are not the same
    if(verbose)
         cat("On average", round(100*sum(is.na(matfilt))/prod(dim(matfilt)), 2),"% of the probes (",nrow(matfilt),") were set to NA in the probe filtering step! \n")
    mid <- match(rownames(matfilt), rownames(matnorm))
    matnorm <- matnorm[mid,]
    matnorm[is.na(matfilt)] <- NA

    ##set NA from detectionP
    ##order seems OK just to be sure
    mid <- match(rownames(matnorm), rownames(idPvalmat))
    idPvalmat <- idPvalmat[mid,]
    matnorm[idPvalmat] <- NA

    ##Replaced by gap_hunting
    ##set chen CpGs/probes NA
    ##if(verbose)
    ##    message("Removing cross-reactive or polymorphic probes...")

    ##data("chen", package="Leiden450K")
    ##matnorm[rownames(matnorm) %in% names(chenProbes),] <- NA

    ##calculate success rates and reduce
    if(verbose)
        cat("Calculate success rates and reduce... \n")
    
    ##srCols <- apply(matnorm, 2, function(x) sum(!is.na(x))/(length(x) - 30969)) ##chen CpGs excluded
    srCols <- apply(matnorm, 2, function(x) sum(!is.na(x))/(length(x)))
    srRows <- apply(matnorm, 1, function(x) sum(!is.na(x))/length(x))

    if(verbose){
        cat("Percentage of samples having success rate above", cutsamples, "is", round(100*sum(srCols > cutsamples)/length(srCols),2),"% \n")
        cat("Percentage of CpGs having success rate above", cutcpgs, "is", round(100*sum(srRows > cutcpgs)/length(srRows),2),"% \n")
    }
    invisible(matnorm[srRows > cutcpgs,  srCols > cutsamples])
}
