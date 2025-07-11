##' Scree plot of the control probe intensities (450k or EPIC)
##'
##' Scree-plot to detect the optimal number of PCs to use in Functional Normalization
##' @title Scree plot of the control probe intensities
##' @param RGset object of class 'RGChannelSet'
##' @param nmax maximal number of PCs to show
##' @return output of prcomp using the control probe matrix
##' @importFrom stats prcomp
##' @importFrom graphics barplot
##' @export
screeplot <- function(RGset, nmax = 10) {
  controlMatrix <- minfi:::.buildControlMatrix(minfi:::.extractFromRGSet(RGset))
  pc <- prcomp(controlMatrix)
  
  nmax <- min(nmax, nrow(controlMatrix))
  barplot(summary(pc)$importance[2, 1:nmax], 
          ylab = "Proportion of Variance", 
          main = "Scree Plot", 
          col = "#7cb4c9")
  
  invisible(pc)
}

##' preprocessFunnorm from minfi with additional return options
##'
##' preprocessFunnorm from minfi with additional return options: keepCN and what
##' @title preprocessFunnormv2
##' @param rgSet object of class 'RGChannelSet'
##' @param nPCs number of PC to use for normalization
##' @param sex vector of sexes of each sample
##' @param bgCorr logical TRUE/FALSE
##' @param dyeCorr logical TRUE/FALSE
##' @param verbose logical TRUE/FALSE
##' @param keepCN logical TRUE/FALSE
##' @param ... optional arguments to ratioConvert
##' @return object of class 'GenomicRatioSet'
##' @author mvaniterson
##' @importFrom SummarizedExperiment updateObject assay assay<-
##' @export
preprocessFunnorm.DNAmArray <- function(rgSet, nPCs=2, sex = NULL, bgCorr = TRUE, dyeCorr = TRUE, verbose = TRUE, keepCN=TRUE, ...) {
    
    rgSet <- updateObject(rgSet) ## FIXM: might not KDH: technically, this should not be needed, but might be nice

    if(!keepCN)
        cat("[preprocessFunnorm] Modified version not returning CN! \n")
    
    # Background correction and dye bias normalization:
    if (bgCorr){
        if(verbose && dyeCorr) {
            cat("[preprocessFunnorm] Background and dye bias correction with noob \n") 
        } else {
            cat("[preprocessFunnorm] Background correction with noob \n") 
        }
        gmSet <- preprocessNoob(rgSet, dyeCorr = dyeCorr)
        if(verbose) cat("[preprocessFunnorm] Mapping to genome\n")
        gmSet <- mapToGenome(gmSet)
    } else {
        if(verbose) cat("[preprocessFunnorm] Mapping to genome\n")
        gmSet <- mapToGenome(rgSet)
    }
  
    subverbose <- max(as.integer(verbose) - 1L, 0)
    
    if(verbose) cat("[preprocessFunnorm] Quantile extraction\n")
    extractedData <- minfi:::.extractFromRGSet450k(rgSet)
    rm(rgSet)

    if (is.null(sex)) {
        gmSet <- addSex(gmSet, getSex(gmSet, cutoff = -3))
        sex <- rep(1L, length(gmSet$predictedSex))
        sex[gmSet$predictedSex == "F"] <- 2L
    }
    if(verbose) cat("[preprocessFunnorm] Normalization\n")

    if(keepCN)
        CN <- getCN(gmSet)
    
    gmSet <- minfi:::.normalizeFunnorm450k(object = gmSet, extractedData = extractedData,
                                            sex = sex, nPCs = nPCs, verbose = subverbose)
    
    grSet <- ratioConvert(gmSet, keepCN=keepCN, type="Illumina")
    
    if(keepCN)
        assay(grSet, "CN") <- CN
    
    grSet@preprocessMethod <- c(preprocessMethod(gmSet),
                                mu.norm = sprintf("Funnorm, nPCs=%s", nPCs))
    return(grSet)
 }	

