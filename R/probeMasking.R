#' This function masks cross-reactive and polymorphic probes as specified
#' in Zhou et al. "Comprehensive characterization, annotation and innovative
#' use of Infinium DNA methylation BeadChip probes" (2017).
#'
#' Specification for both EPIC and HM450 arrays is possible and masking is
#' genome dependent; different probes are filtered depending
#' on whether hg38 or hg19 is specified as a reference.
#'
#' @title Cross-reactive and Polymorphic Probe Masking
#' @param values beta or M values
#' @param array EPIC or 450 array
#' @param genome hg19 or hg38 as reference
#' @param verbose default is TRUE
#' @return beta or M values
#' @import htm2txt
#' @export
#' @author ljsinke

probeMasking <- function(values, array=c("EPIC","450K"), genome=c("hg19","hg38"), verbose=TRUE){
  
  # testUrl <- substring(gettxt("http://zwdzwd.github.io/InfiniumAnnotation#current"),121,130)
  # if(testUrl != "Jul-4-2020"){
  #  warning("This function appears to be out-of-date. Please contact the maintainer.", immediate.=TRUE)
  # }
  
  if(verbose==TRUE) {
    cat("[probeFilterDNAmArray] Extracting probe filter... \n")
    numValues <- nrow(values)
  }
  
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
  
  maskProbes <- names(maskProbes[maskProbes$MASK_general])
  
  if(verbose==TRUE) {
    numMask <- length(maskProbes)
    cat("[DNAmArray]", numMask, "probes considered for filtering... \n")
  }
  
  values <- values[!(rownames(values) %in% maskProbes),]
  
  if(verbose==TRUE) {
    numGone <- numValues - nrow(values)
    cat("[DNAmArray]", numGone,"/",numValues, "probes removed from the dataset \n")
  }
  return(values)
}
