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

probeMasking <- function(values, array=c("EPIC","450"), genome=c("hg19","hg38"), verbose=TRUE){
  
  testUrl <- gettxt("http://zwdzwd.github.io/InfiniumAnnotation#current")
  cat(testUrl)
  
  if(verbose==TRUE) {
    cat("[probeFilterDNAmArray] Extracting probe filter... \n")
    numValues <- nrow(values)
  }
  
  if(array=="EPIC" & genome=="hg19") {
    maskProbes <- DNAmArray::maskEPIChg19
  }
  if(array=="EPIC" & genome=="hg38") {
    maskProbes <- DNAmArray::maskEPIChg38
  }
  if(array=="450" & genome=="hg19") {
    maskProbes <- DNAmArray::mask450Khg19
  }
  if(array=="450" & genome=="hg38") {
    maskProbes <- DNAmArray::mask450Khg38
  }
  
  if(verbose==TRUE) {
    numMask <- length(maskProbes)
    cat("[probeFilterDNAmArray]", numMask, "probes considered for filtering... \n")
  }
  
  values <- values[!(rownames(values) %in% maskProbes),]
  
  if(verbose==TRUE) {
    numGone <- numValues - nrow(values)
    cat("[probeFilterDNAmArray]", numGone,"/",numValues, "probes removed from the dataset \n")
  }
  return(values)
}
