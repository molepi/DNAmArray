##DMR finder
.DMR <- function(regions, mismatches, icd){

    coord <- start(regions)
    crit <- mcols(regions)$crit
    crit <- crit[order(coord)]
    coord <- coord[order(coord)]
   
    MAXIMUM_REGION_LENGTH <- icd
    last_coordinate <- length(regions)
    next_coordinate <- 0 # so that a region that has been called will be skipped
    current_row <- 0
    for (i in 1:(last_coordinate-1)) {
        if ( i >= next_coordinate ) {
            if (crit[i] == 1) {
                start_location <- coord[i]
                last_visited_crit_loc <- start_location
                sum_of_ones <- 1
                number_of_items <- 1

                ## start crawling loop
                for (j in (i+1):last_coordinate ) {
                    if (coord[j] > (last_visited_crit_loc + MAXIMUM_REGION_LENGTH)) {
                        break
                    }
                    if((number_of_items-sum_of_ones) > mismatches) {
                        break
                    }   #Number of mismatches
                    number_of_items <- number_of_items + 1
                    if (crit[j] == 1) {
                        last_visited_crit_loc <- coord[j]
                        sum_of_ones <- sum_of_ones + 1
                    }
                }
                ## now check if the result is good enough
                if (sum_of_ones >= 3) {
                    last_one <- i + number_of_items-1
                    for (k in (i + number_of_items-1):1) {
                        if ( crit[k] == 0 ) {
                            last_one <- last_one - 1
                            number_of_items <- number_of_items - 1
                        }else{
                            break
                        }
                    }
                    current_row <- current_row + 1
                    out <- rbind(out, c(start_location, coord[last_one], signif(100*sum_of_ones/number_of_items))) ##slow and memory inefficient
                    next_coordinate <- last_one + 1
                }
            }
        }
    }

    if(is.null(out))
        return(NULL)
    else {
        out <- out[order(out[,1]), ]
        return(invisible(GRanges(seqnames=unique(as.character(seqnames(regions))),
                                 IRanges(start=out[,1], end=out[,2]),
                                 percentage=out[,3])))
    }
}

##' DMR finder
##'
##' Algorithm is described in Identification and systematic annotation
##' of tissue-specific differentially methylated regions using the
##' Illumina 450k array Slieker et al. Epigentics and Chromatin 2013,
##' 6,26.
##' @title DMR finder
##' @param regions GRanges-object, named vector or matrix/data.frame
##' @param mismatches number of non-significant DMPs default <= 3
##' @param icd inter-CpG distance default 1000bp
##' @return GRanges object with DMRs
##' @author R Sliecker, E.W Lameijer and M. van Iterson
##' @importFrom IRanges IRanges
##' @importFrom GenomeInfoDb keepSeqlevels
##' @export
##' @examples
##' require(FDb.InfiniumMethylation.hg19)
##' feats <- features(FDb.InfiniumMethylation.hg19)
##' regions <- feats[seqnames(feats) %in% c("chr21", "chr22")]
##' mcols(regions)  <- NULL ##drop unnecessary metadata
##' mcols(regions)$crit <- rbinom(length(regions), 1, prob=0.2)
##'
##' ##input GRanges
##' (dmrsR <- DMRfinder(regions))
##'
##' ##input named-vector
##' dmps <- mcols(regions)$crit
##' names(dmps) <- names(regions)
##' (dmrsV <- DMRfinder(dmps))
##'
##' ##bed-like data.frame
##' regions <- as.data.frame(regions)
##' regions <- regions[, c("seqnames", "start", "crit")]
##' colnames(regions) <- c("chr", "pos", "crit")
##' (dmrsD <- DMRfinder(regions))
##'
##' ##all equal
##' length(dmrsR)
##' length(dmrsV)
##' length(dmrsD)
##'
##' head(dmrsR)
##' head(dmrsV)
##' head(dmrsD)
##'
##' dmrsR == dmrsV
##' dmrsR == dmrsD
##' dmrsV == dmrsD
DMRfinder <- function(regions, mismatches=3, icd=1000){

    if(is.vector(regions)){
        requireNamespace("FDb.InfiniumMethylation.hg19")
        feats <- features(FDb.InfiniumMethylation.hg19)
        mid <- match(names(regions), names(feats))
        feats <- feats[mid]
        mcols(feats)$crit <- as.numeric(regions)
        regions <- feats
    }
    else if(is.data.frame(regions) | is.matrix(regions)) {
        if(any(!(colnames(regions) %in% c("end"))))
            regions$end <- regions[, which(colnames(regions) %in% c("start","pos", "bp"))]
        regions <- makeGRangesFromDataFrame(regions, start.field=c("start","pos", "bp"),
                                            keep.extra.columns=TRUE,
                                            ignore.strand=TRUE)
    }

    stopifnot(class(regions) != "GRanges" | class(regions) != "GRangesList")

    if(!any(colnames(mcols(regions)) %in% "crit"))
        stop("No `crit` column provided with significant DMPs!")

    ##dropping used chromosomes before splitting on chromosome
    regions <- keepSeqlevels(regions, unique(as.character(seqnames(regions))))
    regions <- split(regions, seqnames(regions))

    ##apply DMRfinder on each chromosome separately
    out <- lapply(regions, .DMR, mismatches, icd)
    invisible(unlist(GRangesList(out)))
}
