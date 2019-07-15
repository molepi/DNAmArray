##' DMR finder
##'
##' Algorithm to identify differentially methylated regions (DMRs)
##' from a list of differentially methylated probes (DMPs). This was
##' developed by Roderick Slieker and is described in "Identification 
##' and systematic annotation of tissue-specific differentially methylated 
##' regions using the Illumina 450k array" Slieker et al. (2013)
##' @title DMR finder
##' @param regions GRanges-class object, named vector, matrix, or data.frame
##' @param mismatches Number of non-significant DMPs allowed (default: 3)
##' @param icd Inter-CpG distance (default: 1000bp)
##' @return GRanges-class object containing DMRs
##' @author R Sliecker, E.W Lameijer, and M. van Iterson
##' @importFrom IRanges IRanges
##' @importFrom GenomeInfoDb keepSeqlevels
##' @export
##' @examples
##' dmrData <- data.frame(dmp = attr(padj.cate, "names"), crit = ifelse(padj.cate<0.05, 1, 0))
##' DMRs <- DMRfinder(dmrData, mismatches=3, icd=1000)

DMRfinder = function(data, chromosome=c(1:22,"X","Y"), mismatches=3, icd=1000, illumina=TRUE){
  
  #Check for annotation  
  if(illumina==TRUE){
    check  = length(grep("cg", data[,1]))> 1
    check2  = length(grep("0", data[,1]))> 1
  
    if(check==FALSE | check2==FALSE){
      print("Please check your data, is the format right?")
      }else{
        library("FDb.InfiniumMethylation.hg19")
        InfiniumMethylation <- features(FDb.InfiniumMethylation.hg19)
        probesselect= as.data.frame(data, stringsAsFactors = F)
        IM = InfiniumMethylation[names(InfiniumMethylation) %in% probesselect[,1],]
        IMx  =   cbind(as.character(seqnames(IM)),start(IM),end(IM))
        rownames(IMx) = names(IM)
        IMx = as.data.frame(IMx[probesselect[,1],],stringsAsFactors=F)
        probesselect$chr = IMx$V1
        probesselect$coord  = IMx$V2
        probesselect = probesselect[,c(3,4,2)]
        }
  }else{
    probesselect = data
  }
  
  pb = txtProgressBar(min=1,max=length(chromosome)+1,style=1)
  
  #DMR finder loop per chromosome
  DMR = function(chromosome,data,mismatches,icd){
	clist = which(c(1:22,"X","Y") %in% chromosome)
    setTxtProgressBar(pb,clist)
    close(pb)
    tryCatch.W.E <- function(expr)
    {
      W <- NULL
      w.handler <- function(w){ # warning handler
        W <<- w
        invokeRestart("muffleWarning")
      }
      list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                       warning = w.handler),
           warning = W)
    }
	chromosomex = paste("chr",chromosome,sep="")
    chr1 = probesselect[probesselect[,1]==chromosomex,]
	check.zero = sum(chr1[,3]==1) == 0
	if(check.zero==TRUE){
	}else{
    order = order(as.numeric(as.character(chr1[,2])))
    chr.sorted=chr1[order,]
    chr.final = data.frame(coord = as.numeric(chr.sorted[,2]) , crit = as.numeric(chr.sorted[,3]))
    file = "temp.txt"
    sink(file)
    MAXIMUM_REGION_LENGTH = icd # constant, can be adjusted
    last_coordinate = length( chr.final$crit )
    next_coordinate = 0 # so that a region that has been called will be skipped
    
    for (i in 1:(last_coordinate-1)) {
      if ( i>=next_coordinate ) {
        if (chr.final$crit[ i ]==1) {
          start_location = chr.final$coord[ i ]
          last_visited_crit_loc = start_location
          sum_of_ones = 1
          number_of_items = 1
          
          # start crawling loop
          for (j in (i+1):last_coordinate ) {
            if (chr.final$coord[ j ] > (last_visited_crit_loc + MAXIMUM_REGION_LENGTH)) { break }
            if((number_of_items-sum_of_ones)>mismatches) { break }   #Number of mismatches
            number_of_items = number_of_items + 1
            if (chr.final$crit[j]==1) { 
              last_visited_crit_loc = chr.final$coord[ j ]
              sum_of_ones = sum_of_ones + 1 
            }
          }
          
          # now check if the result is good enough
          if (sum_of_ones>=3) {
            last_one=i+number_of_items-1
            for (k in (i+number_of_items-1):1) {
              if ( chr.final$crit[k] == 0 ) {
                last_one = last_one - 1
                number_of_items = number_of_items - 1
              }
              else {
                break
              }
            }
            cat(start_location,";",chr.final$coord[last_one],";",sum_of_ones/number_of_items,"\n")
            next_coordinate = last_one + 1
          }
        }
      }
    }
    sink()
    check = tryCatch.W.E(read.table("temp.txt" , sep=";"))
    check.for.error = class(check$value)[1] == "simpleError"
    if(check.for.error==TRUE){
      
    }else{
      
      dmr = as.matrix(check$value)
      dmr.chr = cbind(rep(chromosomex),dmr[,1:2,drop=F])
      dmr.chr
    }
  }
  }  
  #Run DMR finder for selected chromosomes
  dmr.func = lapply(chromosome,DMR,probesselect,mismatches,icd)
  dmr.allchr = do.call(rbind , dmr.func)
  check.class = class(dmr.allchr) == "NULL"
  if(check.class==TRUE){
	  print("No DMRs identified!")
  }else{
  dmrs.retGR = GRanges(seqnames=dmr.allchr[,1],IRanges(as.numeric(dmr.allchr[,2]),as.numeric(dmr.allchr[,3])))
  dmrs.retGR
}

}

