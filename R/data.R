#' Dataframe with overlaps GoNL variants and 450K probes
#'
#' Dataframe containing all SNPs and short INDELS from GoNLv5 that
#' overlap with 450K probes.  This release does not include X and Y
#' chromosomes, so only information for autosomal probes is
#' available. For each overlap there is an unique row. Consequently,
#' some probes are duplicated (probes that overlap with multiple
#' variants) and some variants are duplicated (some variants overlap
#' with more than one probe).
#'
#'
#' @format A data frame with 207866 rows and 19 variables:
#' \describe{
#'   \item{CHROM}{chromosome, X and Y chromosomes are not available,
#'                since they are not included in this GoNL release}
#'   \item{probe}{probe ID}
#'   \item{type}{Infinium probedesign}
#'   \item{strand}{orientation of the probe}
#'   \item{probeType}{whether the probe measures a CpG-site (cg) or
#'                    a non-CpG site (ch)}
#'   \item{location_c}{Location of the queried 'C' of the CpG dinucleotide. Note
#'                     that this is the location of the C that is actually measured.
#'                     For probes that interrogate the reverse strand (plus-strand probes) this
#'                     is one base downstream of the C nucleotide on the forward strand}
#'   \item{location_g}{Location of the G nucleotide of the CpG dinucleotide. Note that this
#'                     is the location of the queried G. For probes that interrogate the reverse strand
#'                     (plus-strand probes) this is one base upstream of the G nucleotide on the forward strand}
#'   \item{ID}{SNP ID}
#'   \item{snpBeg}{Start coordinate of the variant. Identical to snpEnd for SNPs.}
#'   \item{snpEnd}{End coordinate of the variant. Identical to snpBeg for SNPs}
#'   \item{AF}{Allele frequency of alternative allele}
#'   \item{REF}{Reference allele}
#'   \item{ALT}{Alternative allele}
#'   \item{FILTER}{Filter information from GoNL.}
#'   \item{MAF}{Minor allele frequency}
#'   \item{variantType}{SNP or INDEL}
#'   \item{distance_3end}{Distance between SNP and 3'end of the probe. For type I probes
#'                        the 3'end of the probe coincides with the queried C nucleotide.
#'                        For type II probes the 3'end of the probe coincides with the G
#'                        nucleotide directly after the C nucleotide.}
#'   \item{distance_c}{Distance from queried C nucleotide. A distance of -1 indicates
#'                     that the SNPs overlaps the SBE-position for type I probes.}
#'   \item{channel_switch}{Indicates whether a variant in the SBE-location of type I probes
#'                         causes a color-channel-switch or overlap with an INDEL. For plus-strand probes C/T, C/A and C/G
#'                         SNPs are expected to cause a color-channel switch. For min-strand probes
#'                         A/G, G/T and C/G SNPs are expected to cause a color-channel switch.}
#'                     }
#'
#' @usage data(hg19.GoNLsnps)
#' 
#' @examples
#'     data(hg19.GoNLsnps)
#'     
#'     # Select variants that overlap with queried C nucleotide
#'     snps_c <- hg19.GoNLsnps[hg19.GoNLsnps$distance_c == 0, ] 
#'     
#'     # Select all INDELS
#'     indels <- hg19.GoNLsnps[hg19.GoNLsnps$variantType == "INDEL",] 
#'     
#'     # Select SNPs that cause a channel-switch
#'     channel_switch <- hg19.GoNLsnps[!is.na(hg19.GoNLsnps$channel_switch)
#' & hg19.GoNLsnps$channel_switch == "Yes",]
#'
#' @source 
#'         \url{http://zwdzwd.github.io/InfiniumAnnotation}
#'         
#'         \url{https://molgenis26.target.rug.nl/downloads/gonl_public/variants/release5/}
"hg19.GoNLsnps"

#' HM450 population-specific probe-masking recommendations
#' 
#' Adapted version of the annotation file provided by Zhou et al. (see
#' source, Mar-13-2017 release).  This annotation file contains
#' population-specific probe-masking recommendations based on SNPs
#' within 5 bases from the 3'end of the probe, mapping issues,
#' non-unique 3' 30bp subsequence and channel-switching SNPs in the
#' single-base-extension for type I probes. We added
#' population-specific masking recommendations for the Dutch
#' population using GoNL release 5. This release does not include X
#' and Y chromosomes, so for the Dutch population, only masking
#' information for the autosomal probes is available.
#' 
#' Note: Zhou et al. identified several probes that match to a
#' different location than annotated in the original Illumina manifest
#' file. The authors have used the 'updated' location in their
#' annotation file. Therefore, a handful of probes in this annotation
#' file are annotated to a different location than the original
#' Illumina manifest file. For the identification of the overlaps we
#' used the locations as annotated in the original Illumina file.
#' Therefore, a few of the identified overlaps do not match the
#' locations as specified in this annotation file.  This also explains
#' why GoNL masking information is available for a couple of probes
#' that are located in the X- and Y-chromosome in this annotation:
#' these probes map to autosomal probes in the original Illumina
#' file. All probes that map to a different location than originally
#' annotated are recommended to be masked (in the MASK.mapping and
#' MASK.general column), so generally they won't be included in
#' further analyses.
#'
#' @format A GRanges object with 485577 ranges and 65 metadata columns:
#' \describe{
#'   \item{MASK.general.pop}{Recommended general purpose masking merged from "MASK.sub30.copy", "MASK.mapping" (in either
#'         the hg38 or hg19 genome), "MASK.extBase", "MASK.typeINextBaseSwitch" and "MASK.snp5.pop" from the 
#'         "hm450.manifest" file and the "hm450.manifest.pop" file (see source).
#'         For GoNL, "MASK.typeINextBaseSwitchandINDEL.GoNL" is used instead of "MASK.typeINextBaseSwitch"}
#'   \item{MASK.snp5.pop}{Whether the 5bp 3'-subsequence (including extension for type II) overlap with a SNP with
#'                        population-specific AF > 0.01}
#'   \item{MASK.typeINextBaseSwitchandINDEL.GoNL}{SNPs (that cause a color-channel switch) and INDELS with AF > 0.01 in GoNL. 
#'         In contrast, "MASK.typeINextBaseSwitch" column is based on all SNPs in 1000 genomes and dbSNP,
#'         regardless of population or allele frequency}
#' }
#'
#' @source \url{http://zwdzwd.github.io/InfiniumAnnotation}
#' 
#' @usage data(hm450.manifest.pop.GoNL)
#' 
#' @examples 
#' # Select probes that should be masked in Dutch population (note that X and Y chromosomes are not included)
#' hm450.manifest.pop.GoNL <- hm450.manifest.pop.GoNL[!is.na(hm450.manifest.pop.GoNL$MASK.general.GoNL) &
#'                                                      hm450.manifest.pop.GoNL$MASK.general.GoNL == TRUE, ]
#'   
#' # Select probes that should be masked in Dutch population because there is 
#' # a SNP within 5 bases of the 3'end of the probe (note that X and Y chromosomes are not included)
#' hm450.manifest.pop.GoNL <- hm450.manifest.pop.GoNL[!is.na(hm450.manifest.pop.GoNL$MASK.snp5.GoNL) &
#'                                                        hm450.manifest.pop.GoNL$MASK.snp5.GoNL == TRUE, ]
#' # When studying a Dutch population and one wants to include X and Y chromosomal probes,
#' # the EUR or CEU population can be used.                                                      
#' # Select probes that should be masked in European population (these include X and Y chromosomes)                                                      
#' hm450.manifest.pop.GoNL <- hm450.manifest.pop.GoNL[hm450.manifest.pop.GoNL$MASK.general.EUR == TRUE,]
#'
#'
#' @references
#'        Zhou W, Laird PW and Shen H: Comprehensive characterization, annotation and innovative use of Infinium DNA Methylation BeadChip probes.
#'        Nucleic Acids Research 2016
"hm450.manifest.pop.GoNL"
