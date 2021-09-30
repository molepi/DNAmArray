# TSV file with hg19 mask information for EPIC probes
#
# Table containing all SNPs and short INDELS that overlap with 
# EPIC probes.For each overlap there is an unique row. Consequently,
# some probes are duplicated (probes that overlap with multiple
# variants) and some variants are duplicated (some variants overlap
# with more than one probe).
#
#
# @format A data frame with 865918 rows and 57 variables:
# \describe{
#   \item{CpG_chrm}{chromosome location of the target}
#   \item{CpG_beg}{0-based co-ordinate of the target. The 
#   co-ordinates should have a span of 2 nucleotides for
#   CpG probes, or 1 nucleotide for CpH and SNP probes. Some
#   erroneous CpH probe co-ordinates mapping information in 
#   the manufacturer's manifest have been corrected.}
#   \item{CpG_end}{1-based co-ordinate of the target.The 
#   co-ordinates should have a span of 2 nucleotides for
#   CpG probes, or 1 nucleotide for CpH and SNP probes. Some
#   erroneous CpH probe co-ordinates mapping information in 
#   the manufacturer's manifest have been corrected.}
#   \item{probe_strand}{strand orientation of the actual
#   probe. '+' is for all the up-probes positioned in 
#   smaller co-ordinates and '-' for all the down-probes
#   positioned in greater co-ordinates with respect to
#   the target CpGs. '*' is used for unmapped probes.}
#   \item{probeID}{probe ID}
#   \item{address_A}{address of probe A on the chip
#   designated by the original manifest}
#   \item{address_B}{address of probe B on the chip
#   designated by the original manifest}
#   \item{channel}{'Both' for type II probes and 'Grn' or 
#   'Red' for Type I probes}
#   \item{designType}{Type of probe, either 'I' or 'II'}
#   \item{nextBase}{the actual extension base (on the probe
#   strand) after bisulfite conversion ('A' or 'C' or 'T').
#   Unmapped probes have extension base labeled in the original
#   manifest.}
#   \item{nextBaseRef}{the extension base (on the hybridized /
#   template DNA) before bisulfite conversion ('A', 'C', 'G',
#   or 'T'). Unmapped probes have 'NA'.}
#   \item{probeType}{either 'cg', 'ch', or 'rs'}
#   \item{orientation}{either 'up or 'down', specifying 
#   whether the probe is positioned upstream (in smaller
#   co-ordinates) or downstream (in greater co-ordinates)
#   from the target}
#   \item{probeCpGcnt}{the number of additional CpGs in the
#   probe (not counting the interrogated CpG)}
#   \item{context35}{the number of CpGs in the [-35bp, +35bp]
#   window}
#   \item{probeBeg}{the mapped start position of the probe,
#   which is always 50bp long}
#   \item{probeEnd}{the mapped end position of the probe,
#   which is always 50bp long}
#   \item{ProbeSeq_A}{the probe sequence for allele A}
#   \item{ProbeSeq_B}{the probe sequence for allele B}
#   \item{gene}{comma separated list of gene annotations
#   (unique and alphabetically sorted). Gene models follow
#   GENCODE version 22 (hg38)}
#   \item{gene_HGNC}{comma separated list of gene annotations
#   (unique and alphabetically sorted). Genes are checked
#   using HGNChelper for compatibility with HGNC. Gene models
#   follows GENCODE version 22 (hg38)}
#   \item{chrm_A}{the mapping info for probe A excluding 
#   decoy chromosomes}
#   \item{beg_A}{the mapping info for probe A excluding 
#   decoy chromosomes}
#   \item{flag_A}{the mapping info for probe A excluding 
#   decoy chromosomes}
#  \item{mapQ_A}{the mapping quality score for probe A excluding 
#   decoy chromosomes, with 60 being the best}
#   \item{cigar_A}{the mapping info for probe A excluding 
#   decoy chromosomes}
#   \item{NM_A}{the mapping info for probe A excluding 
#   decoy chromosomes}
#   \item{chrm_B}{the mapping info for probe B excluding 
#   decoy chromosomes}
#   \item{beg_B}{the mapping info for probe B excluding 
#   decoy chromosomes}
#   \item{flag_B}{the mapping info for probe B excluding 
#   decoy chromosomes}
#   \item{mapQ_B}{the mapping quality score for probe B excluding 
#   decoy chromosomes, with 60 being the best}
#   \item{cigar_B}{the mapping info for probe B excluding 
#   decoy chromosomes}
#   \item{NM_B}{the mapping info for probe B excluding 
#   decoy chromosomes}
#   \item{wDecoy_chrm_A}{the mapping info for probe A including 
#   decoy chromosomes}
#   \item{wDecoy_beg_A}{the mapping info for probe A including 
#   decoy chromosomes}
#   \item{wDecoy_flag_A}{the mapping info for probe A including 
#   decoy chromosomes}
#   \item{wDecoy_mapQ_A}{the mapping quality score for probe A including 
#   decoy chromosomes, with 60 being the best}
#   \item{wDecoy_cigar_A}{the mapping info for probe A including 
#   decoy chromosomes}
#   \item{wDecoy_NM_A}{the mapping info for probe A including 
#   decoy chromosomes}
#   \item{wDecoy_chrm_B}{the mapping info for probe B including 
#   decoy chromosomes}
#   \item{wDecoy_beg_B}{the mapping info for probe B including 
#   decoy chromosomes}
#   \item{wDecoy_flag_B}{the mapping info for probe B including 
#   decoy chromosomes}
#   \item{wDecoy_mapQ_B}{the mapping quality score for probe B including 
#   decoy chromosomes, with 60 being the best}
#   \item{wDecoy_cigar_B}{the mapping info for probe B including 
#   decoy chromosomes}
#   \item{wDecoy_NM_B}{the mapping info for probe B including 
#   decoy chromosomes}
#   \item{posMatch}{whether the mapping matches the original
#   manifest, it only applies to hg19 and will be NA under 
#   hg38}
#   \item{MASK_mapping}{whether the probe is masked for mapping
#   reasons. Probes retained should have high quality (>=40)
#   consistent (with designed MAPINFO) mapping (for both in the
#   cast of type I) without INDELS}
#   \item{MASK_typeINextBaseSwitch}{whether the probe has a SNP
#   in the extension base that causes a color channel switch from
#   the official annotation (described as color channel switching
#   or CCS SNP in the reference). These probes should be processed
#   differently than designed (by summin up both color channels
#   instead of just the annotated color channel)}
#   \item{MASK_rmsk15}{whetehr the 15bp 3'-subsequence of the
#   probe overlaps with repeat masker, this MASK is NOT 
#   recommended}
#   \item{MASK_sub40_copy}{whether the 40bp 3'subsequence of the
#   probe is non-unique}
#   \item{MASK_sub35_copy}{whether the 35bp 3'subsequence of the
#   probe is non-unique}
#   \item{MASK_sub30_copy}{whether the 30bp 3'subsequence of the
#   probe is non-unique}
#   \item{MASK_sub25_copy}{whether the 25bp 3'subsequence of the
#   probe is non-unique}
#   \item{MASK_snp5_common}{whether 5bp 3'-subsequence (including
#   extension for tyep II) overlaps with any of the common SNPs
#   from dbSNP (global MAF can be under 1%)}
#   \item{MASK_snp5_GMAF1p}{whether 5bp 3'-subsequence (including
#   extension for type II) overlaps with any of the SNPs with 
#   global MAF > 1%}
#   \item{MASK_extBase}{probes masked for extension base inconsistent
#   with specified color channel (type I) or CpG (Type II) based
#   on mapping}
#   \item{MASK_general}{recommended general ppurpose masking
#   merged from MASK_sub30_copy, MASK_mapping, MASK_extBase,
#   MASK_typeINExtBaseSwitch and MASK_snp5_GMAF1p}
#
# @usage read_tsv('data/EPIC.hg19.manifest.tsv.gz')
#
# @source 
#         \url{http://zwdzwd.github.io/InfiniumAnnotation}
#         
"EPIC.hg19.manifest.txt"

# TSV file with hg38 mask information for EPIC probes
#
# Table containing all SNPs and short INDELS that overlap with 
# EPIC probes.For each overlap there is an unique row. Consequently,
# some probes are duplicated (probes that overlap with multiple
# variants) and some variants are duplicated (some variants overlap
# with more than one probe).
#
#
# @format A data frame with 865918 rows and 57 variables:
# \describe{
#   \item{CpG_chrm}{chromosome location of the target}
#   \item{CpG_beg}{0-based co-ordinate of the target. The 
#   co-ordinates should have a span of 2 nucleotides for
#   CpG probes, or 1 nucleotide for CpH and SNP probes. Some
#   erroneous CpH probe co-ordinates mapping information in 
#   the manufacturer's manifest have been corrected.}
#   \item{CpG_end}{1-based co-ordinate of the target.The 
#   co-ordinates should have a span of 2 nucleotides for
#   CpG probes, or 1 nucleotide for CpH and SNP probes. Some
#   erroneous CpH probe co-ordinates mapping information in 
#   the manufacturer's manifest have been corrected.}
#   \item{probe_strand}{strand orientation of the actual
#   probe. '+' is for all the up-probes positioned in 
#   smaller co-ordinates and '-' for all the down-probes
#   positioned in greater co-ordinates with respect to
#   the target CpGs. '*' is used for unmapped probes.}
#   \item{probeID}{probe ID}
#   \item{address_A}{address of probe A on the chip
#   designated by the original manifest}
#   \item{address_B}{address of probe B on the chip
#   designated by the original manifest}
#   \item{channel}{'Both' for type II probes and 'Grn' or 
#   'Red' for Type I probes}
#   \item{designType}{Type of probe, either 'I' or 'II'}
#   \item{nextBase}{the actual extension base (on the probe
#   strand) after bisulfite conversion ('A' or 'C' or 'T').
#   Unmapped probes have extension base labeled in the original
#   manifest.}
#   \item{nextBaseRef}{the extension base (on the hybridized /
#   template DNA) before bisulfite conversion ('A', 'C', 'G',
#   or 'T'). Unmapped probes have 'NA'.}
#   \item{probeType}{either 'cg', 'ch', or 'rs'}
#   \item{orientation}{either 'up or 'down', specifying 
#   whether the probe is positioned upstream (in smaller
#   co-ordinates) or downstream (in greater co-ordinates)
#   from the target}
#   \item{probeCpGcnt}{the number of additional CpGs in the
#   probe (not counting the interrogated CpG)}
#   \item{context35}{the number of CpGs in the [-35bp, +35bp]
#   window}
#   \item{probeBeg}{the mapped start position of the probe,
#   which is always 50bp long}
#   \item{probeEnd}{the mapped end position of the probe,
#   which is always 50bp long}
#   \item{ProbeSeq_A}{the probe sequence for allele A}
#   \item{ProbeSeq_B}{the probe sequence for allele B}
#   \item{gene}{comma separated list of gene annotations
#   (unique and alphabetically sorted). Gene models follow
#   GENCODE version 22 (hg38)}
#   \item{gene_HGNC}{comma separated list of gene annotations
#   (unique and alphabetically sorted). Genes are checked
#   using HGNChelper for compatibility with HGNC. Gene models
#   follows GENCODE version 22 (hg38)}
#   \item{chrm_A}{the mapping info for probe A excluding 
#   decoy chromosomes}
#   \item{beg_A}{the mapping info for probe A excluding 
#   decoy chromosomes}
#   \item{flag_A}{the mapping info for probe A excluding 
#   decoy chromosomes}
#   \item{mapQ_A}{the mapping quality score for probe A excluding 
#   decoy chromosomes, with 60 being the best}
#   \item{cigar_A}{the mapping info for probe A excluding 
#   decoy chromosomes}
#   \item{NM_A}{the mapping info for probe A excluding 
#   decoy chromosomes}
#   \item{chrm_B}{the mapping info for probe B excluding 
#   decoy chromosomes}
#   \item{beg_B}{the mapping info for probe B excluding 
#   decoy chromosomes}
#   \item{flag_B}{the mapping info for probe B excluding 
#   decoy chromosomes}
#   \item{mapQ_B}{the mapping quality score for probe B excluding 
#   decoy chromosomes, with 60 being the best}
#   \item{cigar_B}{the mapping info for probe B excluding 
#   decoy chromosomes}
#   \item{NM_B}{the mapping info for probe B excluding 
#   decoy chromosomes}
#   \item{wDecoy_chrm_A}{the mapping info for probe A including 
#   decoy chromosomes}
#   \item{wDecoy_beg_A}{the mapping info for probe A including 
#   decoy chromosomes}
#   \item{wDecoy_flag_A}{the mapping info for probe A including 
#   decoy chromosomes}
#   \item{wDecoy_mapQ_A}{the mapping quality score for probe A including 
#   decoy chromosomes, with 60 being the best}
#   \item{wDecoy_cigar_A}{the mapping info for probe A including 
#   decoy chromosomes}
#   \item{wDecoy_NM_A}{the mapping info for probe A including 
#   decoy chromosomes}
#   \item{wDecoy_chrm_B}{the mapping info for probe B including 
#   decoy chromosomes}
#   \item{wDecoy_beg_B}{the mapping info for probe B including 
#   decoy chromosomes}
#   \item{wDecoy_flag_B}{the mapping info for probe B including 
#   decoy chromosomes}
#   \item{wDecoy_mapQ_B}{the mapping quality score for probe B including 
#   decoy chromosomes, with 60 being the best}
#   \item{wDecoy_cigar_B}{the mapping info for probe B including 
#   decoy chromosomes}
#   \item{wDecoy_NM_B}{the mapping info for probe B including 
#   decoy chromosomes}
#   \item{posMatch}{whether the mapping matches the original
#   manifest, it only applies to hg19 and will be NA under 
#   hg38}
#   \item{MASK_mapping}{whether the probe is masked for mapping
#   reasons. Probes retained should have high quality (>=40)
#   consistent (with designed MAPINFO) mapping (for both in the
#   cast of type I) without INDELS}
#   \item{MASK_typeINextBaseSwitch}{whether the probe has a SNP
#   in the extension base that causes a color channel switch from
#   the official annotation (described as color channel switching
#   or CCS SNP in the reference). These probes should be processed
#   differently than designed (by summin up both color channels
#   instead of just the annotated color channel)}
#   \item{MASK_rmsk15}{whetehr the 15bp 3'-subsequence of the
#   probe overlaps with repeat masker, this MASK is NOT 
#   recommended}
#   \item{MASK_sub40_copy}{whether the 40bp 3'subsequence of the
#   probe is non-unique}
#   \item{MASK_sub35_copy}{whether the 35bp 3'subsequence of the
#   probe is non-unique}
#   \item{MASK_sub30_copy}{whether the 30bp 3'subsequence of the
#   probe is non-unique}
#   \item{MASK_sub25_copy}{whether the 25bp 3'subsequence of the
#   probe is non-unique}
#   \item{MASK_snp5_common}{whether 5bp 3'-subsequence (including
#   extension for tyep II) overlaps with any of the common SNPs
#   from dbSNP (global MAF can be under 1%)}
#   \item{MASK_snp5_GMAF1p}{whether 5bp 3'-subsequence (including
#   extension for type II) overlaps with any of the SNPs with 
#   global MAF > 1%}
#   \item{MASK_extBase}{probes masked for extension base inconsistent
#   with specified color channel (type I) or CpG (Type II) based
#   on mapping}
#   \item{MASK_general}{recommended general ppurpose masking
#   merged from MASK_sub30_copy, MASK_mapping, MASK_extBase,
#   MASK_typeINExtBaseSwitch and MASK_snp5_GMAF1p}
#
# @usage read_tsv('data/EPIC.hg19.manifest.tsv.gz')
#
# @source 
#         \url{http://zwdzwd.github.io/InfiniumAnnotation}
#         
NULL

# TSV file with hg19 mask information for 450K probes
#
# Table containing all SNPs and short INDELS that overlap with 
# EPIC probes.For each overlap there is an unique row. Consequently,
# some probes are duplicated (probes that overlap with multiple
# variants) and some variants are duplicated (some variants overlap
# with more than one probe).
#
#
# @format A data frame with 485577 rows and 57 variables:
# \describe{
#   \item{CpG_chrm}{chromosome location of the target}
#   \item{CpG_beg}{0-based co-ordinate of the target. The 
#   co-ordinates should have a span of 2 nucleotides for
#   CpG probes, or 1 nucleotide for CpH and SNP probes. Some
#   erroneous CpH probe co-ordinates mapping information in 
#   the manufacturer's manifest have been corrected.}
#   \item{CpG_end}{1-based co-ordinate of the target.The 
#   co-ordinates should have a span of 2 nucleotides for
#   CpG probes, or 1 nucleotide for CpH and SNP probes. Some
#   erroneous CpH probe co-ordinates mapping information in 
#   the manufacturer's manifest have been corrected.}
#   \item{probe_strand}{strand orientation of the actual
#   probe. '+' is for all the up-probes positioned in 
#   smaller co-ordinates and '-' for all the down-probes
#   positioned in greater co-ordinates with respect to
#   the target CpGs. '*' is used for unmapped probes.}
#   \item{probeID}{probe ID}
#   \item{address_A}{address of probe A on the chip
#   designated by the original manifest}
#   \item{address_B}{address of probe B on the chip
#   designated by the original manifest}
#   \item{channel}{'Both' for type II probes and 'Grn' or 
#   'Red' for Type I probes}
#   \item{designType}{Type of probe, either 'I' or 'II'}
#   \item{nextBase}{the actual extension base (on the probe
#   strand) after bisulfite conversion ('A' or 'C' or 'T').
#   Unmapped probes have extension base labeled in the original
#   manifest.}
#   \item{nextBaseRef}{the extension base (on the hybridized /
#   template DNA) before bisulfite conversion ('A', 'C', 'G',
#   or 'T'). Unmapped probes have 'NA'.}
#   \item{probeType}{either 'cg', 'ch', or 'rs'}
#   \item{orientation}{either 'up or 'down', specifying 
#   whether the probe is positioned upstream (in smaller
#   co-ordinates) or downstream (in greater co-ordinates)
#   from the target}
#   \item{probeCpGcnt}{the number of additional CpGs in the
#   probe (not counting the interrogated CpG)}
#   \item{context35}{the number of CpGs in the [-35bp, +35bp]
#   window}
#   \item{probeBeg}{the mapped start position of the probe,
#   which is always 50bp long}
#   \item{probeEnd}{the mapped end position of the probe,
#   which is always 50bp long}
#   \item{ProbeSeq_A}{the probe sequence for allele A}
#   \item{ProbeSeq_B}{the probe sequence for allele B}
#   \item{gene}{comma separated list of gene annotations
#   (unique and alphabetically sorted). Gene models follow
#   GENCODE version 22 (hg38)}
#   \item{gene_HGNC}{comma separated list of gene annotations
#   (unique and alphabetically sorted). Genes are checked
#   using HGNChelper for compatibility with HGNC. Gene models
#   follows GENCODE version 22 (hg38)}
#   \item{chrm_A}{the mapping info for probe A excluding 
#   decoy chromosomes}
#   \item{beg_A}{the mapping info for probe A excluding 
#   decoy chromosomes}
#   \item{flag_A}{the mapping info for probe A excluding 
#   decoy chromosomes}
#   \item{mapQ_A}{the mapping quality score for probe A excluding 
#   decoy chromosomes, with 60 being the best}
#   \item{cigar_A}{the mapping info for probe A excluding 
#   decoy chromosomes}
#   \item{NM_A}{the mapping info for probe A excluding 
#   decoy chromosomes}
#   \item{chrm_B}{the mapping info for probe B excluding 
#   decoy chromosomes}
#   \item{beg_B}{the mapping info for probe B excluding 
#   decoy chromosomes}
#   \item{flag_B}{the mapping info for probe B excluding 
#   decoy chromosomes}
#   \item{mapQ_B}{the mapping quality score for probe B excluding 
#   decoy chromosomes, with 60 being the best}
#   \item{cigar_B}{the mapping info for probe B excluding 
#   decoy chromosomes}
#   \item{NM_B}{the mapping info for probe B excluding 
#   decoy chromosomes}
#   \item{wDecoy_chrm_A}{the mapping info for probe A including 
#   decoy chromosomes}
#   \item{wDecoy_beg_A}{the mapping info for probe A including 
#   decoy chromosomes}
#   \item{wDecoy_flag_A}{the mapping info for probe A including 
#   decoy chromosomes}
#   \item{wDecoy_mapQ_A}{the mapping quality score for probe A including 
#   decoy chromosomes, with 60 being the best}
#   \item{wDecoy_cigar_A}{the mapping info for probe A including 
#   decoy chromosomes}
#   \item{wDecoy_NM_A}{the mapping info for probe A including 
#   decoy chromosomes}
#   \item{wDecoy_chrm_B}{the mapping info for probe B including 
#   decoy chromosomes}
#   \item{wDecoy_beg_B}{the mapping info for probe B including 
#   decoy chromosomes}
#   \item{wDecoy_flag_B}{the mapping info for probe B including 
#   decoy chromosomes}
#   \item{wDecoy_mapQ_B}{the mapping quality score for probe B including 
#   decoy chromosomes, with 60 being the best}
#   \item{wDecoy_cigar_B}{the mapping info for probe B including 
#   decoy chromosomes}
#   \item{wDecoy_NM_B}{the mapping info for probe B including 
#   decoy chromosomes}
#   \item{posMatch}{whether the mapping matches the original
#   manifest, it only applies to hg19 and will be NA under 
#   hg38}
#   \item{MASK_mapping}{whether the probe is masked for mapping
#   reasons. Probes retained should have high quality (>=40)
#   consistent (with designed MAPINFO) mapping (for both in the
#   cast of type I) without INDELS}
#   \item{MASK_typeINextBaseSwitch}{whether the probe has a SNP
#   in the extension base that causes a color channel switch from
#   the official annotation (described as color channel switching
#   or CCS SNP in the reference). These probes should be processed
#   differently than designed (by summin up both color channels
#   instead of just the annotated color channel)}
#   \item{MASK_rmsk15}{whetehr the 15bp 3'-subsequence of the
#   probe overlaps with repeat masker, this MASK is NOT 
#   recommended}
#   \item{MASK_sub40_copy}{whether the 40bp 3'subsequence of the
#   probe is non-unique}
#   \item{MASK_sub35_copy}{whether the 35bp 3'subsequence of the
#   probe is non-unique}
#   \item{MASK_sub30_copy}{whether the 30bp 3'subsequence of the
#   probe is non-unique}
#   \item{MASK_sub25_copy}{whether the 25bp 3'subsequence of the
#   probe is non-unique}
#   \item{MASK_snp5_common}{whether 5bp 3'-subsequence (including
#   extension for tyep II) overlaps with any of the common SNPs
#   from dbSNP (global MAF can be under 1%)}
#   \item{MASK_snp5_GMAF1p}{whether 5bp 3'-subsequence (including
#   extension for type II) overlaps with any of the SNPs with 
#   global MAF > 1%}
#   \item{MASK_extBase}{probes masked for extension base inconsistent
#   with specified color channel (type I) or CpG (Type II) based
#   on mapping}
#   \item{MASK_general}{recommended general ppurpose masking
#   merged from MASK_sub30_copy, MASK_mapping, MASK_extBase,
#   MASK_typeINExtBaseSwitch and MASK_snp5_GMAF1p}
#
# @usage read_tsv('data/EPIC.hg19.manifest.tsv.gz')
#
# @source 
#         \url{http://zwdzwd.github.io/InfiniumAnnotation}
#         
NULL


# TSV file with hg38 mask information for 450K probes
#
# Table containing all SNPs and short INDELS that overlap with 
# EPIC probes.For each overlap there is an unique row. Consequently,
# some probes are duplicated (probes that overlap with multiple
# variants) and some variants are duplicated (some variants overlap
# with more than one probe).
#
#
# @format A data frame with 485577 rows and 57 variables:
# \describe{
#   \item{CpG_chrm}{chromosome location of the target}
#   \item{CpG_beg}{0-based co-ordinate of the target. The 
#   co-ordinates should have a span of 2 nucleotides for
#   CpG probes, or 1 nucleotide for CpH and SNP probes. Some
#   erroneous CpH probe co-ordinates mapping information in 
#   the manufacturer's manifest have been corrected.}
#   \item{CpG_end}{1-based co-ordinate of the target.The 
#   co-ordinates should have a span of 2 nucleotides for
#   CpG probes, or 1 nucleotide for CpH and SNP probes. Some
#   erroneous CpH probe co-ordinates mapping information in 
#   the manufacturer's manifest have been corrected.}
#   \item{probe_strand}{strand orientation of the actual
#   probe. '+' is for all the up-probes positioned in 
#   smaller co-ordinates and '-' for all the down-probes
#   positioned in greater co-ordinates with respect to
#   the target CpGs. '*' is used for unmapped probes.}
#   \item{probeID}{probe ID}
#   \item{address_A}{address of probe A on the chip
#   designated by the original manifest}
#   \item{address_B}{address of probe B on the chip
#   designated by the original manifest}
#   \item{channel}{'Both' for type II probes and 'Grn' or 
#   'Red' for Type I probes}
#   \item{designType}{Type of probe, either 'I' or 'II'}
#   \item{nextBase}{the actual extension base (on the probe
#   strand) after bisulfite conversion ('A' or 'C' or 'T').
#   Unmapped probes have extension base labeled in the original
#   manifest.}
#   \item{nextBaseRef}{the extension base (on the hybridized /
#   template DNA) before bisulfite conversion ('A', 'C', 'G',
#   or 'T'). Unmapped probes have 'NA'.}
#   \item{probeType}{either 'cg', 'ch', or 'rs'}
#   \item{orientation}{either 'up or 'down', specifying 
#   whether the probe is positioned upstream (in smaller
#   co-ordinates) or downstream (in greater co-ordinates)
#   from the target}
#   \item{probeCpGcnt}{the number of additional CpGs in the
#   probe (not counting the interrogated CpG)}
#   \item{context35}{the number of CpGs in the [-35bp, +35bp]
#   window}
#   \item{probeBeg}{the mapped start position of the probe,
#   which is always 50bp long}
#   \item{probeEnd}{the mapped end position of the probe,
#   which is always 50bp long}
#   \item{ProbeSeq_A}{the probe sequence for allele A}
#   \item{ProbeSeq_B}{the probe sequence for allele B}
#   \item{gene}{comma separated list of gene annotations
#   (unique and alphabetically sorted). Gene models follow
#   GENCODE version 22 (hg38)}
#   \item{gene_HGNC}{comma separated list of gene annotations
#   (unique and alphabetically sorted). Genes are checked
#   using HGNChelper for compatibility with HGNC. Gene models
#   follows GENCODE version 22 (hg38)}
#   \item{chrm_A}{the mapping info for probe A excluding 
#   decoy chromosomes}
#   \item{beg_A}{the mapping info for probe A excluding 
#   decoy chromosomes}
#   \item{flag_A}{the mapping info for probe A excluding 
#   decoy chromosomes}
#   \item{mapQ_A}{the mapping quality score for probe A excluding 
#   decoy chromosomes, with 60 being the best}
#   \item{cigar_A}{the mapping info for probe A excluding 
#   decoy chromosomes}
#   \item{NM_A}{the mapping info for probe A excluding 
#   decoy chromosomes}
#   \item{chrm_B}{the mapping info for probe B excluding 
#   decoy chromosomes}
#   \item{beg_B}{the mapping info for probe B excluding 
#   decoy chromosomes}
#   \item{flag_B}{the mapping info for probe B excluding 
#   decoy chromosomes}
#   \item{mapQ_B}{the mapping quality score for probe B excluding 
#   decoy chromosomes, with 60 being the best}
#   \item{cigar_B}{the mapping info for probe B excluding 
#   decoy chromosomes}
#   \item{NM_B}{the mapping info for probe B excluding 
#   decoy chromosomes}
#   \item{wDecoy_chrm_A}{the mapping info for probe A including 
#   decoy chromosomes}
#   \item{wDecoy_beg_A}{the mapping info for probe A including 
#   decoy chromosomes}
#   \item{wDecoy_flag_A}{the mapping info for probe A including 
#   decoy chromosomes}
#   \item{wDecoy_mapQ_A}{the mapping quality score for probe A including 
#   decoy chromosomes, with 60 being the best}
#   \item{wDecoy_cigar_A}{the mapping info for probe A including 
#   decoy chromosomes}
#   \item{wDecoy_NM_A}{the mapping info for probe A including 
#   decoy chromosomes}
#   \item{wDecoy_chrm_B}{the mapping info for probe B including 
#   decoy chromosomes}
#   \item{wDecoy_beg_B}{the mapping info for probe B including 
#   decoy chromosomes}
#   \item{wDecoy_flag_B}{the mapping info for probe B including 
#   decoy chromosomes}
#   \item{wDecoy_mapQ_B}{the mapping quality score for probe B including 
#   decoy chromosomes, with 60 being the best}
#   \item{wDecoy_cigar_B}{the mapping info for probe B including 
#   decoy chromosomes}
#   \item{wDecoy_NM_B}{the mapping info for probe B including 
#   decoy chromosomes}
#   \item{posMatch}{whether the mapping matches the original
#   manifest, it only applies to hg19 and will be NA under 
#   hg38}
#   \item{MASK_mapping}{whether the probe is masked for mapping
#   reasons. Probes retained should have high quality (>=40)
#   consistent (with designed MAPINFO) mapping (for both in the
#   cast of type I) without INDELS}
#   \item{MASK_typeINextBaseSwitch}{whether the probe has a SNP
#   in the extension base that causes a color channel switch from
#   the official annotation (described as color channel switching
#   or CCS SNP in the reference). These probes should be processed
#   differently than designed (by summin up both color channels
#   instead of just the annotated color channel)}
#   \item{MASK_rmsk15}{whetehr the 15bp 3'-subsequence of the
#   probe overlaps with repeat masker, this MASK is NOT 
#   recommended}
#   \item{MASK_sub40_copy}{whether the 40bp 3'subsequence of the
#   probe is non-unique}
#   \item{MASK_sub35_copy}{whether the 35bp 3'subsequence of the
#   probe is non-unique}
#   \item{MASK_sub30_copy}{whether the 30bp 3'subsequence of the
#   probe is non-unique}
#   \item{MASK_sub25_copy}{whether the 25bp 3'subsequence of the
#   probe is non-unique}
#   \item{MASK_snp5_common}{whether 5bp 3'-subsequence (including
#   extension for tyep II) overlaps with any of the common SNPs
#   from dbSNP (global MAF can be under 1%)}
#   \item{MASK_snp5_GMAF1p}{whether 5bp 3'-subsequence (including
#   extension for type II) overlaps with any of the SNPs with 
#   global MAF > 1%}
#   \item{MASK_extBase}{probes masked for extension base inconsistent
#   with specified color channel (type I) or CpG (Type II) based
#   on mapping}
#   \item{MASK_general}{recommended general ppurpose masking
#   merged from MASK_sub30_copy, MASK_mapping, MASK_extBase,
#   MASK_typeINExtBaseSwitch and MASK_snp5_GMAF1p}
#
# @usage read_tsv('data/EPIC.hg19.manifest.tsv.gz')
#
# @source 
#         \url{http://zwdzwd.github.io/InfiniumAnnotation}
#         
NULL


# Data frame with coefficients for estimating blood counts
#
# Contains coefficients derived from the WBCC predictor
#
# @format A data frame with 481392 rows and 5 variables:
# \describe{
#   \item{baso_perc}{coefficient for estimating basophils}
#   \item{eos_perc}{coefficient for estimating eosinophils}
#   \item{lymph_perc}{coefficient for estimating lymphocytes}
#   \item{mono_perc}{coefficient for estimating monocytes}
#   \item{neut_perc}{coefficient for estimating neutrophils}
#
# @usage data(DNAmPredictorCoef)
#         
NULL