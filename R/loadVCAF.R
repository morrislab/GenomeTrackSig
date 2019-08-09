# loadVCAF.R
# author: Cait Harrigan
# R functions to depricate make_corrected_vaf.py and call_scripts.R

#' \code{vcfToCounts} Take an input vcf file and annotation and generate the counts data
#'
#' @rdname callScripts
#' @name vcfToCounts
#'
#' @param vcfFile path to variant calling format (vcf) file
#' @param cnaFile path to copy number abberation (cna) file
#' @param purityFile path to sample purity file
#'
#' @export
vcfToCounts <- function(vcfFile, cnaFile = NULL, purity = 1, binSize = 100,
                        context = generateContext(c("CG", "TA")),
                        refGenome = BSgenome.Hsapiens.UCSC.hg19,
                        saveIntermediate = F, intermediateFile = NULL) {

  # load CNA and purity dataframe (not loaded with VCF for parallelization memory saving)
  # could be done as a single annotation load.... one function to load each file
  # loads the following - all shared between all VCF's, all optional (but not necessarily independent)
  # cna, purity, tumortypes, signatures (alex, cosmic), trinucleotide, sigactivities
  # tumortype_file = "", signature_file = "", trinucleotide_file = "", active_signatures_file = ""

  # input checking and path expansion
  vcfFile <- path.expand(vcfFile)
  stopifnot(file.exists(vcfFile))

  if (!is.null(cnaFile)){
    cnaFile <- path.expand(cnaFile)
    stopifnot(file.exists(cnaFile))
  }

  if (!is.null(purityFile)){
    purityFile <- path.expand(purityFile)
    stopifnot(file.exists(purityFile))
  }

  # get vcf
  vcf <- parseVcfFile(vcfFile)

  # get cna reconstruction
  cna <- parseCnaFile(cnaFile)

  # annotate the vcf with copy number
  vcf <- annotateCn(vcf, cna)

  # vcaf has vcf and vaf data concatenated
  vcaf <- getVcaf(vcf, purity, refGenome)
  vcaf <- getTrinuc(vcaf, refGenome, saveIntermediate, intermediateFile)

  return( getBinCounts(vcaf, binSize, context) )


}


#' @rdname loadVCAF
#' @name parseVcf

parseVcfFile <- function(vcfFile, cutoff = 10000, refGenome = BSgenome.Hsapiens.UCSC.hg19){

  vcf <- VariantAnnotation::readVcf(vcfFile, genome = providerVersion(refGenome))

  # TODO: remove any duplicates

  # TODO: remove samples with missing ref or alt counts


  # implement cutoff if too many variants present in sample
  if (dim(vcf)[1] > cutoff){

    vcf <- sample(vcf, cutoff)

  }



  return(vcf)
}

#' @rdname loadVCAF
#' @name parseCna

parseCnaFile <- function(cnaFile){

  cnaGR <- read.table(cnaFile, header = T)
  cnaGR <- GRanges(cnaGR$chromosome, IRanges(cnaGR$start, cnaGR$end), cn = cnaGR$total_cn)

  return(cnaGR)
}

#' @rdname loadVCAF
#' @name annotateCn

parsePurityFile <- function(purityFile){

  purities <- read.table(purityFile, header = T)

  return(purities)
}

#' @rdname loadVCAF
#' @name annotateCn

annotateCn <- function(vcf, cnaGR = NULL){

  vcfGR <- rowRanges(vcf)

  # input type checking, cnaGR colname checking
  assertthat::assert_that(class(vcfGR) == "GRanges")
  if( !is.null(cnaGR) ){
    assertthat::assert_that(class(cnaGR) == "GRanges")
    if (is.null(cnaGR$cn)){ stop("cnaGR$cn not found. Please assure that cnaGR$cn is a valid indexing of cnaGR") }
  }

  # set all cn to 2
  vcfGR$cn <- 2

  # if cn reconstruction available, add cna metadata to vcf
  if (!is.null(cnaGR)){

    # look up cna for vcf regions
    overlaps <- GenomicRanges::findOverlaps(vcfGR, cnaGR)

    vcfGR$cn[to(overlaps)] <- cnaGR$cn[from(overlaps)]
  }

  rowRanges(vcf) <- vcfGR

  return(vcf)

}

#' \code{getVcaf} Take an input vcf file and annotation and make vaf data
#'
#' @rdname callScripts
#' @name getVcaf
#'
#' @param vcf CollapsedVCF object with cna annotation
#' @param purityFile path to sample purity file
#' @param refGenome reference BSgenome to use
#' @return A vcaf dataframe that has vcf and vaf data concatenated

getVcaf <- function(vcfVA, purity, refGenome){
  #replaces make_corrected_vaf.py

  print("making vaf")

  # formatting - vcf and vaf concatenated and dataframe hold strings
  vcaf <- as.data.frame(rowRanges(vcfVA))[c("seqnames", "start", "REF", "ALT", "cn")]

  # want alt and ref counts for all loci
  assertthat::assert_that(dim(vcaf)[1] == length(info(vcfVA)$t_alt_count), dim(vcaf)[1] == length(info(vcfVA)$t_ref_count))



  colnames(vcaf) <- c("chr", "pos", "ref", "alt", "phi", "phi2", "vi", "ri", "purity")
  vcaf$phi <- as.numeric(vcaf$phi)
  vcaf$phi2 <- as.numeric(vcaf$phi2)
  vcaf$vi <- as.numeric(vcaf$vi)
  vcaf$ri <- as.numeric(vcaf$ri)
  vcaf$pos <- as.numeric(vcaf$pos)
  vcaf$purity <- as.numeric(vcaf$purity)

  # multiallelic hits keep only the first allele
  vcaf$alt <- substr(vcaf$alt, 2, 2)

  # order mutations by phi
  vcaf <- vcaf[order(vcaf$phi, decreasing = T),]

  # prelim formatting check
  vcaf <- checkVcaf(vcaf, refGenome)

  return(vcaf)
}

#' \code{checkVcaf} Perform some shallow input checks on a vcaf data frame. \cr
#' Check for SNP criteria, and remove instances where reference allele matches alt allele.\cr
#' Check chromosome and position is valid in reference genome.
#'
#' @rdname callScripts
#' @name checkVcaf
#'
#' @param vcaf vcaf data frame
#' @param refGenome reference BSgenome to use
#' @return A vcaf dataframe that has vcf and vaf data concatenated

checkVcaf <- function(vcaf, refGenome = BSgenome.Hsapiens.UCSC.hg19){

  # input checking
  assertthat::assert_that(class(refGenome) == "BSgenome")
  assertthat::assert_that(class(vcaf) == "data.frame")
  assertthat::assert_that(all(colnames(vcaf) == c("chr", "pos", "ref", "alt", "phi", "phi2", "vi", "ri", "purity")))

  # some VCF formatting checks, filter for SNP's
  # no read quality filtering performed.

  # ref should not match alt in a mutation
  rmSet <- vcaf$ref == vcaf$alt
  if (sum(rmSet) > 0){
    warning(sprintf("%s mutations dropped for refrence allele matching alt", length(rmSet)))
    vcaf <- vcaf[!rmSet,]
  }

  # mutation should be a SNP
  rmSet <- c()

  # lists of >2 alt alleles not SNP
  # don't count incluce python list characters - [,]
  rmSet <- union(rmSet, which(nchar(vcaf$alt) > 5))

  # lists of >1 ref allele not SNP
  # don't count incluce python list characters - [,]
  rmSet <- union(rmSet, which(nchar(vcaf$ref) > 3))

  # chromosome should be valid in refrence genome
  # don't load genome - use BSgenome previewe accessing
  # "chr" is stripped by make_vcaf so return for matching with BSgenome
  rmSet <- union(rmSet, which(!(paste0("chr", vcaf$chr) %in% seqnames(refGenome))))

  # postition should be valid in refrence genome
  # not less than 1
  rmSet <- union(rmSet, which( vcaf$pos < 1 ) )

  #and less than the maximum for that chromosome
  rmSet <- union(rmSet, which ( ! ( vcaf$pos < seqlengths(refGenome)[paste0("chr", vcaf$chr)] ) ))

  if (length(rmSet) > 0){
    warning( sprintf("%s mutations dropped for not meeting SNP cirteria" , length(rmSet) ) )
    vcaf <- vcaf[-rmSet,]
  }

  return ( vcaf )
}

#' @rdname loadVCAF
#' @name annotateCn
#' \code{annotateCn} add copy number aberation information as metadata to the data read from VCF file
#' @param vcfGR GRanges object with vcf data
#' @param cnaGR GRanges object with copy number as metadata, stored in cnaGR$cn. If none provided, copy number is assumed to be two for all loci.
#' @return passed vcfGR, augmented with copy number metadata



# [END]
