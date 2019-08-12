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

  # vcaf has vcf and vaf data concatenated
  vcaf <- getVcaf(vcf, purity, cna, refGenome)
  vcaf <- getTrinuc(vcaf, refGenome, saveIntermediate, intermediateFile)

  return( getBinCounts(vcaf, binSize, context) )


}


#' @rdname loadVCAF
#' @name parseVcfFile

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
#' @name parseCnaFile

parseCnaFile <- function(cnaFile){

  cnaGR <- read.table(cnaFile, header = T)
  cnaGR <- GRanges(cnaGR$chromosome, IRanges(cnaGR$start, cnaGR$end), cn = cnaGR$total_cn)

  return(cnaGR)
}

#' @rdname loadVCAF
#' @name parsePurityFile

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

  # update header information with cn
  cnInfoHeader <- data.frame(row.names = c("cn"), Number = 1, Type = "Integer",
                             Description = "Locus copy number as provided to TrackSig",
                             stringsAsFactors = F)
  info(header(vcf)) <- rbind(info(header(vcf)), cnInfoHeader)

  info(vcf)$cn <- vcfGR$cn

  return(vcf)

}

#' \code{getVcaf} Take an input vcf file and annotation and make vaf data
#'
#' @rdname callScripts
#' @name getVcaf
#'
#' @param vcf CollapsedVCF object
#' @param purity sample purity percentage between 0 and 1
#' @param cna GRanges object with cna information for the sample
#' @param refGenome reference BSgenome to use
#' @return A vcaf dataframe that has vcf and vaf data concatenated

getVcaf <- function(vcf, purity, cna, refGenome = BSgenome.Hsapiens.UCSC.hg19){
  #replaces make_corrected_vaf.py

  # annotate the vcf with copy number
  vcf <- annotateCn(vcf, cna)

  # prelim formatting check
  vcaf <- vcafConstruction(vcf, refGenome)

  vcaf$vi <- info(vcf)$t_alt_count
  vcaf$ri <- info(vcf)$t_ref_count
  vcaf$purity <- purity
  vcaf$phat1 <- rbeta(dim(vcaf)[1], vcaf$vi + 1, vcaf$ri + 1)
  vcaf$phat2 <- rbeta(dim(vcaf)[1], vcaf$vi + 1, vcaf$ri + 1)
  vcaf$phi <- (2 + vcaf$purity * (vcaf$cn - 2)) * vcaf$phat1   #phi = ccf * purity
  vcaf$phi2 <- (2 + vcaf$purity * (vcaf$cn - 2)) * vcaf$phat2

  # TODO: depricate phi
  # sort on phi
  vcaf <- vcaf[order(vcaf$phi, decreasing = T), ]

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

vcafConstuction <- function(vcf, refGenome = BSgenome.Hsapiens.UCSC.hg19){
  # some VCF formatting checks, filter for SNP's
  # no read quality filtering performed.

  # input checking
  assertthat::assert_that(class(refGenome) == "BSgenome")
  assertthat::assert_that(class(vcf) == "CollapsedVCF")
  assertthat::assert_that("REF" %in% colnames(fixed(vcf)))
  assertthat::assert_that("ALT" %in% colnames(fixed(vcf)))

  # mutations in vcf should be SNPs => one ref, one alt allele
  # drop those that are not SNVs
  rmSel <- !isSNV(vcf, singleAltOnly = F)

  if (sum(rmSel) > 0){
    warning( sprintf("%s mutations dropped for not meeting SNP cirteria" , sum(rmSel) ) )
    vcaf <- vcaf[!rmSel,]
  }

  # formatting vcaf - vcf and vaf concatenated
  # using CharacterList for alt, ref is 400x faster than casting CollapsedVCF object directly
  # subset with [0,] to avoid casting metadata columns
  vcaf <- as.data.frame(rowRanges(vcf)[,0])[c("seqnames", "start")]

  # take first variant of those with multiallelic alternate variant
  allAlts <- CharacterList(rowRanges(vcf)$ALT)
  vcaf$alt <- unlist(lapply(allAlts, function(x){return(x[[1]][1])}))

  # want alt and ref counts for all loci
  assertthat::assert_that(dim(vcaf)[1] == length(info(vcf)$t_alt_count), dim(vcaf)[1] == length(info(vcf)$t_ref_count))
  colnames(vcaf) <- c("chr", "pos", "ref", "alt", "cn")

  # TODO: vcf formats may provide vi and depth, rather than vi and ri

  # drop secondary alleles in multiallelic alt hits
  oneCharAllele <- function(alt){alt[[1]][1]}

  # ref should not match alt in a mutation
  rmSel <- vcaf$ref == vcaf$alt
  if (sum(rmSel) > 0){
    warning(sprintf("%s mutations dropped for refrence allele matching alt", length(rmSet)))
    vcaf <- vcaf[!rmSel,]
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


# [END]
