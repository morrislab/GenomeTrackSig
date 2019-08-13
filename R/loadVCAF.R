# loadVCAF.R
# author: Cait Harrigan
# R functions to depricate make_corrected_vaf.py and call_scripts.R

#' \code{vcfToCounts} Take an input vcf file and annotation and generate the counts data
#'
#' @rdname callScripts
#' @name vcfToCounts
#'
#' @param vcfFile path to variant calling format (VCF) file
#' @param cnaFile path to copy number abberation (CNA) file
#' @param purity sample purity percentage between 0 and 1
#' @param binSize size of TrackSig timeline bins
#' @param context list of mutation types and their trinucleotide context
#' @param refGenome refrence genome used to create VCF file
#'
#' @export
vcfToCounts <- function(vcfFile, cnaFile = NULL, purity = 1, binSize = 100,
                        context = generateContext(c("CG", "TA")),
                        refGenome = BSgenome.Hsapiens.UCSC.hg19) {


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
  cnaRanges <- parseCnaFile(cnaFile)

  # vcaf has vcf and vaf data concatenated
  vcaf <- getVcaf(vcf, purity, cnaRanges, refGenome)
  vcaf <- getTrinuc(vcaf, refGenome)

  return( getBinCounts(vcaf, binSize, context) )


}


#' @export
vcfToCounts_simulation <- function(vcfFile, mutTypesFile, cnaFile = NULL, purityFile = NULL,
                                   context = generateContext(c("CG", "TA")), refGenome = BSgenome.Hsapiens.UCSC.hg19, binSize = 100,
                                   saveIntermediate = F, intermediateFile = NULL) {

  # load CNA and purity dataframe (not loaded with VCF for parallelization memory saving)
  # could be done as a single annotation load.... one function to load each file
  # loads the following - all shared between all VCF's, all optional (but not necessarily independent)
  # cna, purity, tumortypes, signatures (alex, cosmic), trinucleotide, sigactivities
  # tumortype_file = "", signature_file = "", trinucleotide_file = "", active_signatures_file = ""

  # input checking and path expansion


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
  cnaRanges <- parseCnaFile(cnaFile)

  # vcaf has vcf and vaf data concatenated
  vcaf <- getVcaf(vcf, purity, cnaRanges, refGenome)
  mutTypes <- read.delim(mutTypesFile, stringsAsFactors = F)

  # strip chr if present
  mutTypes$chromosome <- unlist(strsplit(mutTypes$chromosome, "chr"))[c(F, T)]

  # match ordering
  vcaf <- vcaf[order(vcaf$pos, vcaf$chr),]
  mutTypes <- mutTypes[order(mutTypes$start, mutTypes$chromosome),]

  # sanity check
  stopifnot( all(vcaf$chr == mutTypes$chromosome) & all(vcaf$pos == mutTypes$start) )

  # set trinucs and re-order by phi
  vcaf$mutType <- mutTypes$tri
  vcaf <- vcaf[order(vcaf$phi),]

  # save intermediate is necessary
  if (saveIntermediate == TRUE){
    mut_types <- vcaf[,c("chr", "pos", "phi", "ref", "alt", "mutType")]
    write.table(mut_types, file = intermediateFile, quote = F, col.names = F, row.names = F, sep = "\t")
  }

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

  # populate vcaf with phi calculations
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

vcafConstruction <- function(vcf, refGenome = BSgenome.Hsapiens.UCSC.hg19){
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
    vcf <- vcf[!rmSel,]
  }

  # formatting vcaf - vcf and vaf concatenated
  # using CharacterList for alt, ref is 400x faster than casting CollapsedVCF object directly
  # subset with [0,] to avoid casting metadata columns
  vcaf <- as.data.frame(rowRanges(vcf)[,0])[c("seqnames", "start")]
  colnames(vcaf) <- c("chr", "pos")

  vcaf$ref <- unlist(CharacterList(list(rowRanges(vcf)$REF)))

  # drop secondary alleles in multiallelic alt hits
  allAlts <- CharacterList(rowRanges(vcf)$ALT)
  vcaf$alt <- unlist(lapply(allAlts, function(x){return(x[[1]][1])}))

  # check - names and dimenions should match
  assertthat::assert_that( all( rownames(vcaf) == names(vcf) ) )
  assertthat::assert_that( dim(vcaf)[1] == dim(info(vcf))[1] )

  # get copy number for all loci
  vcaf$cn <- info(vcf)$cn

  # get alt and ref counts for all loci
  # TODO: vcf formats may provide vi and depth, rather than vi and ri
  vcaf$vi <- info(vcf)$t_alt_count
  vcaf$ri <- info(vcf)$t_ref_count

  # check - ref should not match alt in a mutation
  rmSel <- vcaf$ref == vcaf$alt
  if (sum(rmSel) > 0){
    warning(sprintf("%s mutations dropped for refrence allele matching alt", sum(rmSel)))
    vcaf <- vcaf[!rmSel,]
  }

  # check - mutations should be SNP
  rmSet <- c()

  # check - chromosome should be valid in refrence genome
  # don't load genome - use BSgenome preview accessing
  rmSet <- union(rmSet, which(!(vcaf$chr %in% seqnames(refGenome))))

  # postition should be valid in refrence genome
  # not less than 1
  rmSet <- union(rmSet, which( vcaf$pos < 1 ) )
  #and less than the maximum for that chromosome
  rmSet <- union(rmSet, which ( ! ( vcaf$pos < seqlengths(refGenome)[paste0("chr", vcaf$chr)] ) ))

  if (length(rmSet) > 0){
    warning( sprintf("%s mutations dropped for not appearing in reference genome" , length(rmSet) ) )
    vcaf <- vcaf[-rmSet,]
  }

  # strip "chr" for downstream
  vcaf$chr <- as.character(vcaf$chr)
  vcaf$chr <- unlist(strsplit(vcaf$chr, "chr"))[c(F, T)]

  return ( vcaf )
}

#' \code{getTrinuc} Get the trinucleotide context for each mutation in a vcaf data frame
#' @rdname callScripts
#' @name getTrinuc
#'
#' @param vcaf vcaf data frame
#' @param refGenome reference BSgenome to use
#' @param intermediateFile file where to save intermediate results if saveIntermediate is True
#' @return An updated vcaf data frame with trinucleotide context added for each mutation
getTrinuc <- function(vcaf, refGenome){
  # replaces getMutationTypes.pl

  print("making mutation types")

  # input checking
  assertthat::assert_that(class(refGenome) == "BSgenome")

  # get trinucleotide context in refrence
  # strandedness should be forward
  # concat to GRanges object
  mutRanges <- GRanges( paste0("chr", vcaf$chr, ":", vcaf$pos - 1, "-", vcaf$pos + 1, ":+") )

  # look up trinucleotide context
  triNuc <- getSeq(refGenome, mutRanges)
  vcaf$mutType <- as.character(triNuc)

  # context matches ref?
  # perl script ignored this and grabbed trinuc context regardless.
  # Here will drop these rows and throw warning
  mismatchedRef <- which(!(vcaf$ref == substr(vcaf$mutType, 2, 2)))

  if (length(mismatchedRef) > 0){

    #warning( sprintf("%s mutations dropped for vcf refrence allele not matching the selected reference genome" , length(mismatchedRef) ) )
    #vcaf <- vcaf[-mismatchedRef,]
    #context <- context[-mismatchedRef]

    warning( sprintf("%s mutations have vcf refrence allele mismatch with the selected reference genome" , length(mismatchedRef) ) )
    substr(vcaf$mutType, 2, 2) <- vcaf$ref
  }

  # remove mutations with "N" in refrence context
  rmSet <- sapply(triNuc, FUN = BSgenome::hasOnlyBaseLetters)
  if (sum(rmSet) > 0){

    warning( sprintf("%s mutations dropped for uncertain identity in reference genome" , sum(!rmSet)) )
    vcaf <- vcaf[rmSet,]
  }

  # take reverse complement of ref purines for context format
  complementSel <- (vcaf$ref == "G" | vcaf$ref == "A")
  vcaf$mutType[complementSel] <- as.character(reverseComplement(DNAStringSet(vcaf$mutType))[complementSel])

  # complement alt and ref where ref is a purine
  vcaf$alt[vcaf$ref == "G"] <- as.character(complement(DNAStringSet(vcaf$alt[vcaf$ref == "G"])))
  vcaf$alt[vcaf$ref == "A"] <- as.character(complement(DNAStringSet(vcaf$alt[vcaf$ref == "A"])))
  vcaf$ref[vcaf$ref == "G"] <- "C"
  vcaf$ref[vcaf$ref == "A"] <- "T"


  return (vcaf)
}

#' \code{getBinCounts} Get the mutation type counts data for a vcaf dataframe
#'
#' @rdname callScripts
#' @name getBinCounts
#'
#' @param vcaf vcaf data frame
#' @param binSize number of mutations per bin
#' @param context trinucleotide combinations possible
#' @return A data frame of summary statistics and mutation type counts for each bin.
getBinCounts <- function(vcaf, binSize, context){
  # replaces make_hundreds.py script

  print("making counts")

  nMut <- dim(vcaf)[1]
  assertthat::assert_that(nMut > binSize, msg = "number of mutations may not be less than specified bin size")
  assertthat::assert_that(dim(unique(vcaf[c("ref","alt","mutType")]))[1] <= dim(context)[1], msg = sprintf("too many mutation types (%s) for context (%s)",
                                                                                                           dim(unique(vcaf[c("ref","alt","mutType")]))[1],  dim(context)[1]) )

  #nBins <- (nMut / binSize) + (nMut %/% binSize > 0)
  nBins <- floor( (nMut / binSize) )

  # only filling as many complete bins as we can
  # up to the last (binSize - 1) mutations of smallest phi may be excluded.
  vcaf$binAssignment <- c(rep(1:nBins, each = binSize), rep(NA, nMut - nBins * binSize))

  # aggregate on bins
  binCounts <- data.frame(row.names = 1:nBins)

  # counts for each bin
  binCounts <- cbind (binCounts, aggregate(paste(vcaf$ref, vcaf$alt, vcaf$mutType, sep = "_"), by = list(vcaf$binAssignment), FUN = function(x){return(as.array(table(x)))})$x )

  # check that all mutation types have a count
  missingTypes <- setdiff(paste(context$V1, context$V2, context$V3, sep = "_"), names(binCounts))
  for (col in missingTypes){
    binCounts[col] <- 0
  }

  return ( list(vcaf, t(binCounts)) )

}

# [END]
