# call_scripts.R

#' \code{callScripts} Call the supporting (non-R) scipts in TrackSig
#'
#' \code{vcf_to_counts} From a supplied VCF file run the necessary scripts to get corrected_vaf and make_counts output. Return a multi-slotted object with useful dataframes.
#'
#' \code{run_simulation} depricated
#'


#' @rdname callScripts
#' @name vcfToCounts
#'
#' @param vcfFile path to variant calling format (vcf) file
#' @param cnaFile path to copy number abberation (cna) file
#' @param purityFile path to sample purity file
#'
#' @export

vcfToCounts <- function(vcfFile, cnaFile = NULL, purityFile = NULL,
                        context = trinucleotide_internal, refGenome = Hsapians, binSize = 100) {

  # load CNA and purity dataframe (not loaded with VCF for parallelization memory saving)
  # could be done as a single annotation load.... one function to load each file
  # loads the following - all shared between all VCF's, all optional (but not necessarily independent)
  # cna, purity, tumortypes, signatures (alex, cosmic), trinucleotide, sigactivities
  # tumortype_file = "", signature_file = "", trinucleotide_file = "", active_signatures_file = ""

  # input checking and path expansion

  stopifnot(file.exists(vcfFile))
  vcfFile <- path.expand(vcfFile)

  if (!is.null(cnaFile)){
    stopifnot(file.exists(cnaFile))
    cnaFile <- path.expand(cnaFile)
  }

  if (!is.null(purityFile)){
    stopifnot(file.exists(purityFile))
    purityFile <- path.expand(purityFile)
  }

  vcaf <- getVcaf(vcfFile, cnaFile, purityFile, refGenome)
  vcaf <- getMutTypes(vcaf, refGenome)

  return( getBinCounts(vcaf, binSize, context) )


}


getVcaf <- function(vcfFile, cnaFile, purityFile, refGenome){
  #replaces make_corrected_vaf.py

  # call python with reticulate
  reticulate::source_python(system.file("python/make_corrected_vaf.py", package = "TrackSig"))

  # formatting - vcf and vaf concatenated and dataframe hold strings
  vcaf <- make_vcaf(vcfFile, cnaFile, purityFile)
  colnames(vcaf) <- c("chr", "pos", "ref", "alt", "phi")
  vcaf$phi <- as.numeric(vcaf$phi)
  vcaf$pos <- as.numeric(vcaf$pos)

  # multiallelic hits keep only the first allele
  vcaf$alt <- substr(vcaf$alt, 2, 2)

  # order mutations by phi
  vcaf <- vcaf[order(vcaf$phi, decreasing = T),]

  # prelim formatting check
  vcaf <- checkVcaf(vcaf, refGenome)

  return(vcaf)
}

checkVcaf <- function(vcaf, refGenome){

  # input checking
  assertthat::assert_that(class(refGenome) == "BSgenome")
  assertthat::assert_that(class(vcaf) == "data.frame")
  assertthat::assert_that(all(colnames(vcaf) == c("chr", "pos", "ref", "alt", "phi")))

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

#' @rdname callScripts
#' @name getMutTypes
#'
#' @param vcaf vcaf data frame
#' @param refGenome reference BSgenome to use
#' @param saveIntermediate boolean whether to save intermediate results (mutation types)
#' @param intermediateFile file where to save intermediate results if saveIntermediate is True
#' @return An updated vcaf data frame with trinucleotide context added for each mutation

getMutTypes <- function(vcaf, refGenome, saveIntermediate = F, intermediateFile){
  # replaces getMutationTypes.pl

  # input checking
  assertthat::assert_that(class(refGenome) == "BSgenome")
  assertthat::assert_that(is.logical(saveIntermediate))

  if(missing(intermediateFile)){
    assertthat::assert_that(saveIntermediate == F, msg = "please specify an intermediate file to save to, or set saveIntermediate = FALSE")
  }
  else{
    assertthat::assert_that(file.exists(intermediateFile))
  }


  # get trinucleotide context in refrence
  # strandedness should be forward
  # concat to GRanges object
  mutRanges <- GRanges( paste0("chr", vcaf$chr, ":", vcaf$pos - 1, "-", vcaf$pos + 1, ":+") )

  # look up trinucleotide context
  context <- getSeq(refGenome, mutRanges)
  vcaf$mutType <- as.character(context)

  # context matches ref?
  # perl script ignored this and grabbed trinuc context regardless.
  # Here will drop these rows and throw warning
  mismatchedRef <- which(!(vcaf$ref == substr(vcaf$mutType, 2, 2)))
  if (length(mismatchedRef) > 0){
    warning( sprintf("%s mutations dropped for vcf refrence allele not matching the selected reference genome" , length(mismatchedRef) ) )
    vcaf <- vcaf[-mismatchedRef,]
    context <- context[-mismatchedRef]
  }

  # take reverse complement of ref purines for context format
  complementSel <- (vcaf$ref == "G" | vcaf$ref == "A")
  vcaf$mutType[complementSel] <- as.character(reverseComplement(context)[complementSel])

  # complement alt and ref where ref is a purine
  vcaf$alt[vcaf$ref == "G"] <- as.character(complement(DNAStringSet(vcaf$alt[vcaf$ref == "G"])))
  vcaf$alt[vcaf$ref == "A"] <- as.character(complement(DNAStringSet(vcaf$alt[vcaf$ref == "A"])))
  vcaf$ref[vcaf$ref == "G"] <- "C"
  vcaf$ref[vcaf$ref == "A"] <- "T"

  if (saveIntermediate == TRUE){
    write.table(vcaf, file = intermediateFile)
  }

  return (vcaf)
}


#' @rdname callScripts
#' @name getBinCounts
#'
#' @param vcaf vcaf data frame
#' @param binSize number of mutations per bin
#' @param context trinucleotide combinations possible
#' @return A data frame of summary statistics and mutation type counts for each bin.
getBinCounts <- function(vcaf, binSize, context){
  # replaces make_hundreds.py script

  nMut <- dim(vcaf)[1]
  assertthat::assert_that(nMut > binSize, msg = "number of mutations may not be less than specified bin size")
  assertthat::assert_that(dim(unique(vcaf[c("ref","alt","mutType")]))[1] <= dim(context)[1], msg = "too many mutation types for context")

  #nBins <- (nMut / binSize) + (nMut %/% binSize > 0)
  nBins <- floor( (nMut / binSize) )

  # only filling as many complete bins as we can
  # up to the last (binSize - 1) mutations of smallest phi may be excluded.
  vcaf$binAssignment <- c(rep(1:nBins, each = binSize), rep(NA, nMut - nBins * binSize))

  # aggregate on bins
  binCounts <- data.frame(row.names = 1:nBins)

  # mean phis
  binCounts["phi"] <- aggregate(vcaf$phi, by = list(vcaf$binAssignment), FUN = mean)$x

  # meansq phis
  binCounts["quadPhi"] <- aggregate(vcaf$phi, by = list(vcaf$binAssignment), FUN = function(x){return(mean(x ^ 2))})$x

  # counts for each bin
  binCounts <- cbind( binCounts, aggregate(paste(vcaf$ref, vcaf$alt, vcaf$mutType, sep = "_"), by = list(vcaf$binAssignment), FUN = function(x){return(as.array(table(x)))})$x )

  # check that all mutation types have a count
  missingTypes <- setdiff(paste(context$V1, context$V2, context$V3, sep = "_"), names(binCounts))
  for (col in missingTypes){
    binCounts[col] <- 0
  }

  return (binCounts)

}





#' @rdname callScripts
#' @name run_simulation
#'
#' @export

run_simulation <- function(simName, dataDir){

  fileName <- system.file("scripts", "run_simulations.sh", package = "TrackSig")
  packagePath <- system.file(package = "TrackSig")

  system(sprintf("%s %s %s %s/%s %s %s", fileName, packagePath, dataDir, dataDir, simName, simName, TrackSig.options()$bin_size))

}

run_sample <- function(sampleName){

  #fileName <- system.file("scripts", "run_simulations.sh", package = "TrackSig")
  #packagePath <- system.file(package = "TrackSig")

  #system(sprintf("%s %s data/mut_types/ data/%s %s", fileName, packagePath, simName, simName))

  warning("run_sample not implemented")
  return(NULL)

}




# [END]
