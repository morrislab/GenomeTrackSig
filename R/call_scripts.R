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
                        context = trinucleotide_internal, refGenome = Hsapian
                        ) {

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

  vcaf <- getVcaf(vcfFile, cnaFile, purityFile)
  vcaf <- checkVcaf(vcaf)
  vcaf <- getMutTypes(vcaf)

  # bin mutations

}


getVcaf <- function(vcfFile, cnaFile, purityFile){
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

  return(vcaf)
}

checkVcaf <- function(vcaf, refGenome = Hsapiens){
  # some VCF formatting checks, filter for SNP's
  # no read quality filtering performed.

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
    warning( sprintf("%s loci were removed from the vcf data (they may have not met the SNP cirteria)" , length(rmSet) ) )
    vcaf <- vcaf[-rmSet,]
  }

  return ( vcaf )
}

getMutTypes <- function(vcaf, refGenome = Hsapiens){
  # replaces getMutationTypes.pl

  # get trinucleotide context in refrence
  # strandedness should be forward
  # concat to GRanges object
  mutRanges <- GRanges( paste0("chr", vcaf$chr, ":", vcaf$pos - 1, "-", vcaf$pos + 1, ":+") )

  # look up trinucleotide context
  context <- getSeq(refGenome, mutRanges)
  vcaf$context <- as.character(context)

  # context matches ref?
  # perl script ignored this and grabbed trinuc context regardless.
  # Here will drop these rows and throw warning
  mismatchedRef <- which(!(vcaf$ref == substr(vcaf$context, 2, 2)))
  if (length(mismatchedRef) > 0){
    warning( sprintf("%s loci were removed from the vcf data for not matching the selected reference" , length(mismatchedRef) ) )
    vcaf <- vcaf[-mismatchedRef,]
    context <- context[-mismatchedRef]
  }

  # take reverse complement of ref purines for context format
  complementSel <- (vcaf$ref == "G" | vcaf$ref == "A")
  vcaf$context[complementSel] <- as.character(reverseComplement(context)[complementSel])

  # swap ref purines to pyrimidines
  vcaf$ref[vcaf$ref == "G"] <- "C"
  vcaf$ref[vcaf$ref == "A"] <- "T"

  return (vcaf)
}

countBins <- function(vcaf, binSize, ){
  # calls make_hundreds script

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
