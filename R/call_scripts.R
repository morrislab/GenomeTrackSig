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

vcfToCounts <- function(vcfFile, cnaFile = NULL, purityFile = NULL) {

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

  # call python with reticulate

  reticulate::source_python(system.file("python/make_corrected_vaf.py", package = "TrackSig"))
  a <- make_vaf(vcfFile, cnaFile, purityFile)


  # get mutation types


  # bin mutations

}


getMutTypes <- function(vcf, vaf){
  # replaces perl script

  vcf


}

getBins <- function(){
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
