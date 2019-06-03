# call_scripts.R

#' \code{call_scripts} Call the supporting (non-R) scipts in TrackSig
#'
#' \code{vcf_to_counts} From a supplied VCF file run the necessary scripts to get corrected_vaf and make_counts output. Return a multi-slotted object with useful dataframes.
#'
#' \code{run_simulation} depricated
#'


#' @rdname call_scripts
#' @name vcf_to_counts
#'
#' @export

vcf_to_counts <- function(vcf_file, cna_file = NULL, purity_file = NULL){



}


#' @rdname call_scripts
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
