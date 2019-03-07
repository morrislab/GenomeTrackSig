# call_scripts.R

#' \code{call_scripts} Call the supporting (non-R) scipts in TrackSig
#'
#' @rdname call_scripts
#' @name run_simulation
#'
#' @export

run_simulation <- function(simName){

  fileName <- system.file("scripts", "run_simulations.sh", package = "TrackSig")
  packagePath <- system.file(package = "TrackSig")

  system(sprintf("%s %s data/mut_types/ data/%s %s %s", fileName, packagePath, simName, simName, TrackSig.options()$bin_size))

}

run_sample <- function(sampleName){

  #fileName <- system.file("scripts", "run_simulations.sh", package = "TrackSig")
  #packagePath <- system.file(package = "TrackSig")

  #system(sprintf("%s %s data/mut_types/ data/%s %s", fileName, packagePath, simName, simName))

  warning("run_sample not implemented")
  return(NULL)

}

# [END]
