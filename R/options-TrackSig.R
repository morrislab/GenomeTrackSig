# options-TrackSig.R
# Author: Cait Harrigan
# Replicated from CRAN package MNF nmf.options()

#' \code{options-TrackSig} TrackSig Package specific options
#'
#' @section Available options:
#' \describe{
#'
#' \item{compute_bootstrap}{Whether to perform bootstrapping for uncertainty calculation}
#'
#' } % end description
#'
#' @param ... option specifications. For TrackSig.options this can be named\cr
#' arguments or a single named list (see \link{options})\cr
#' For TrackSig.resetOptions, this must be the names of the options to reset. \cr
#' Note that pkgmaker version >= 0.9.1 is required for this to work correctly,\cr
#' when options other than the default have been set after the package is loaded.
#'
#' @rdname options
#' @name options-TrackSig

NULL
.allOPTS <- pkgmaker::setupPackageOptions(
  #defaults from header.R

  # Select fitting all signatures or only signatures for the particular cancer type
  #sig_amount <- "onlyKnownSignatures" # recommended
  sig_amount = "onlyKnownSignatures" # not recommended, time-consuming

  # if the signatures are specified per cancer type or per sample
  , cancer_type_signatures = TRUE

  # if signatures trajectories need to be computed on bootstrapped signatures as well
  # bootstrapping provides the uncertainty estimations on the trajectories
  # warning: by default, mutations are bootstrapped 30 times and the script will run 30 time longer
  , compute_bootstrap = TRUE
  , sliding_window = FALSE
  , noise_sig = NULL
  , simulated_data = FALSE
  , postfix = ""

  # specifies the changepoint detection algorithm.
  , changepoint_method = "PELT"

  # file with cancer types of each sample
  , tumortype_file = "data/tumortypes.txt"

  # folders with mutation counts, mutation order and bootstrapped mutations
  # don't need to be changed unless different folder were specified in make_counts.sh
  , DIR_COUNTS = "data/counts/"
  , mutation_order = "data/mut_order/"
  , BOOTSTRAP_COUNTS = "data/bootstrap/"
  , purity_file = "data/example_purity.txt"

  # folder to write results to
  , DIR_RESULTS = "results_signature_trajectories/"

  # file with signatures definitions
  , signature_file = "annotation/alexSignatures.txt"

  # file with trinucleotide context
  , trinucleotide_file = "annotation/trinucleotide.txt"

  # specifies active signatures in TCGA cancer types
  , active_signatures_file = "annotation/active_signatures_transposed.txt"

  # specifies active signatures in each sample. Contains the active signatures for the example
  # active_signatures_file = "annotation/active_in_samples.txt"

  , SAVED_SAMPLES_DIR = "saved_data/"

  , PLOT_FULL_NAME = FALSE
  , mutation_assignments = ""

)

#check DIR_RESULTS exists
ResultDir_check <- function(){
  if(!file.exists(TrackSig.options()$DIR_RESULTS)){
    stop(sprintf("DIR_RESULTS file, %s doesn't exist.
                 See ?TrackSig.options for more information"),
         TrackSig.options()$DIR_RESULTS)
  }
}


#' \code{TrackSig.options} sets/get single or multiple options, that are \cr
#' specific to the TrackSig package. It behaves in the same way as \link{options}.
#'
#' @inheritParams base::options
#' @rdname options
#' @export

TrackSig.options <- .allOPTS$options

#' \code{TrackSig.getOption} returns the value of a single option, that is \cr
#' specific to the TrackSig package. It behaves in the same way as \cr
#' \link{getOption}.
#'
#' @inheritParams base::getOption
#' @rdname options
#' @export
TrackSig.getOption <- .allOPTS$getOption

#' \code{TrackSig.resetOptions} reset all TrackSig specific options to their \cr
#' default values.
#'
#' @param ALL logical that indicates if options that are not part of the default set
#' of options should be removed.
#' Note that in \pkg{pkgmaker <= 0.9} this argument is only taken into account when
#' no other argument is present. This is fixed in version 0.9.1.
#'
#' @rdname options
#' @export
TrackSig.resetOptions <- .allOPTS$resetOptions

#' \code{TrackSig.printOptions} prints all TrackSig specific options along with \cr
#' their default values, in a relatively compact way.
#' @rdname options
#' @export
TrackSig.printOptions <- .allOPTS$printOptions
# [END]
