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


  # from helper_functions.R
  COLORS2 = c("#000000","#1CE6FF","#FF34FF", "#FFFF00", "#FF4A46","#008941",
                "#006FA6","#A30059","#FFDBE5","#7A4900","#0000A6","#63FFAC",
                "#B79762","#004D43","#8FB0FF","#997D87","#5A0007","#809693",
                "#FEFFE6","#1B4400","#4FC601","#3B5DFF","#4A3B53","#FF2F80",
                "#61615A","#BA0900","#6B7900","#00C2A0","#FFAA92","#FF90C9",
                "#B903AA","#D16100","#DDEFFF","#000035","#7B4F4B","#A1C299",
                "#300018","#0AA6D8","#013349","#00846F","#372101","#FFB500",
                "#C2FFED","#A079BF","#CC0744","#C0B9B2","#C2FF99","#001E09",
                "#00489C","#6F0062","#0CBD66","#EEC3FF","#456D75","#B77B68",
                "#7A87A1","#788D66","#885578","#FAD09F","#FF8A9A","#D157A0",
                "#BEC459","#456648","#0086ED","#886F4C","#34362D","#B4A8BD",
                "#00A6AA","#452C2C","#636375","#A3C8C9","#FF913F","#938A81",
                "#575329","#00FECF","#B05B6F","#8CD0FF","#3B9700","#04F757",
                "#C8A1A1","#1E6E00","#7900D7","#A77500","#6367A9","#A05837",
                "#6B002C","#772600","#D790FF","#9B9700","#549E79","#FFF69F",
                "#201625","#72418F","#BC23FF","#99ADC0","#3A2465","#922329",
                "#5B4534","#FDE8DC","#404E55","#0089A3","#CB7E98","#A4E804",
                "#324E72","#6A3A4C","#83AB58","#001C1E","#D1F7CE","#004B28",
                "#C8D0F6","#A3A489","#806C66","#222800","#BF5650","#E83000",
                "#66796D","#DA007C","#FF1A59","#8ADBB4","#1E0200","#5B4E51",
                "#C895C5","#320033","#FF6832","#66E1D3","#CFCDAC","#D0AC94",
                "#7ED379","#012C58","#7A7BFF","#D68E01","#353339","#78AFA1",
                "#FEB2C6","#75797C","#837393","#943A4D","#B5F4FF","#D2DCD5",
                "#9556BD","#6A714A","#001325","#02525F","#0AA3F7","#E98176",
                "#DBD5DD","#5EBCD1","#3D4F44","#7E6405","#02684E","#962B75",
                "#8D8546","#9695C5","#E773CE","#D86A78","#3E89BE","#CA834E",
                "#518A87","#5B113C","#55813B","#E704C4","#00005F","#A97399",
                "#4B8160","#59738A","#FF5DA7","#F7C9BF","#643127","#513A01",
                "#6B94AA","#51A058","#A45B02","#1D1702","#E20027","#E7AB63",
                "#4C6001","#9C6966","#64547B","#97979E","#006A66","#391406",
                "#F4D749","#0045D2","#006C31","#DDB6D0","#7C6571","#9FB2A4",
                "#00D891","#15A08A","#BC65E9","#FFFFFE","#C6DC99","#203B3C",
                "#671190","#6B3A64","#F5E1FF","#FFA0F2","#CCAA35","#374527",
                "#8BB400","#797868","#C6005A","#3B000A","#C86240","#29607C",
                "#402334","#7D5A44","#CCB87C","#B88183","#AA5199","#B5D6C3",
                "#A38469","#9F94F0","#A74571","#B894A6","#71BB8C","#00B433",
                "#789EC9","#6D80BA","#953F00","#5EFF03","#E4FFFC","#1BE177",
                "#BCB1E5","#76912F","#003109","#0060CD","#D20096","#895563",
                "#29201D","#5B3213","#A76F42","#89412E","#1A3A2A","#494B5A",
                "#A88C85","#F4ABAA","#A3F3AB","#00C6C8","#EA8B66","#958A9F",
                "#BDC9D2","#9FA064","#BE4700","#658188","#83A485","#453C23",
                "#47675D","#3A3F00","#061203","#DFFB71","#868E7E","#98D058",
                "#6C8F7D","#D7BFC2","#3C3E6E","#D83D66","#2F5D9B","#6C5E46",
                "#D25B88","#5B656C","#00B57F","#545C46","#866097","#365D25",
                "#252F99","#00CCFF","#674E60","#FC009C","#92896B")

  # Select fitting all signatures or only signatures for the particular cancer type
  , sig_amount = "onlyKnownSignatures" # "onlyKnownSignatures" recommended, "all" not recommended, time-consuming

  # if the signatures are specified per cancer type or per sample
  , cancer_type_signatures = TRUE

  # if signatures trajectories need to be computed on bootstrapped signatures as well
  # bootstrapping provides the uncertainty estimations on the trajectories
  # warning: by default, mutations are bootstrapped 30 times and the script will run 30 time longer
  , compute_bootstrap = TRUE
#  , sliding_window = FALSE
  , noise_sig = NULL
#  , simulated_data = FALSE
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

  # from compute_mutational_signatures.R
  , group = 0
  , EXAMPLES_PER_GROUP = 500

  # signature reference format (for loading annotation)
  , pcawg_format = TRUE

  # bin size for timepoints
  , bin_size = 100

  # likelihoods
  , pelt_penalty = expression((n_sigs - 1) * log(n_bins))
  , pelt_score_fxn = TrackSig:::log_likelihood_mixture_multinomials

  # save space on simulation output
  , simulation_pdf = TRUE

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
