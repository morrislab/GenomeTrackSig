# SpaceTrack.R
# Defines main functions for user to interact with package TrackSig, when calculating activity trajectory over genome.
# Author: Caitlin Timmons


#' Determine genomic trajectory
#'
#' @description
#' \code{SpaceTrack} will take an input CSV file of mutation counts, and
#' determine a trajectory over the genome for the sample, based on changepoints
#' found using the PELT segmentation algorithm.
#'
#' @param path path to mutation counts (csv) file.
#' @param chr_level logical whether TrackSig should be run separately on
#' each chromosome; default is FALSE
#' @param bootstrapSamples integer denoting number of bootstrap samples to take;
#' default is 0
#' @param parallelize logical whether functions should be run in parallel;
#' default is FALSE. Cannot be TRUE if chr_level is FALSE and bootstrapSamples is
#' 0.
#' @param activeInSample list of active signatures to fit
#' (must match signatures present in referenceSignatures)
#' @param binSize number of mutations per bin
#' @param referenceSignatures dataframe containing definitions of mutational
#'   signatures.
#' @param refGenome BSgenome to use as reference
#'
#' @return trajectory object containing one trajectory per bootstrap sample.
#' @export

SpaceTrack <- function(path, chr_level = FALSE, bootstrapSamples = 0,
                       parallelize = FALSE, activeInSample, binSize,
                       referenceSignatures = TrackSig::alex_merged,
                       refGenome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19) {

  # read in mutation counts file
  counts <- readr::read_csv(path)

  # run TrackSig either on entire genome or on individual chromosomes
  if (chr_level == TRUE) {
    traj <- trackChromosome(master = counts, activeInSample = activeInSample,
                            binSize = binSize, nSamples = bootstrapSamples,
                            parallelize = parallelize)
  }

  else {
    traj <- trackGenome(master = counts, activeInSample = activeInSample,
                        binSize = binSize, nSamples = bootstrapSamples,
                        parallelize = parallelize)
  }

  return (traj)
}
