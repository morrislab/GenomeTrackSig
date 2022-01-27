# GenomeTrackSig.R
# Defines main functions for user to interact with package GenomeTrackSig, when calculating activity trajectory over genome.
# Author: Caitlin Timmons


#' Determine genomic trajectory
#'
#' @description
#' \code{GenomeTrackSig} will take an input data frame of mutation counts, and
#' determine a profile over the genome for the sample, based on changepoints
#' found using the PELT segmentation algorithm.
#'
#' @param counts data frame containing mutation counts.
#' @param chr_level logical whether GenomeTrackSig should be run separately on
#' each chromosome; default is FALSE
#' @param bootstrapMethod string indicating desired bootstrapping method; default 'None'.
#' 'Mutation' samples with replacement from mutations in each bin. 'Shuffle' samples a random
#' order of chromosomes (not applicable when chr_level = TRUE).
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

GenomeTrackSig <- function(counts, chr_level = FALSE, bootstrapMethod = 'None',
                       bootstrapSamples = 0,
                       parallelize = FALSE, activeInSample, binSize,
                       referenceSignatures = GenomeTrackSig:::cosmicV3,
                       refGenome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19) {

  # run GenomeTrackSig either on entire genome or on individual chromosomes
  if (chr_level == TRUE) {
    traj <- trackChromosome(master = counts, activeInSample = activeInSample,
                            binSize = binSize,
                            bootstrapMethod = bootstrapMethod,
                            nSamples = bootstrapSamples,
                            parallelize = parallelize)
  }

  else {
    traj <- trackGenome(master = counts, activeInSample = activeInSample,
                        binSize = binSize, bootstrapMethod = bootstrapMethod,
                        nSamples = bootstrapSamples,
                        parallelize = parallelize)
  }

  if (!is.null(traj)) {
    traj <- fixIndexing(traj, activeInSample)
  }

  return (traj)
}


fixIndexing <- function(traj, activeInSample) {
  # corrects trajectory format when no changepoints are identified
  # plotting functions expect trajectory components to be at certain indices
  changepoints <- list(changepoints = NULL)
  j <- 1
  new_traj <- list()
  for (i in 1:length(traj)) {
    if (typeof(traj[[i]])=='double') {
      if (!is.null(dim(traj[[i]])[1])) {
        if (as.numeric(dim(traj[[i]])[1])==length(activeInSample)) {
          if (typeof(traj[[i+1]])=='character') {
            new_traj <- c(new_traj, traj[i], changepoints)
            j <- j + 3
          }
          else {
            new_traj <- c(new_traj, traj[i])
          }
        }
        else {
          new_traj <- c(new_traj, traj[i])
        }
      }
      else {
        new_traj <- c(new_traj, traj[i])
      }
    }
    else {
      new_traj <- c(new_traj, traj[i])
    }
  }
  return (new_traj)
}
