## \code{trackParallel} Run TrackSig on each chromosome separately.
##
## @param counts binned counts dataframe
## @param i integer denoting which chromosome will be run
## @param activeInSample list of active signatures to tracj
##
## @examples
## breastcancer_results <- foreach(i = c(1:23), .combine = 'c') %dopar% trackParallel(breastcancer_counts, i, breastcancer_sigs)
##
## @name trackParallel

trackParallel <- function (counts, i, activeInSample) {
  if (i == 23) {
    temp <- counts[counts$start_chrom >= i, ]
    temp$temp_bin <- rep(1:nrow(temp))
    traj <- TrackSigCopy(temp, binSize = 200, activeInSample = cervical_sigs, sampleID = "test")

    # assigning changepoints to chr_pos
    # if (!is.null(traj[['changepoints']])) {
    #   for (i in length(traj[['changepoints']])) {
    #     traj[['changepoints']][i] <- traj[['binData']]$chr_pos[traj[['binData']]$temp_bin==traj[['changepoints']][i]]
    #   }
    #
    # }

    # assigning changepoints to bin
    if (!is.null(traj[['changepoints']])) {
      for (i in length(traj[['changepoints']])) {
        traj[['changepoints']][i] <- traj[['binData']]$bin[traj[['binData']]$temp_bin==traj[['changepoints']][i]]
      }

    }
  }
  else {
    temp <- counts[counts$start_chrom == i, ]
    temp$temp_bin <- rep(1:nrow(temp))
    traj <- TrackSigCopy(temp, binSize = 200, activeInSample = activeInSample, sampleID = "test")

    # assigning changepoints to chr_pos
    # if (!is.null(traj[['changepoints']])) {
    #   for (i in length(traj[['changepoints']])) {
    #     traj[['changepoints']][i] <- traj[['binData']]$chr_pos[traj[['binData']]$temp_bin==traj[['changepoints']][i]]
    #   }
    #
    # }

    # assigning changepoints to bin
    if (!is.null(traj[['changepoints']])) {
      for (i in length(traj[['changepoints']])) {
        traj[['changepoints']][i] <- traj[['binData']]$bin[traj[['binData']]$temp_bin==traj[['changepoints']][i]]
      }

    }
  }
  return (traj)
}

## \code{combineTraj} Merge trajectories for all chromosomes into a single trajectory.
##
## @param traj list object containing results for all chromosomes
##
## @examples
## breastcancer_traj <- combineTraj(breastcancer_results)
##
## @name combineTraj

combineTraj <- function (traj) {

  mixtures <- traj[[1]]
  sampleID <- traj[[3]]
  changepoints <- traj[[2]]
  binData <- traj[[4]]

  for (i in 5:length(traj)) {
    if (typeof(traj[[i]]) == 'double') {
      if (is.null(dim(traj[[i]])[2])) {
        changepoints <- c(changepoints, traj[[i]])
      }
      else {
        mixtures <- cbind(mixtures, traj[[i]])
      }
    }
    else if (typeof(traj[[i]]) == 'list') {
      binData <- rbind(binData, traj[[i]])
    }
  }

  combined_traj <- list(mixtures = mixtures, changepoints = changepoints, sampleID = sampleID, binData = binData)
  return (combined_traj)
}
