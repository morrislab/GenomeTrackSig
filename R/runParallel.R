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

trackParallel <- function (master, i, activeInSample, binSize) {
  if (i == 23) {
    temp <- master[master$start_chrom >= i, ]
    counts <- binningNmut(temp, binSize)
    traj <- TrackSigCopy(counts, binSize = binSize, activeInSample = activeInSample, sampleID = "sample")
  }
  else {
    temp <- master[master$start_chrom == i, ]
    counts <- binningNmut(temp, binSize)
    traj <- TrackSigCopy(counts, binSize = binSize, activeInSample = activeInSample, sampleID = "test")
  }
  return (traj)
}

## \code{cleanTraj} Define bin numbers—and by extension—activity locations and changepoint locations—in
## relation to the entire genome instead of a single chromosome
##
## @param traj output of trackParallelBootstrap(); large list containing mixtures, changepoints, sampleID,
## and binData for each chromosome run through TrackSig
##
## @examples
## breastcancer_results <- foreach(i = c(1:23), .combine = 'c') %dopar% trackParallel(breastcancer_counts, i, breastcancer_sigs)
##
## @name trackParallel

cleanTraj <- function(traj) {
  # add variable to each trajectory's binData representing its bin number
  # in relation to the whole genome
  start_index <- 1
  for (i in 1:length(traj)) {
    if (typeof(traj[[i]]) == 'list') {
      traj[[i]]$actual_bin <- rep(start_index:((nrow(traj[[i]])+start_index)-1))
      start_index <- start_index + nrow(traj[[i]])
    }
  }


  for (i in 1:length(traj)) {
    # change colnames of mixtures to be actual bin number (in relation to whole genome)
    # instead of bin number (in relation to single chromosome)
    if (typeof(traj[[i]]) == 'double') {
      if (!is.null(dim(traj[[i]])[2])) {
        colnames(traj[[i]]) <- c(traj[[i+3]]$actual_bin)
      }
      else {
        # change changepoint locations to actual bin numbers
        cp <- c()
        if (!is.null(traj[[i]])) {
          for (j in traj[[i]]) {
            bin <- traj[[i+2]]$actual_bin[j]
            cp <- c(cp, as.numeric(bin))
          }
          traj[[i]] <- cp
        }
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

  combined_traj[['changepoints']] <- sort(combined_traj[['changepoints']])
  combined_traj[['binData']] <- combined_traj[['binData']] %>%
    dplyr::arrange(dplyr::desc(combined_traj[['binData']]$actual_bin))

  return (combined_traj)
}

