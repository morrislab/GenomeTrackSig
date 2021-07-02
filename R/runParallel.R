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

    # if (!is.null(traj[['changepoints']])) {
    #   cp <- c()
    #   for (i in traj[['changepoints']]) {
    #     cp <- c(cp, as.numeric(colnames(traj[['mixtures']])[i]))
    #   }
    #   traj[['changepoints']] <- cp
    # }
    #traj[['mixtures']] <- traj[['mixtures']][, ncol(traj[['mixtures']]):1]
  }
  else {
    temp <- master[master$start_chrom == i, ]
    counts <- binningNmut(temp, binSize)
    traj <- TrackSigCopy(counts, binSize = binSize, activeInSample = activeInSample, sampleID = "test")

    # if (!is.null(traj[['changepoints']])) {
    #   cp <- c()
    #   for (i in traj[['changepoints']]) {
    #     cp <- c(cp, as.numeric(colnames(traj[['mixtures']])[i]))
    #   }
    #   traj[['changepoints']] <- cp
    # }
    #traj[['mixtures']] <- traj[['mixtures']][, ncol(traj[['mixtures']]):1]
  }
  return (traj)
}

cleanTraj <- function(traj) {

  start_index <- 1
  for (i in 1:length(traj)) {
    if (typeof(traj[[i]]) == 'list') {
      traj[[i]]$actual_bin <- rep(start_index:((nrow(traj[[i]])+start_index)-1))
      start_index <- start_index + nrow(traj[[i]])
    }
  }

  for (i in 1:length(traj)) {
    if (typeof(traj[[i]]) == 'double') {
      if (!is.null(dim(traj[[i]])[2])) {
        colnames(traj[[i]]) <- c(traj[[i+3]]$actual_bin)
      }
      else {
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
    dplyr::arrange(dplyr::desc(actual_bin))

  return (combined_traj)
}

