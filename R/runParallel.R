## \code{trackChromosome} Run TrackSig on each chromosome.
##
## @param master dataframe of mutation counts for your sample/s
## @param activeInSample list of active signatures to fit
## @param binSize desired number of mutations per bin
## @param nSamples number of bootstrap samples to take
## @param parallelize logical indicating whether functions should be run in parallel
##
## @examples
## breastcancer_results <- trackChromosome(pooled_breastcancer_counts,
##                                        breastcancer_sigs, 200, 20, TRUE)
##
## @name trackChromosome

trackChromosome <- function(master, activeInSample, binSize, nSamples, parallelize) {
  j <- NULL
  `%dopar%` <- foreach::`%dopar%` # define parallelization functions
  `%do%` <- foreach::`%do%`
  trajectories <- c() # empty list of trajectories to add each bootstrap sample to

  if (parallelize == TRUE) {
    # intialize cluster
    cores <- base::floor(0.9*(parallel::detectCores()))
    myCluster <- parallel::makeCluster(cores, type="FORK")
    doParallel::registerDoParallel(myCluster)
  }

  # bootstrapping = True
  if (nSamples > 0) {
    for (i in 1:nSamples) {
      # parallelize = True
      if (parallelize == TRUE) {
        traj <- foreach::foreach(j = c(1:23), .combine = 'c') %dopar% bootstrapChromosomes(master, j, activeInSample, binSize)
      }
      # parallelize = False
      else {
        traj <- foreach::foreach(j = c(1:23), .combine = 'c') %do% bootstrapChromosomes(master, j, activeInSample, binSize)
      }
      # helper functions to organize trajectories and combine trajectories for each chromosome
      # into a whole genome trajectory
      cleaned_traj <- cleanTraj(traj)
      combined_traj <- combineTraj(cleaned_traj)
      # add to list of bootstrap trajectories
      trajectories <- c(trajectories, base::assign(paste("traj", as.character(i), sep = ""), combined_traj))
    }
  }

  # bootstrapping = False
  else {
    # parallelize = True
    if (parallelize == TRUE) {
      traj <- foreach::foreach(j = c(1:23), .combine = 'c') %dopar% chrLevel(master, j, activeInSample, binSize)
    }
    # parallelize = False
    else {
      traj <- foreach::foreach(j = c(1:23), .combine = 'c') %do% chrLevel(master, j, activeInSample, binSize)
    }
    cleaned_traj <- cleanTraj(traj)
    combined_traj <- combineTraj(cleaned_traj)
  }

  if (parallelize==TRUE) {
    parallel::stopCluster(myCluster)
  }

  if (nSamples > 0) {
    return (trajectories)
  }
  else {
    return (combined_traj)
  }

    # # parallelizing with bootstrap samples
    # if (nSamples > 0) {
    #   for (i in 1:nSamples) {
    #     traj <- foreach::foreach(j = c(1:23), .combine = 'c') %dopar% bootstrapChromosomes(master, j, activeInSample, binSize)
    #     # helper functions to organize trajectories and combine trajectories for each chromosome
    #     # into a whole genome trajectory
    #     cleaned_traj <- cleanTraj(traj)
    #     combined_traj <- combineTraj(cleaned_traj)
    #     # add to list of bootstrap trajectories
    #     trajectories <- c(trajectories, base::assign(paste("traj", as.character(i), sep = ""), combined_traj))
    #   }
    # }
}

## \code{trackGenome} Run TrackSig on the entire genome.
##
## @param master dataframe of mutation counts for your sample/s
## @param activeInSample list of active signatures to fit
## @param binSize desired number of mutations per bin
## @param nSamples number of bootstrap samples to take
## @param parallelize logical indicating whether functions should be run in parallel
##
## @examples
## breastcancer_results <- trackGenome(pooled_breastcancer_counts,
##                                    breastcancer_sigs, 200, 30, TRUE)
## @name trackGenome

trackGenome <- function(master, activeInSample, binSize, nSamples, parallelize) {
  j <- NULL
  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`

  if (parallelize == TRUE) {
    # intialize cluster
    cores <- base::floor(0.9*(parallel::detectCores()))
    myCluster <- parallel::makeCluster(cores, type="FORK")
    doParallel::registerDoParallel(myCluster)
  }

  # bootstrapping = True
  if(nSamples > 0) {
    # parallelize = True
    if (parallelize == TRUE) {
      traj <- foreach::foreach(j = c(1:nSamples), .combine = 'c') %dopar% bootstrapGenome(master, j, activeInSample, binSize)
    }
    # parallelize = False
    else {
      traj <- foreach::foreach(j = c(1:nSamples), .combine = 'c') %do% bootstrapGenome(master, j, activeInSample, binSize)
    }
  }

  # bootstrapping = False
  else {
    # print error message if user tries to parallelize a single task
    if (parallelize == TRUE) {
      print("Error: parallelization requires that a function must be run multiple times. Set parallelize = FALSE or increase nSamples.")
    }
    # parallelize = False
    else {
      traj <- genomeLevel(master, activeInSample, binSize)
    }
  }

  if (parallelize == TRUE) {
    parallel::stopCluster(myCluster)
  }

  return (traj)
}


## \code{genomeLevel} Run TrackSig on the whole genome; no bootstrapping
##
## @param master binned dataframe of mutation counts
## @param activeInSample list of active signatures to fit
## @param binSize desired number of mutations per bin
##
## @examples
## breastcancer_traj <- genomeLevel(pooled_breastcancer_counts,
##                               breastcancer_sigs, 200)
##
## @name genomeLevel

genomeLevel <- function (master, activeInSample, binSize) {
  counts <- binByChrom(master, binSize)
  traj <- TrackSig(counts, binSize, activeInSample)
  return (traj)
}

## \code{chrLevel} Run TrackSig on each chromosome; no bootstrapping
##
## @param master binned dataframe of mutation counts
## @param i integer denoting chromosome number to run
## @param activeInSample list of active signatures to fit
## @param binSize desired number of mutations per bin
##
## @examples
## breastcancer_traj <- chrLevel(pooled_breastcancer_counts, 14,
##                               breastcancer_sigs, 200)
##
## @name chrLevel

chrLevel <- function (master, i, activeInSample, binSize) {
  if (i == 23) { # run sex chromosomes in tandem
    temp <- master[master$start_chrom >= i, ]
    counts <- binningNmut(temp, binSize)

    traj <- TrackSig(counts, binSize = binSize, activeInSample = activeInSample, sampleID = "sample")
  }
  else { # run autosomal chromosomes individually
    temp <- master[master$start_chrom == i, ]
    counts <- binningNmut(temp, binSize)

    traj <- TrackSig(counts, binSize = binSize, activeInSample = activeInSample, sampleID = "test")
  }
  return (traj)
}


## \code{cleanTraj} Define bin numbers—and by extension—activity locations and changepoint locations—in
## relation to the entire genome instead of a single chromosome
##
## @param traj large list containing mixtures, changepoints, sampleID,
## and binData for each chromosome run through TrackSig
##
## @examples
## breastcancer_results <- foreach(i = c(1:23), .combine = 'c') %dopar% trackParallel(breastcancer_counts, i, breastcancer_sigs)
##
## @name cleanTraj

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
  changepoints <- traj[[2]]
  sampleID <- traj[[3]]
  binData <- traj[[4]]

  # add all chromosomal trajectories into a single trajectory for the whole genome
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

