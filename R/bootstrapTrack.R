## \code{bootstrapSample} Sample with replacement from the mutations in each bin
##
## @param master binned dataframe of mutation counts
## @param binSize desired number of mutations per bin
##
## @examples
## breast_bootstrap1 <- bootstrapSample(binCounts, 200)
##
## @name bootstrapSample

bootstrapSample <- function(counts, binSize) {
  index <- NULL
  # initialize empty dataframe to store the bootstrapped counts for each bin
  new_counts <- data.frame()
  # sample with replacement from each bin
  for (i in 1:nrow(counts)) {
    samples <- data.frame(index = sample(rep(c(1:96), times=counts[i,7:102]),
                                         size=binSize, replace=TRUE, prob=NULL)) %>%
      dplyr::group_by(index) %>%
      dplyr::summarize(n = dplyr::n())
    # add in zero values for mutation types that are not represented in the bootstrap sample
    samples_full <- dplyr::left_join(data.frame(index = c(1:96)), samples, by='index') %>%
      tidyr::replace_na(replace = list(n = 0))
    # add bootstrapped counts to dataframe
    new_counts <- base::rbind(new_counts, samples_full$n)
  }
  # set counts to bootstrapped counts
  counts[,7:102] <- new_counts
  base::remove(new_counts, samples, samples_full)
  return (counts)
}

## \code{bootstrapChromosomes} Framework to take bootstrap samples and run
## TrackSig at the chromosomal level
##
## @param master binned dataframe of mutation counts
## @param i integer denoting chromosome number to run
## @param activeInSample list of active signatures to fit
## @param binSize desired number of mutations per bin
##
## @name bootstrapChromosomes

bootstrapChromosomes <- function (master, i, activeInSample, binSize) {
  if (i == 23) { # run sex chromosomes in tandem
    temp <- master[master$start_chrom >= i, ]
    counts <- binningNmut(temp, binSize)
    bootstrap_counts <- bootstrapSample(counts)

    traj <- TrackSig(bootstrap_counts, binSize = binSize, activeInSample = activeInSample, sampleID = "sample")
  }
  else { # run autosomal chromosomes individually
    temp <- master[master$start_chrom == i, ]
    counts <- binningNmut(temp, binSize)
    bootstrap_counts <- bootstrapSample(counts)

    traj <- TrackSig(bootstrap_counts, binSize = binSize, activeInSample = activeInSample, sampleID = "test")

  }
  return (traj)
}

## \code{bootstrapGenome} Framework to take bootstrap samples and run TrackSig at
## the whole genome level.
##
## @param master dataframe of mutation counts for single or pooled sample/s
## @param activeInSample list of active signatures to fit
## @param binSize desired number of mutations per bin
##
## @examples
## breastcancer_results <- bootstrapGenome(breast_master, activeInSample = c('SBS1', 'SBS3', 'SBS5'), binSize=200, nSamples = 30)
##
## @name trackParallel

bootstrapGenome <- function (master, i, activeInSample, binSize) {
  set.seed(i)
  counts <- binByChrom(master, binSize)
  bootstrap_counts <- bootstrapSample(counts)
  traj <- TrackSig(bootstrap_counts, activeInSample = activeInSample, binSize = binSize, sampleID = "sample")

  return (traj)
}
