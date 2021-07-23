
readFormat <- function(path) {
  seqnames <- strand <- NULL
  # read in dataframe of file IDs with corresponding cancer types


  # initialize counts dataframe
  master <- readr::read_csv(path,
                            col_types = readr::cols("seqnames" = "c",
                                                    "strand" = "c",
                                                    .default = "d"))

  # make data type consistent in seqnames column
  # split seqnames column into start chrom and end chrom
  # remove unnecessary variables
  master <- master %>%
    dplyr::mutate(seqnames = dplyr::case_when(seqnames == "X" ~ 23,
                                              seqnames == "Y" ~ 24,
                                              TRUE ~ as.numeric(seqnames)),
                  start_chrom = seqnames,
                  end_chrom = seqnames) %>%
    dplyr::select(-seqnames, -strand) %>%
    base::subset(select = c(100,1,101,2:99))

  # save output into csv and return dataframe
  return (master)
}

indivSamples <- function(get_files, i) {
  counts <- readFormat(get_files[i])
  sample_counts <- c(sample_counts, sum(colSums(counts[,6:101])))
  traj <- SpaceTrack(counts, activeInSample = colon_sigs, binSize = 300)
  colorect_trajectories <- c(colorect_trajectories, traj)
  base::saveRDS(colorect_trajectories, file = "colorect_trajectories3.rds")

  return (colorect_trajectories)
}

chromosomalActivities <- function(traj, mixtures, sigs, sample) {
  change_bins <- assignChromosomeBounds(traj, chr_level=FALSE)
  activity_means <- data.frame()
  for (i in 1:length(change_bins)) {
    if (i != length(change_bins)) {
      if (change_bins[i+1]-change_bins[i]<=1) {
        activity_means <- rbind(activity_means, mixtures[,i])
      }
      else {
        rowmeans <- rowMeans(mixtures[,change_bins[i]:change_bins[i+1]-1])
        activity_means <- rbind(activity_means, rowmeans)
      }
    }
    else {
      if (change_bins[i] != ncol(mixtures)) {
        rowmeans <- rowMeans(mixtures[,floor(change_bins[i]):ncol(mixtures)])
        activity_means <- rbind(activity_means, rowmeans)
      }
      else {
        activity_means <- rbind(activity_means, mixtures[,ncol(mixtures)])
      }
    }
  }
  colnames(activity_means) <- sigs
  activity_means <- activity_means %>%
    dplyr::mutate("chromosome" = c(1:24),
                  "sample" = as.integer(sample))

  return (activity_means)
}
