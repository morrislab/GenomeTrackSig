## \code{poolSamples} Generate dataframe summarizing mutation counts across all samples of a
## given cancer type.
##
## @param archivePath file path to folder containing all cancer samples
## @param typesPath file path to dataframe (csv) containing file IDs and corresponding cancer type
## @param cancerType string denoting cancer type for which all samples will be pooled
##
## @examples
## pooled_breastcancer_counts <- poolSamples(archivePath = '~/Desktop/archive',
##                                           typesPath = '~/Desktop/pcawg_cancer_types.csv',
##                                           cancerType = 'Breast-AdenoCA')
##
## @name poolSamples

poolSamples <- function(archivePath, typesPath, cancerType) {
  # read in dataframe of file IDs with corresponding cancer types
  types <- readr::read_csv(typesPath,
                           col_names = c('type', 'guid'))
  types <- types[2:nrow(types), ]

  # list of sample filenames for desired cancer type
  get_files <- c(types$guid[types$type == cancerType])

  # initialize counts dataframe
  master <- readr::read_csv(as.character(paste(archivePath, "/", get_files[1], ".MBcounts.csv", sep = "")),
                            col_types = readr::cols("seqnames" = "c",
                                                    "strand" = "c",
                                                    .default = "d"))

  # add counts from remaining samples into master file and delete remaining files from memory
  for (i in 2:length(get_files)) {
    temp <- readr::read_csv(as.character(paste(archivePath, "/", get_files[i], ".MBcounts.csv", sep = "")),
                            col_types = readr::cols(seqnames = "c",
                                                    strand = "c",
                                                    .default = "d"))
    master[, 6:101] <- master[, 6:101] + temp[, 6:101]
    base::remove(temp)
  }

  master <- master %>%
    dplyr::mutate(seqnames = dplyr::case_when(seqnames == "X" ~ 23,
                                              seqnames == "Y" ~ 24,
                                              TRUE ~ as.numeric(seqnames)),
                  start_chrom = seqnames,
                  end_chrom = seqnames) %>%
    dplyr::select(-seqnames, -strand) %>%
    base::subset(select = c(100,1,101,2:99))

  readr::write_csv(master, paste(cancerType, "_pooled.csv", sep=""))
  return (master)
}

## \code{binningNmut} Group mutation counts dataframe into bins with specified # of mutations per bin.
##
## @param path file path to sample counts dataframe
## @param binSize desired number of mutations per bin (should be >= 100 for optimal results)
##
## @examples
## binCounts <- binningNmut(path = '~/Desktop/Breast-AdenoCA_pooled.csv',
##                          binSize = 100)
##
## @name binningNmut

assignBins <- function(data, binSize) {
  data <- data %>%
    dplyr::mutate(rowsum = rowSums(data[,6:101]),
                  bin = 1) %>%
    base::subset(select = c(1:5,102,103,6:101))

  bin <- 1
  sums <- 0

  for (i in 1:nrow(data)) {
    sums <- sums + data$rowsum[i]
    data$bin[i] <- bin
    if (sums >= binSize) {
      sums <- 0
      bin <- bin + 1
    }
  }

  return (data)
}

miniBoot <- function(data, bin) {
  minidata <- data[data$bin==bin, 8:103]
  bootstrap <- minidata[sample(nrow(minidata), nrow(minidata), replace=TRUE, prob=NULL), ]
  return (bootstrap)
}

bootstrapSample <- function(data) {
  bootcounts <- data.frame()
  bins <- unique(data$bin)

  for (i in bins) {
    bootcounts <- rbind(bootcounts, miniBoot(data, i))
  }

  boot_sample <- cbind(data[,1:7], bootcounts)

  return (boot_sample)
}

genBootstrap <- function(data, binSize) {
  binned_data <- assignBins(data, binSize)
  bootstrap_sample <- bootstrapSample(binned_data)
  return (bootstrap_sample)
}

aggBins <- function(bootstrap_sample) {
  dup <- data.table::copy(bootstrap_sample)
  dup <- dup %>%
    dplyr::group_by(bin) %>%
    dplyr::summarize_at(dplyr::vars(C_A_ACA:T_G_TTT), sum) %>%
    dplyr::select(-bin)

  bootstrap_sample <- bootstrap_sample %>%
    dplyr::group_by(bin) %>%
    dplyr::summarize(start_chrom = min(start_chrom),
                     start = min(start),
                     end_chrom = max(end_chrom),
                     end = max(end),
                     width = sum(width)) %>%
    base::cbind(dup)

  return (bootstrap_sample)
}

bootstrapBin <- function(data, binSize) {
  data <- genBootstrap(data, binSize)
  data <- aggBins(data)
}

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

trackParallelBootstrap <- function (master, i, activeInSample, binSize) {
  if (i == 23) {
    temp <- master[master$start_chrom >= i, ]
    counts <- bootstrapBin(temp, binSize)

    traj <- TrackSigCopy(counts, binSize = binSize, activeInSample = activeInSample, sampleID = "sample")
  }
  else {
    temp <- master[master$start_chrom == i, ]
    counts <- bootstrapBin(temp, binSize)
    traj <- TrackSigCopy(counts, binSize = binSize, activeInSample = activeInSample, sampleID = "test")

  }
  return (traj)
}

bootstrapActivities <- function(master, activeInSample, binSize, nSamples) {
  trajectories <- c()
  for (i in 1:nSamples) {
    traj <- foreach(j = c(1:23), .combine = 'c') %dopar% trackParallelBootstrap(master, j, activeInSample, binSize)
    cleaned_traj <- cleanTraj(traj)
    combined_traj <- combineTraj(cleaned_traj)
    trajectories <- c(trajectories, base::assign(paste("traj", as.character(i), sep = ""), combined_traj))
  }
  return (trajectories)
}

#' Plot the evolutionary trajectory of a tumour
#'
#' For each bin in a set of signature mixtures, the mixture is plotted ac ross
#' pseudo-time. Provided changepoints will be highlighted.
#'
#'
#' @param trajectory a list containing named elements "mixtures", "changepoints",
#'   and "binData". See @seealso \link{TrackSig}.
#' @param linearX logical whether to plot with a linearly spaced x-axis grid, or
#'   with binned phi values
#' @param anmac logical whether to plot x-axis restricted to ccf space, or use
#'   estimated average number of mutant alleles per cell (anmac)
#' @param show logical whether to print the plot
#' @return ggplot object
#'
#' @name plotTrajectory
#' @import rlang
#' @export

plotBootstrap <- function(trajectory, show=TRUE){
  i <- 1
  count <- 1
  timelines <- list()
  while (i <= length(trajectory)) {
    if(!is.null(trajectory[[i]])){
      mixtures <- trajectory[[i]]
      changepoints <- trajectory[[i+1]]
      binData <- trajectory[[i+1]]

      # input checking
      assertthat::assert_that(!is.null(mixtures), msg = "Could not find mixtures for timeline, please supply through results or mixtures paramter.\n")

      # set the phis to colnames(mixtures) - note: used when anmac = T
      phis <- as.numeric(colnames(mixtures))

      # mixtures and phis are binned the same way100
      assertthat::assert_that(length(phis) == dim(mixtures)[2],
                              msg = "The mixtures object is mal-specified. Column names should correspond to binned phis.\n")

      # Plotting the change of mutational signature weights during evolution specified as the order of phi
      colnames(mixtures) <- 1:dim(mixtures)[2]
      timeline <- reshape2::melt(mixtures)
      colnames(timeline) <- c("Signatures", "xBin", "exposure")
      timeline$exposure <- as.numeric(timeline$exposure)

      timelines <- append(timelines, timeline)

      i <- i + 4
    }
  }

  ex <- NULL
  for (i in 1:length(timelines)) {
    if (typeof(timelines[[i]])=="double") {
      ex <- rbind(ex, timelines[[i]])
    }
  }

  timelines[[length(timelines)+1]] <- timelines[[1]]
  timelines[[length(timelines)+1]] <- timelines[[2]]
  timelines[[length(timelines)+1]] <- colMeans(ex)


  # assigning chromosome breaks to bin
  change_phis <- c(1)
  binData <- trajectory[[4]]
  for (i in nrow(binData):2) {
    if (binData$start_chrom[i] < binData$start_chrom[i-1]) {
      if (binData$end_chrom[i] > binData$start_chrom[i]) {
        change_phis <- c(change_phis, binData$actual_bin[i]+0.5)
      }
      else {
        change_phis <- c(change_phis, binData$actual_bin[i-1])
      }
    }
  }

  chr_breaks <- change_phis
  chr_labels <- as.character(c(1:22, "X", "Y"))

  avg_df <- as.data.frame(timelines[(length(timelines)-2):length(timelines)])
  colnames(avg_df) <- c('Signatures', "xBin", "exposure")

  g <- (  ggplot2::ggplot(data = avg_df)
          + ggplot2::geom_vline(xintercept = chr_breaks, alpha = 0.3, col = "red")
          + ggplot2::aes(x = .data$xBin, y = .data$exposure * 100, group = .data$Signatures, color = .data$Signatures)
          + ggplot2::scale_x_continuous(breaks = chr_breaks, labels = chr_labels)

          )

  # general ggplot formatting
  g <- (   g
          # + ggplot2::geom_point()
           + ggplot2::geom_line()
           + ggplot2::theme_bw()
           + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                            panel.grid.minor.x = ggplot2::element_blank())
           + ggplot2::ylab("Signature Exposure (%)")
           + ggplot2::xlab("Chromosome")
           + ggplot2::labs(group="Signatures", color = "Signatures")
  )

  indices <- c()
  for (i in 1:(length(timelines)-3)) {
    if (i %% 3 == 0) {
      indices <- c(indices, i)
    }
  }
  # for (i in indices) {
  #   g <- (g
  #         + geom_line(aes(x = .data$xBin, y = (timelines[[i]]*100)), alpha=0.2))
  # }

  g <- (g
        + geom_line(aes(x = .data$xBin, y = (timelines[[]]*100)), alpha=0.3)
        + geom_line(aes(x = .data$xBin, y = (timelines[[6]]*100)), alpha=0.3)
        + geom_line(aes(x = .data$xBin, y = (timelines[[9]]*100)), alpha=0.3)
        + geom_line(aes(x = .data$xBin, y = (timelines[[12]]*100)), alpha=0.3)
        + geom_line(aes(x = .data$xBin, y = (timelines[[15]]*100)), alpha=0.3))

  if (show){print(g)}

  return(g)
}














