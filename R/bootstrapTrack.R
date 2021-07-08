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


# assignBins <- function(data, binSize) {
#   data <- data %>%
#     dplyr::mutate(rowsum = rowSums(data[,6:101]),
#                   bin = 1) %>%
#     base::subset(select = c(1:5,102,103,6:101))
#
#   bin <- 1
#   sums <- 0
#
#   for (i in 1:nrow(data)) {
#     sums <- sums + data$rowsum[i]
#     data$bin[i] <- bin
#     if (sums >= binSize) {
#       sums <- 0
#       bin <- bin + 1
#     }
#   }
#
#   return (data)
# }

## \code{bootstrapSample} Sample with replacement from the mutations in each bin
##
## @param master binned dataframe of mutation counts
##
## @examples
## breast_bootstrap1 <- bootstrapSample(binCounts)
##
## @name bootstrapSample

bootstrapSample <- function(master) {
  master_counts <- master[,7:102]
  new_counts <- data.frame()
  for (i in 1:nrow(master_counts)) {
    vectors <- c()
    zeros <- c()
    for (j in 1:ncol(master_counts)) {
      vector <- rep(j, each=master_counts[i,j])
      vectors <- c(vectors, vector)
    }
    samples <- as.data.frame(sample(vectors, size=sum(base::colSums(master_counts[i,])), replace=TRUE, prob=NULL))
    colnames(samples) <- c('index')
    samples <- samples %>%
      dplyr::group_by(index) %>%
      dplyr::summarize(n = dplyr::n())

    for (i in 1:96) {
      if (i %in% samples$index == FALSE) {
        zeros <- c(zeros, i)
      }
    }

    zeros <- cbind(zeros, rep(0, each=length(zeros)))
    colnames(zeros) <- c('index', 'n')

    samples_full <- base::rbind(samples, zeros) %>%
      dplyr::arrange(index)

    new_counts <- base::rbind(new_counts, samples_full$n)
  }
  master[,7:102] <- new_counts
  base::remove(new_counts, master_counts, samples, samples_full, zeros, vector, vectors)
  return (master)
}

# miniBoot <- function(data, bin) {
#   minidata <- data[data$bin==bin, 8:103]
#   bootstrap <- minidata[sample(nrow(minidata), nrow(minidata), replace=TRUE, prob=NULL), ]
#   return (bootstrap)
# }
#
# bootstrapSample <- function(data) {
#   bootcounts <- data.frame()
#   bins <- unique(data$bin)
#
#   for (i in bins) {
#     bootcounts <- rbind(bootcounts, miniBoot(data, i))
#   }
#
#   boot_sample <- cbind(data[,1:7], bootcounts)
#
#   return (boot_sample)
# }
#
# genBootstrap <- function(data, binSize) {
#   binned_data <- assignBins(data, binSize)
#   bootstrap_sample <- bootstrapSample(binned_data)
#   return (bootstrap_sample)
# }
#
# aggBins <- function(bootstrap_sample) {
#   dup <- data.table::copy(bootstrap_sample)
#   dup <- dup %>%
#     dplyr::group_by(bin) %>%
#     dplyr::summarize_at(dplyr::vars(C_A_ACA:T_G_TTT), sum) %>%
#     dplyr::select(-bin)
#
#   bootstrap_sample <- bootstrap_sample %>%
#     dplyr::group_by(bin) %>%
#     dplyr::summarize(start_chrom = min(start_chrom),
#                      start = min(start),
#                      end_chrom = max(end_chrom),
#                      end = max(end),
#                      width = sum(width)) %>%
#     base::cbind(dup)
#
#   return (bootstrap_sample)
# }
#
# bootstrapBin <- function(data, binSize) {
#   data <- genBootstrap(data, binSize)
#   data <- aggBins(data)
# }

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

# trackParallelBootstrap <- function (master, i, activeInSample, binSize) {
#   if (i == 23) {
#     temp <- master[master$start_chrom >= i, ]
#     counts <- bootstrapBin(temp, binSize)
#
#     traj <- TrackSigCopy(counts, binSize = binSize, activeInSample = activeInSample, sampleID = "sample")
#   }
#   else {
#     temp <- master[master$start_chrom == i, ]
#     counts <- bootstrapBin(temp, binSize)
#     traj <- TrackSigCopy(counts, binSize = binSize, activeInSample = activeInSample, sampleID = "test")
#
#   }
#   return (traj)
# }
#
# bootstrapActivities <- function(master, activeInSample, binSize, nSamples) {
#   trajectories <- c()
#   for (i in 1:nSamples) {
#     traj <- foreach(j = c(1:23), .combine = 'c') %dopar% trackParallelBootstrap(master, j, activeInSample, binSize)
#     cleaned_traj <- cleanTraj(traj)
#     combined_traj <- combineTraj(cleaned_traj)
#     trajectories <- c(trajectories, base::assign(paste("traj", as.character(i), sep = ""), combined_traj))
#   }
#   return (trajectories)
# }

trackParallelBootstrap <- function (master, i, activeInSample, binSize) {
  if (i == 23) {
    temp <- master[master$start_chrom >= i, ]
    counts <- binningNmut(temp, binSize)
    bootstrap_counts <- bootstrapSample(counts)

    traj <- TrackSig(bootstrap_counts, binSize = binSize, activeInSample = activeInSample, sampleID = "sample")
  }
  else {
    temp <- master[master$start_chrom == i, ]
    counts <- binningNmut(temp, binSize)
    bootstrap_counts <- bootstrapSample(counts)

    traj <- TrackSig(bootstrap_counts, binSize = binSize, activeInSample = activeInSample, sampleID = "test")

  }
  return (traj)
}

bootstrapActivities <- function(master, activeInSample, binSize, nSamples) {
  j <- NULL
  `%dopar%` <- foreach::`%dopar%`
  trajectories <- c()
  for (i in 1:nSamples) {
    traj <- foreach::foreach(j = c(1:23), .combine = 'c') %dopar% trackParallelBootstrap(master, j, activeInSample, binSize)
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
#' @param show logical whether to print the plot
#' @return ggplot object
#'
#' @name plotTrajectory
#' @import rlang
#' @export

plotBootstrap <- function(trajectory, show=TRUE){

  # define centromere locations
  # source for locations: https://learn.familytreedna.com/autosomal-ancestry/universal-dna-matching/centromeres-located-autosomal-chromosome/
  chrom <- c(1:23)
  start <- c(121048139,89877778,89392096,48694599,45793062,58276111,57399776,43186096,46276246,38779433,51382454,
             33199289,13500000,13600000,14100000,34375290,22066299,15365878,24354405,25678968,9935312,9600000,58632012)
  end <- c(142185031,95709245,95019980,52401027,50501672,63405808,61385229,48105627,65360007,42104979,
           56400162,36514518,18400000,19100000,18400000,45007538,23203917,17315223,32455280,29267954,
           13290191,16300000,61632012)
  centromeres <- cbind(chrom,start,end)
  centromeres <- as.data.frame(centromeres)
  centromeres <- centromeres %>%
    dplyr::mutate("mid" = (start + ((end-start)/2)),
                  "len" = (end - start))


  # reshape list of trajectories into format that can be used for plotting
  i <- 1
  count <- 1
  timelines <- list()
  cpPos <- c()
  while (i <= length(trajectory)) {
    if(!is.null(trajectory[[i]])){
      mixtures <- trajectory[[i]]
      changepoints <- trajectory[[i+1]]
      binData <- trajectory[[i+1]]

      cpPos <- c(cpPos, changepoints)

      # input checking
      assertthat::assert_that(!is.null(mixtures), msg = "Could not find mixtures for timeline, please supply through results or mixtures paramter.\n")

      # set the bins to colnames(mixtures)
      bins <- as.numeric(colnames(mixtures))

      # mixtures and bins are binned the same way
      assertthat::assert_that(length(bins) == dim(mixtures)[2],
                              msg = "The mixtures object is mal-specified. Column names should correspond to binned bins.\n")

      # Plotting the change of mutational signature weights across the genome specified as the order of bin
      colnames(mixtures) <- 1:dim(mixtures)[2]
      timeline <- reshape2::melt(mixtures)
      colnames(timeline) <- c("Signatures", "xBin", "exposure")
      timeline$exposure <- as.numeric(timeline$exposure)

      timelines <- append(timelines, timeline)

      i <- i + 4
    }
  }

  #calculate mean activity of signatures at each bin over all bootstrap samples
  ex <- NULL
  for (i in 1:length(timelines)) {
    if (typeof(timelines[[i]])=="double") {
      ex <- rbind(ex, timelines[[i]])
    }
  }

  # append mean activity dataframe to list of bootstrap dataframes
  timelines[[length(timelines)+1]] <- timelines[[1]]
  timelines[[length(timelines)+1]] <- timelines[[2]]
  timelines[[length(timelines)+1]] <- colMeans(ex)


  # assigning chromosome breaks to bins
  change_bins <- c(1)
  binData <- trajectory[[4]]
  for (i in nrow(binData):2) {
    if (binData$start_chrom[i] < binData$start_chrom[i-1]) {
      if (binData$end_chrom[i] > binData$start_chrom[i]) {
        change_bins <- c(change_bins, binData$actual_bin[i]+0.5)
      }
      else {
        change_bins <- c(change_bins, binData$actual_bin[i-1])
      }
    }
  }
  if (length(change_bins) != 24) {
    change_bins <- c(change_bins, binData$actual_bin[1])
  }
  chr_breaks <- change_bins
  chr_labels <- as.character(c(1:22, "X", "Y"))

  # assigning centromere locations to bins
  centromere1 <- c()
  centromere2 <- c()
  for (i in centromeres$chrom) {
    binData_subset <- binData[binData$start_chrom==i, ]
    start_id <- base::which.min(abs(c(binData_subset$start)-centromeres$start[i]))
    end_id <- base::which.min(abs(c(binData_subset$end)-centromeres$end[i]))

    if (binData_subset$start[start_id] > centromeres$start[i]) {
      centromere1 <- c(centromere1, binData_subset$actual_bin[start_id+1])
    }
    else {
      centromere1 <- c(centromere1, binData_subset$actual_bin[start_id])
    }

    if (binData_subset$end[end_id] < centromeres$end[i]) {
      centromere2 <- c(centromere2, binData_subset$actual_bin[end_id-1])
    }
    else {
      centromere2 <- c(centromere2, binData_subset$actual_bin[end_id])
    }
  }

  crPos <- base::cbind(centromere1, centromere2)

  # find mean activities at each bin across all bootstrap samples
  avg_df <- as.data.frame(timelines[(length(timelines)-2):length(timelines)])
  colnames(avg_df) <- c('Signatures', "xBin", "exposure")

  g <- (  ggplot2::ggplot(data = avg_df)
          + ggplot2::aes(x = .data$xBin, y = .data$exposure * 100, group = .data$Signatures, color = .data$Signatures)
          + ggplot2::scale_x_continuous(breaks = chr_breaks, labels = chr_labels)

          )

  # general ggplot formatting
  g <- (   g
           + ggplot2::geom_point()
           + ggplot2::geom_line()
           + ggplot2::geom_vline(xintercept = crPos[,1], col = 'lightblue', alpha=0.7)
           + ggplot2::geom_vline(xintercept = crPos[,2], col = "lightblue", alpha=0.7)
           + ggplot2::theme_bw()
           + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                            panel.grid.minor.x = ggplot2::element_blank())
           + ggplot2::ylab("Signature Exposure (%)")
           + ggplot2::xlab("Chromosome")
           + ggplot2::labs(group="Signatures", color = "Signatures")
  )

  # add bootstrapped activities
  indices <- c()
  for (i in 1:(length(timelines)-3)) {
    if (i %% 3 == 0) {
      indices <- c(indices, i)
    }
  }
  for (i in indices) {
    g <- (g
          + ggplot2::geom_line(ggplot2::aes_string(x = avg_df$xBin, y = (timelines[[i]]*100)), alpha=0.3))
  }

  # add centromere locations
  for (i in 1:dim(crPos)[1]) {
    g <- g + ggplot2::annotate("rect", xmax=crPos[i,2], xmin=crPos[i,1],
                                                                                                                                                                   ymin=-Inf, ymax=Inf, alpha=0.3, fill = "lightblue")

  }

  # add changepoints
  cpPos <- as.data.frame(cpPos)
  colnames(cpPos) <- c('cpPos1')
  cpPos <- cpPos %>%
    dplyr::group_by(cpPos1) %>%
    dplyr::summarize(prob = dplyr::n()/length(indices)) %>%
    dplyr::mutate("cpPos2" = cpPos1+1)
    #dplyr::filter(prob >= 0.5)

  #cpPos <- base::cbind(bins[bins %in% changepoints], bins[bins %in% (changepoints+1)])
  if (!is.null(changepoints)) {
    for (i in 1:dim(cpPos)[1]) {
      g <- g + ggplot2::annotate("rect", xmax=cpPos$cpPos2[i], xmin=cpPos$cpPos1[i],
                                 ymin=-Inf, ymax=Inf, alpha=cpPos$prob[i], fill = "red")

    }
  }

  # add stripes to distinguish chromosomes
  for (i in 1:length(change_bins)) {
    if (i %% 2 != 0) {
      g <- g + ggplot2::annotate("rect", xmin=change_bins[i], xmax = change_bins[i+1],
                             ymin=-Inf, ymax=Inf, alpha=0.3, fill='grey')
      }
  }

  if (show){print(g)}

  return(g)
}














