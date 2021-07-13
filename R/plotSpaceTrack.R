# AUTHORS: Yulia Rubanova and Nil Sahin
# Modified for package trackSig by Cait Harrigan


##### Plotting helper functions:
##### reshapeTraj(), assignChromosomeBounds(), assignCentromereBounds(),
##### assignChangepoints()

reshapeTraj <- function(traj) {

  ### reshape list of trajectories into format that can be used for plotting
  i <- 1
  maps <- list()

  while (i <= length(traj)) {
    if(!is.null(traj[[i]])){
      mixtures <- traj[[i]]
      binData <- traj[[i+1]]

      # input checking
      assertthat::assert_that(!is.null(mixtures), msg = "Could not find mixtures for timeline, please supply through results or mixtures paramter.\n")

      # set the bins to colnames(mixtures)
      bins <- as.numeric(colnames(mixtures))

      # mixtures and bins are binned the same way
      assertthat::assert_that(length(bins) == dim(mixtures)[2],
                              msg = "The mixtures object is mal-specified. Column names should correspond to binned bins.\n")

      # Plotting the change of mutational signature weights across the genome specified as the order of bin
      colnames(mixtures) <- 1:dim(mixtures)[2]
      map <- reshape2::melt(mixtures)
      colnames(map) <- c("Signatures", "xBin", "exposure")
      map$exposure <- as.numeric(map$exposure)

      maps <- append(maps, map)

      i <- i + 4
    }
  }

  #calculate mean activity of signatures at each bin over all bootstrap samples
  make_means <- NULL
  for (i in 1:length(maps)) {
    if (typeof(maps[[i]])=="double") {
      make_means <- rbind(make_means, maps[[i]])
    }
  }

  # append mean activity dataframe to list of bootstrap dataframes
  maps[[length(maps)+1]] <- maps[[1]]
  maps[[length(maps)+1]] <- maps[[2]]
  maps[[length(maps)+1]] <- colMeans(make_means)

  return (maps)
}

assignChromosomeBounds <- function (traj, chr_level) {

  # Assign chromosome breaks to bins

  if (chr_level==F) {
    change_bins <- c(1)
    binData <- traj[[4]]
    for (i in 2:nrow(binData)-1) {
      if (binData$start_chrom[i] < binData$start_chrom[i+1]) {
        if (binData$end_chrom[i] > binData$start_chrom[i]) {
          change_bins <- c(change_bins, binData$bin[i]+0.5)
        }
        else {
          change_bins <- c(change_bins, binData$bin[i+1])
        }
      }
    }
    if (length(change_bins) != 24) {
      change_bins <- c(change_bins, binData$bin[nrow(binData)])
    }
  }

  else {
    change_bins <- c(1)
      binData <- traj[[4]]
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
  }
  return (change_bins)
}


assignCentromereBounds <- function (traj, chr_level) {

  # assign centromere locations to bins

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

  binData <- traj[[4]]
  # assigning centromere locations to bins
  centromere1 <- c()
  centromere2 <- c()

  if (chr_level == F) {
    for (i in centromeres$chrom) {
      binData_subset <- binData[binData$start_chrom==i, ]
      start_id <- base::which.min(abs(c(binData_subset$start)-centromeres$start[i]))
      end_id <- base::which.min(abs(c(binData_subset$end)-centromeres$end[i]))

      if (binData_subset$start[start_id] > centromeres$start[i]) {
        centromere1 <- c(centromere1, binData_subset$bin[start_id-1])
      }
      else {
        centromere1 <- c(centromere1, binData_subset$bin[start_id])
      }

      if (binData_subset$end[end_id] < centromeres$end[i]) {
        centromere2 <- c(centromere2, binData_subset$bin[end_id+1])
      }
      else {
        centromere2 <- c(centromere2, binData_subset$bin[end_id])
      }
    }
  }

  else {
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
  }


  crPos <- base::cbind(centromere1, centromere2)

  return (crPos)

}

assignChangepoints <- function (traj, cutoff) {
  # assign changepoint locations to bins and calculate proportion of bootstrap
  # samples that agree on each changepoint

  i <- 2
  count <- 0
  cpPos <- c()
  cpPos1 <- cpPos2 <- prob <- NULL

  while (i <= length(traj)) {
    cpPos <- c(cpPos, traj[[i]])
    i <- i + 4
    count <- count + 1
  }

  cpPos <- as.data.frame(cpPos)
  colnames(cpPos) <- c('cpPos1')
  cpPos <- cpPos %>%
    dplyr::group_by(cpPos1) %>%
    dplyr::summarize(prob = dplyr::n()/count) %>%
    dplyr::mutate("cpPos2" = cpPos1+1) %>%
    dplyr::filter(prob >= cutoff)

  return (cpPos)
}

# Plotting using helper functions


#' Plot the genomic trajectory of a tumor.
#'
#' @description
#' \code{plotSpaceTrajectory} For each bin in a set of signature
#' mixtures, the mixture is plotted across the genome. Provided changepoints
#' will be highlighted.
#'
#' @param trajectory a list containing named elements "mixtures",
#' "changepoints", and 'binData'. See @seealso \link{TrackSig}.
#' @param show logical whether to print the plot.
#' @param chr_level logical whether TrackSig was run on each chromosome
#' separately; default is FALSE.
#' @param cutoff minimum proportion of bootstrap samples that must agree on a
#' changepoint location in order for that changepoint to be plotted; default is
#' 0.
#' @return ggplot object
#' @import rlang
#' @export

plotSpaceTrajectory <- function(trajectory, show=TRUE, chr_level=F, cutoff=0) {

  # reshape list of trajectories into format that can be used for plotting
  maps <- reshapeTraj(trajectory)

  # assign chromosome breaks to bins
  change_bins <- assignChromosomeBounds(trajectory, chr_level)
  chr_breaks <- change_bins
  chr_labels <- as.character(c(1:22, "X", "Y"))

  # assign centromere locations to bins
  crPos <- assignCentromereBounds(trajectory, chr_level)

  # assign changepoint locations to bins
  cpPos <- assignChangepoints(trajectory, cutoff)

  # find mean activities at each bin across all bootstrap samples
  avg_df <- as.data.frame(maps[(length(maps)-2):length(maps)])
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

  # add bootstrapped activities, if applicable
  indices <- c()
  for (i in 1:(length(maps)-3)) {
    if (i %% 3 == 0) {
      indices <- c(indices, i)
    }
  }
  for (i in indices) {
    g <- (g
          + ggplot2::geom_line(ggplot2::aes_string(x = avg_df$xBin, y = (maps[[i]]*100)), alpha=0.3))
  }


  # add centromere locations
  for (i in 1:dim(crPos)[1]) {
    g <- g + ggplot2::annotate("rect", xmax=crPos[i,2], xmin=crPos[i,1],
                               ymin=-Inf, ymax=Inf, alpha=0.3, fill = "lightblue")

  }

  # add changepoints to plot
  for (i in 1:dim(cpPos)[1]) {
    g <- g + ggplot2::annotate("rect", xmax=cpPos$cpPos2[i], xmin=cpPos$cpPos1[i],
                               ymin=-Inf, ymax=Inf, alpha=cpPos$prob[i], fill = "red")
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
