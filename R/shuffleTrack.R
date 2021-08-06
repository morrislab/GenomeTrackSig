# Plotting modifications

# Adapts plotting functions to be compatible with 'shuffled' genomic trajectory.
# Functions are only called when plotting a shuffled trajectory, which is identified
# by a the presence of a column in binData called `genome_bin`.

assignChromosomeBoundsShuffle <- function (traj, chr_level) {

  # Assign chromosome breaks to bins

  if (chr_level==F) {
    binData <- traj[[4]]
    change_bins <- c(binData$genome_bin[1])
    for (i in 2:nrow(binData)-1) {
      if (binData$start_chrom[i] != binData$start_chrom[i+1]) {
        if (binData$end_chrom[i] != binData$start_chrom[i]) {
          change_bins <- c(change_bins, binData$genome_bin[i]+0.5)
        }
        else {
          change_bins <- c(change_bins, binData$genome_bin[i+1])
        }
      }
    }
  }

  else {
    change_bins <- c(1)
    binData <- traj[[4]]
    for (i in nrow(binData):2) {
      if (binData$start_chrom[i] != binData$start_chrom[i-1]) {
        if (binData$end_chrom[i] != binData$start_chrom[i]) {
          change_bins <- c(change_bins, binData$actual_bin[i]+0.5)
        }
        else {
          change_bins <- c(change_bins, binData$actual_bin[i-1])
        }
      }
    }
  }
  return (change_bins)
}



assignCentromereBoundsShuffle <- function (traj, chr_level) {

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
        centromere1 <- c(centromere1, binData_subset$genome_bin[start_id-1])
      }
      else {
        centromere1 <- c(centromere1, binData_subset$genome_bin[start_id])
      }

      if (binData_subset$end[end_id] < centromeres$end[i]) {
        centromere2 <- c(centromere2, binData_subset$genome_bin[end_id+1])
      }
      else {
        centromere2 <- c(centromere2, binData_subset$genome_bin[end_id])
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
plotSpaceTrajectoryShuffle <- function(trajectory, show=TRUE, chr_level=F, cutoff=0) {

  # reshape list of trajectories into format that can be used for plotting
  maps <- reshapeTraj(trajectory)

  # assign chromosome breaks to bins
  change_bins <- assignChromosomeBoundsShuffle(trajectory, chr_level)
  chr_breaks <- sort(change_bins)
  chr_labels <- as.character(c(1:22, "X"))

  # assign centromere locations to bins
  crPos <- assignCentromereBoundsShuffle(trajectory, chr_level)

  # assign changepoint locations to bins
  cpPos <- assignChangepoints(trajectory, cutoff)


  # find mean activities at each bin across all bootstrap samples
  avg_df <- as.data.frame(maps[(length(maps)-2):length(maps)])
  colnames(avg_df) <- c('Signatures', "xBin", "exposure")

  g <- (  ggplot2::ggplot(data = avg_df)
          + ggplot2::aes(x = .data$xBin, y = .data$exposure * 100, group = .data$Signatures, color = .data$Signatures)
          + ggplot2::scale_x_continuous(breaks = chr_breaks, labels = chr_labels)

  )

  # add stripes to distinguish chromosomes
  for (i in 1:length(sort(change_bins))-1) {
    if (i %% 2 != 0) {
      g <- g + ggplot2::annotate("rect", xmin=sort(change_bins)[i], xmax = sort(change_bins)[i+1],
                                 ymin=-Inf, ymax=Inf, alpha=0.3, fill='grey')
    }
  }
  g <- g + ggplot2::annotate("rect", xmin=sort(change_bins)[length(change_bins)], xmax = max(avg_df$xBin),
                             ymin=-Inf, ymax=Inf, alpha=0.3, fill='grey')


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
  if (nrow(cpPos) > 0) {
    for (i in 1:dim(cpPos)[1]) {
      g <- g + ggplot2::annotate("rect", xmax=cpPos$cpPos2[i], xmin=cpPos$cpPos1[i],
                                 ymin=-Inf, ymax=Inf, alpha=cpPos$prob[i], fill = "red")
    }
  }


  if (show){print(g)}

  return(g)
}


