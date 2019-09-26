# AUTHORS: Yulia Rubanova and Nil Sahin
# Modified for package trackSig by Cait Harrigan


# TODO: phiHist plot - can be added on top of trajectory plot or examined alone
# TODO: phiHist plot should be able to stack or excluse >1 ccf if x range is truncated.
phiHist <- function(){
  NULL
}

#' \code{plotTrajectory}
#'
#'
#' @param mixtures mixtures of mixtures output by TrackSig
#' @param phis list of mean phis corresponding to bins in mixtures matrix
#' @param changepoints list of changepoints to mark on the trajectory plot
#' @param linearX logical whether to plot with a linearly spaced x-axis grid, or with ccf values
#' @param anmac logical whether to plot x-axis restricted to ccf space, or use estimated average number of mutant alleles per cell (anmac)
#' @return ggplot object
#'
#' @name plotTrajectory
#' @export
plotTrajectory <- function(mixtures, phis = NULL, changepoints=NULL, linearX = T, anmac = T, ...){

  # mixtures and phis are binned the same way
  assertthat::assert_that(length(phis) == dim(mixtures)[2])

  # phis should be decreasing
  assertthat::assert_that(all(order(phis, decreasing = T) == 1:length(phis)))

  if(!anmac){ # take x-axis as ccf scale

    # ccf is min(1, anmac)
    # truncate x-axis at phi = 1
    truncateSel <- which(phis <= 1)
    phis <- phis[truncateSel]
    mixtures <- mixtures[,truncateSel]

    # change x-axis lable
    xAx <- "Cancer cell fraction"

    # adjust changepoint indexing
    if (!is.null(changepoints)){
      changepoints <- which(truncateSel %in% changepoints)
    }

  }else{ xAx <- "Average number of mutant alleles per cell" }

  # Plotting the change of mutational signature weights during evolution specified as the order of phi
  colnames(mixtures) <- dim(mixtures)[2]:1
  trajectory <- reshape2::melt(mixtures)
  colnames(trajectory) <- c("Signatures", "Bin", "meanPhi")
  trajectory$Bin <- as.numeric(trajectory$Bin)
  trajectory$meanPhi <- as.numeric(trajectory$meanPhi)

  if(!linearX){ # ggplot formatting specific for real scale

    # "real" scale shows ccf densities

    # place labels in a way that depends on bin density
    # take 8 times the smallest spacing (%)
    spacing <- 800 * min(c(NA, phis) - c(phis, NA), na.rm = T)

    ticSel <- seq(1, length(phis), by = spacing)
    ticLab <- rep("", length(phis))
    ticLab[ticSel] <- round(phis, 2)[ticSel]

    # increasing phi by bin
    trajectory$Bin <- phis[trajectory$Bin]
    trajectory$Bin <- trajectory$Bin[length(trajectory$Bin) : 1]

    g <- (  ggplot2::ggplot(data = trajectory)
          + ggplot2::geom_vline(xintercept = phis, alpha = 0.3)
          + ggplot2::aes(x = Bin, y = meanPhi, group = Signatures, color = Signatures)
          + ggplot2::scale_x_reverse(breaks = phis, labels = ticLab)
         )

    # slice changepoints (reverse axis means max to min)
    cpPos <- cbind(phis[changepoints], phis[changepoints + 1])



  }else{ # ggplot formatting specific for linear scale

    ticSel <- seq(1, length(phis), length.out = min(length(phis), 25))
    ticLab <- rep("", length(phis))
    ticLab[ticSel] <- round(phis, 2)[ticSel]

    g <- (  ggplot2::ggplot(data = trajectory)
          + ggplot2::geom_vline(xintercept = 0:(length(phis) + 1), alpha = 0.3)
          + ggplot2::aes(x = Bin, y = meanPhi, group = Signatures, color = Signatures)
          + ggplot2::scale_x_reverse(breaks = length(phis):1, labels = ticLab)
    )

    # slice changepoints (reverse axis means max to min)
    cpPos <- cbind((length(phis):1)[changepoints], (length(phis):1)[changepoints + 1])

  }

  # TODO: have truncate x range as option
  # TODO: adjust text element size, and alpha for repear lines

  # general ggplot formatting
  g <- (   g
           + ggplot2::geom_point()
           + ggplot2::geom_line()
           + ggplot2::theme_bw()
           + ggplot2::theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
           + ggplot2::ylab("Signature Exposure (%)")
           + ggplot2::xlab(xAx)


           # TODO: allow additional input to ggplot
           #+ list(...)
  )


  if (!is.null(changepoints)) {

    for (i in 1:dim(cpPos)[1]) {
      g <- g + ggplot2::annotate("rect", xmax=cpPos[i,1], xmin=cpPos[i,2],
                                 ymin=-Inf, ymax=Inf, alpha=0.3, fill = "black")

    }
  }

  return(g)
}






