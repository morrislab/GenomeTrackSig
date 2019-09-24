# AUTHORS: Yulia Rubanova and Nil Sahin
# Modified for package trackSig by Cait Harrigan


# TODO: phiHist plot - can be added on top of trajectory plot or examined alone
# TODO: phiHist plot should be able to stack or excluse >1 ccf if x range is truncated.

plotTrajectory <- function(mixtures, phis = NULL, changepoints=NULL, linearScale = T, ...){

  # mixtures and phis are binned the same way
  assertthat::assert_that(length(phis) == dim(mixtures)[2])

  # phis should be decreasing
  assertthat::assert_that(all(order(phis, decreasing = T) == 1:length(phis)))

  # Plotting the change of mutational signature weights during evolution specified as the order of phi
  trajectory <- reshape2::melt(mixtures)
  colnames(trajectory) <- c("Signatures", "Bin", "meanPhi")
  trajectory$Bin <- as.numeric(trajectory$Bin)
  trajectory$meanPhi <- as.numeric(trajectory$meanPhi)

  # ggplot formatting speficif for real scale
  if(!linearScale){
    # "real" scale shows ccf densities

    # place labels in a way that depends on bin density
    # take 8 times the smallest spacing (%)
    spacing <- 800 * min(c(NA, phis) - c(phis, NA), na.rm = T)

    ticSel <- seq(1, length(phis), by = spacing)
    ticLab <- rep("", length(phis))
    ticLab[ticSel] <- round(phis, 2)[ticSel]

    trajectory$Bin <- phis[trajectory$Bin]

    g <- (  ggplot2::ggplot(data = trajectory)
          + geom_vline(xintercept = phis, alpha = 0.3)
          + ggplot2::aes(x = Bin, y = meanPhi, group = Signatures, color = Signatures)
          + scale_x_reverse(breaks = phis, labels = ticLab)
         )

  }else{

    ticSel <- seq(1, length(phis), length.out = min(length(phis), 25))
    ticLab <- rep("", length(phis))
    ticLab[ticSel] <- round(phis, 2)[ticSel]

    trajectory$Bin <- trajectory$Bin[length(trajectory$Bin) : 1]

    g <- (  ggplot2::ggplot(data = trajectory)
          + geom_vline(xintercept = 0:(length(phis) + 1), alpha = 0.3)
          + ggplot2::aes(x = Bin, y = meanPhi, group = Signatures, color = Signatures)
          + scale_x_reverse(breaks = length(phis):1, labels = ticLab)
    )

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
           + ggplot2::xlab("Cancer Cell Fraction")


           #+ geom_vline(xintercept = 0:length(phis), alpha=0.3)
           # allow additional input to ggplot
           #+ list(...)
  )

  g

  if (length(changepoints) > 0) {
    for (i in 1:length(changepoints)) {
      g <- g +  ggplot2::annotate("rect", xmax=changepoints[i]-1,
          xmin=changepoints[i], ymin=-Inf, ymax=Inf, alpha=0.3)
    }
  }



  return(list(plot = g, data = df))
}





  }


























