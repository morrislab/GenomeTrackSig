# trajectory <- traj
#
#
#   if(!is.null(trajectory)){
#     mixtures <- trajectory[["mixtures"]]
#     changepoints <- trajectory[["changepoints"]]
#     binData <- trajectory[["binData"]]
#   }
#
#   # input checking
#   assertthat::assert_that(!is.null(mixtures), msg = "Could not find mixtures for timeline, please supply through results or mixtures paramter.\n")
#
#   # set the phis to colnames(mixtures) - note: used when anmac = T
#   phis <- as.numeric(colnames(mixtures))
#
#   # mixtures and phis are binned the same way
#   assertthat::assert_that(length(phis) == dim(mixtures)[2],
#                           msg = "The mixtures object is mal-specified. Column names should correspond to binned phis.\n")
#
#   # phis should be increasing
#   assertthat::assert_that(all(order(phis, decreasing = F) == 1:length(phis)),
#                           msg = "The mixtures object is mal-specified. Binned phis (column names) should be in increasing order.\n")
#
#   if(!anmac){ # take x-axis as ccf scale
#
#     # ccf is min(1, anmac)
#     # truncate x-axis at phi = 1
#     # truncateSel <- which(phis <= 1)
#     # phis <- phis[truncateSel]
#     # mixtures <- mixtures[,truncateSel,drop = FALSE]
#
#     # change x-axis label
#     xAx <- "Locus (chromosome.position)"
#
#     # adjust changepoint indexing
#     # if (!is.null(changepoints)){
#     #   changepoints <- which(truncateSel %in% changepoints)
#     }
#
#   }else{ xAx <- "Average number of mutant alleles per cell" }
#
#   # Plotting the change of mutational signature weights during evolution specified as the order of phi
#   colnames(mixtures) <- 1:dim(mixtures)[2]
#   timeline <- reshape2::melt(mixtures)
#   colnames(timeline) <- c("Signatures", "xBin", "exposure")
#   timeline$xBin <- as.numeric(timeline$xBin)
#   timeline$exposure <- as.numeric(timeline$exposure)
#
#
#   if(!linearX){ # ggplot formatting specific for non-linear scale
#
#     # non-linear scale shows ccf densities
#
#     # place labels in a way that depends on bin density
#     # take 8 times the smallest spacing (%)
#     spacing <- 800 * min(c(NA, phis) - c(phis, NA), na.rm = T)
#
#     ticSel <- seq(1, length(phis), by = spacing)
#     ticLab <- rep("", length(phis))
#     ticLab[ticSel] <- round(phis, 2)[ticSel]
#
#     # increasing phi by bin
#     timeline$xBin <- phis[timeline$xBin]
#     timeline$xBin <- timeline$xBin[length(timeline$xBin) : 1]
#
#     g <- (  ggplot2::ggplot(data = timeline)
#             + ggplot2::geom_vline(xintercept = phis, alpha = 0.3)
#             + ggplot2::aes(x = .data$xBin, y = .data$exposure * 100,
#                            group = .data$Signatures, color = .data$Signatures)
#             + ggplot2::scale_x_reverse(breaks = phis, labels = ticLab)
#     )
#
#     # slice changepoints (reverse axis means max to min)
#     cpPos <- base::cbind(phis[changepoints], phis[changepoints + 1])
#
#   }else{ # ggplot formatting specific for linear scale
#
#     ticSel <- seq(1, length(phis), length.out = min(length(phis), 25))
#     ticLab <- rep("", length(phis))
#     ticLab[ticSel] <- round(phis, 2)[ticSel]
#
#     g <- (  ggplot2::ggplot(data = timeline)
#             + ggplot2::geom_vline(xintercept = 0:(length(phis) + 1), alpha = 0.3)
#             + ggplot2::aes(x = .data$xBin, y = .data$exposure * 100, group = .data$Signatures, color = .data$Signatures)
#             + ggplot2::scale_x_reverse(breaks = length(phis):1, labels = ticLab)
#     )
#
#     # slice changepoints (reverse axis means max to min)
#     cpPos <- base::cbind((length(phis):1)[changepoints], (length(phis):1)[changepoints + 1])
#
#   }
#
#   # TODO: have truncate x range as option
#   # TODO: adjust text element size, and alpha for repear lines
#
#   # general ggplot formatting
#   g <- (   g
#            + ggplot2::geom_point()
#            + ggplot2::geom_line()
#            + ggplot2::theme_bw()
#            + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
#                             panel.grid.minor.x = ggplot2::element_blank())
#            + ggplot2::ylab("Signature Exposure (%)")
#            + ggplot2::xlab(xAx)
#   )
#
#   # add changepoints
#   if (!is.null(changepoints)) {
#
#     for (i in 1:dim(cpPos)[1]) {
#       g <- g + ggplot2::annotate("rect", xmax=cpPos[i,1], xmin=cpPos[i,2],
#                                  ymin=-Inf, ymax=Inf, alpha=0.3, fill = "black")
#
#     }
#   }
#
#   if (show){print(g)}
#
#   return(g)
# }
#
#
# #TODO: ggplot version of this funciton that can be passed to addPhiHist
# plotChangepointChoice <- function(trajectory){
#
#   nMut <- dim(trajectory$binData)[1]
#   nBin <- trajectory$binData$bin[nMut]
#   binSize <- sum(trajectory$binData$bin == 1)
#
#   potentialCps <- binSize * 1:(nMut/binSize)
#
#   graphics::plot(trajectory$binData$phi, ylab = "Empirical Phi", xlab = "Mutation Index",
#                  main = "Potential changepoints in black, selected changepoints in red")
#
#   graphics::abline(v = potentialCps)
#
#   if(!is.null(trajectory$changepoints)){
#     graphics::abline(v = potentialCps[trajectory$changepoints], col = 2)
#   }
#
# }
