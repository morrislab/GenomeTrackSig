# AUTHORS: Yulia Rubanova and Nil Sahin
# Modified for package trackSig by Cait Harrigan

truncate_to_range <- function(mixtures, range_) {
  warning("Called a depricated function.")
  min = range_[1]
  max = range_[2]

  x <- mixtures
  col_names <- as.numeric(colnames(x))
  to_leave <- which(col_names <= max+0.01 & col_names >= min-0.01)

  x2 <- x[,to_leave, drop=F]
  colnames(x2) <- col_names[to_leave]
  return(list(x2,to_leave))
}


plot_signatures <- function (mixtures, phis = NULL,
                             changepoints=NULL,
                             ytitle = "Signature exposure (%)",
                             xtitle = "Avg number of mutant alleles per cancer cell",

                             sig_colors = NULL) {

  # order phis if not passed in order
  phis <- phis[order(phis)]

  # Weight matrix is edited with appropriate row and column names
  signatures <- rownames(mixtures)
  df <- data.frame(signatures,mixtures)
  col_names <- c("Signatures", phis)


  colnames(df) <- col_names

  # Plotting the change of mutational signature weights during evolution specified as the order of phi
  # The signature with the maximum change is printed in the plot with the annotate function on ggplot (can be removed if unnecessary)

  df.m <- reshape2::melt(df, id.vars = "Signatures")
  maxx <- apply(mixtures, 2, which.max)
  maxy <- apply(mixtures ,2,max)


  alpha <- 1


  g <- ggplot2::ggplot(data = df.m, ggplot2::aes(x = variable, y = value , group = Signatures, color = Signatures)) +
    ggplot2::geom_line(alpha=alpha, size=size) +
    ggplot2::geom_point(alpha=alpha, size=size) +
    ggplot2::theme_bw() + ggplot2::theme(text = ggplot2::element_text(size = 20))
    ggplot2::theme(axis.title = ggplot2::element_text(size = 20)) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 15))
    #scale_color_manual(values=COLORS2[-c(4,9)])

  if (!is.null(sig_colors)) {
    g <- g + scale_colour_manual(values=sig_colors)
  }

  if (!is.null(phis))
  {
    breaks = as.numeric(col_names[-1])
    labels=paste(round(phis,2),sep = "\n")

    # if there are too many time points, display only every other value
    if (length(col_names) > 50) {
        # Display every forth label
        labels[seq(2,length(labels),4)] <- ""
        labels[seq(3,length(labels),4)] <- ""
        labels[seq(4,length(labels),4)] <- ""
    } else if (length(col_names) > 25) {
        # Display every third label
        labels[seq(2,length(labels),3)] <- ""
        labels[seq(3,length(labels),3)] <- ""
    }
    g <- g + ggplot2::scale_x_discrete(breaks = breaks, labels=labels)
  }


  if (length(changepoints) > 0) {
    for (i in 1:length(changepoints)) {
      g <- g +  ggplot2::annotate("rect", xmax=changepoints[i]-1,
          xmin=changepoints[i], ymin=-Inf, ymax=Inf, alpha=0.3)
    }
  }



  return(list(plot = g, data = df))
}


plot_signatures_real_scale <- function (dd, plot_name, phis = NULL, fitted_data = NULL,
                                        mark_max_signature=F, mark_change_points=F,
                                        change_points=NULL, error_bars = NULL, save=T,
                                        #ytitle = "Signature exposure (%)",
                                        #xtitle = "Avg number of mutant alleles per cancer cell",
                                        ytitle = NULL,
                                        xtitle = "",
                                        assigns_phylo_nodes = NULL,  transition_points = NULL,
                                        remove_sigs_below = 0, cut_at_range = NULL, sig_colors = NULL) {


  if (!is.null(error_bars) & sum(dim(dd) == dim(error_bars)) != 2) {
    stop("Dimentions of error bar matrix should be the same as dimentions of mixture matrix")
  }

  # order phis if not
  phis <- phis[order(phis)]

  sigs_to_remove <- apply(dd, 1, mean) < remove_sigs_below
  dd <- dd[!sigs_to_remove, ]
  if (!is.null(cut_at_range)) {
    dd <- truncate_to_range(dd, cut_at_range)
  }

  # Weight matrix is edited with appropriate row and column names
  signatures <- rownames(dd)
  df <- data.frame(signatures,dd)
  col_names <- c("Signatures")

  decrement <- as.integer(150/ncol(dd))
  for (n in 1:ncol(dd)) {
    col_names <- append(col_names, 150 - decrement*(n-1))
  }

  #colnames(df) <- col_names
  colnames(df) <- c("Signatures", round(phis,3 ))

  # Plotting the change of mutational signature weights during evolution specified as the order of phi
  # The signature with the maximum change is printed in the plot with the annotate function on ggplot (can be removed if unnecessary)

  df.m <- reshape2::melt(df, id.vars = "Signatures")
  df.m$variable <- sapply(df.m$variable, function(x) as.numeric(toString(x)))
  maxx <- apply(dd, 2, which.max)
  maxy <- apply(dd ,2,max)

  if (!is.null(fitted_data))
  {
    if (is.vector(fitted_data))
    {
      all_phis <- c()
      for (i in 1:length(fitted_data))
      {
        all_phis <- c(all_phis, round(sapply(colnames(fitted_data[[i]]), function(x) as.numeric(toString(x))),3))
      }
      all_phis <- unique(all_phis)
    }
  }

  if (!is.null(assigns_phylo_nodes))
  {
    names(assigns_phylo_nodes) <- round(phis, 3) #col_names[-1]
    tree_clusters = factor(assigns_phylo_nodes[sapply(df.m$variable, toString)])

    order_of_clusters <- c("1", "2", "3", "4")

    for (i in 1:length(levels(tree_clusters)))
    {
      tree_clusters <- factor(tree_clusters)
      if (order_of_clusters[i] %in% levels(tree_clusters))
      {

        if (i == 1)
        {
          current_boundary <- min(as.numeric(names(tree_clusters[tree_clusters == order_of_clusters[i]])))
          phis_to_add <- all_phis[all_phis > current_boundary]
          phis_to_add.t <- rep(order_of_clusters[i], length(phis_to_add))
          names(phis_to_add.t) <- phis_to_add
          tree_clusters <- c(tree_clusters,  phis_to_add.t)
        }
        if (i == length(levels(tree_clusters)))
        {
          phis_to_add <- all_phis[all_phis < current_boundary]
          phis_to_add.t <- rep(order_of_clusters[i], length(phis_to_add))
          names(phis_to_add.t) <- phis_to_add
          tree_clusters <- c(tree_clusters,  phis_to_add.t)
        } else {
          current_boundary_upper <- min(as.numeric(names(tree_clusters[tree_clusters == order_of_clusters[i]])))
          phis_to_add <- all_phis[all_phis > current_boundary]
          phis_to_add.t <- rep(order_of_clusters[i], length(phis_to_add))
          names(phis_to_add.t) <- phis_to_add
          tree_clusters <- c(tree_clusters,  phis_to_add.t)

          current_boundary <- current_boundary_upper
        }
      }
    }

    tree_clusters <- factor(tree_clusters)

    if ("Branch" %in% levels(tree_clusters) &
          "Trunk" %in% levels(tree_clusters))
    {
      tree_clusters <- factor(tree_clusters,
                              levels=c("Trunk", "Branch"), order=T)
    }
    tree_clusters <- tree_clusters[order(as.numeric(names(tree_clusters)), decreasing = T)]
    names(tree_clusters) <- round(as.numeric(names(tree_clusters)),3)
    df.m$tree_clusters <- tree_clusters[sapply(round(df.m$variable,3), toString)]
  }

  alpha <- 1
  #   if(!is.null(fitted_data))
  #   {
  #     alpha <- 0.3
  #   }

  size = 1
  if (!is.null(fitted_data)) {
    if (is.vector(fitted_data)) {
      size = 1.5
    }
  }

  g <- ( ggplot2::ggplot(data = df.m, ggplot2::aes(x = variable, y = value, group = Signatures, color = Signatures))
        + ggplot2::geom_line(alpha=alpha, size=1.7)
        + ggplot2::xlab(xtitle)
        + ggplot2::ylab(ytitle)
        + ggplot2::geom_point(alpha=alpha, size=1.7)
        + ggplot2::theme_bw()
        + ggplot2::theme(text = ggplot2::element_text(size = 5))
        + ggplot2::theme(axis.title = ggplot2::element_text(size = 5))
        + ggplot2::theme(axis.text = ggplot2::element_text(size = 5))
        + ggplot2::theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
  )

  if (!is.null(sig_colors)) {
    g <- g + scale_colour_manual(values=sig_colors)
  }

  #   if (!is.null(phis))
  #   {
  #     g <- g + ggplot2::scale_x_discrete(breaks = as.numeric(col_names[-1]), labels=paste(round(phis,2),sep = "\n"))
  #   }

  g <- g + scale_x_reverse()

  if (mark_max_signature)
  {
    g <- g + ggplot2::annotate("text", x=col_names[-1], y=maxy*1.05, label=rownames(dd)[maxx])
  }

  if (!is.null(fitted_data))
  {
    if (is.vector(fitted_data))
    {
      for (i in 1:length(fitted_data))
      {
        alpha = 0.3

          if (!is.null(cut_at_range)) {
            fitted_data[[i]] <- truncate_to_range(fitted_data[[i]], cut_at_range)
          }
          if (ncol(fitted_data[[i]]) == 0) {
            next
          }

        fitted_data.m <- data.frame(signatures,fitted_data[[i]][!sigs_to_remove,])
        #colnames(fitted_data.m) <- colnames(df)
        colnames(fitted_data.m) <- c("Signatures", round(as.numeric(colnames(fitted_data[[i]])), 3))
        fitted_data.m <- reshape2::melt(cbind(Signatures=df[,1],fitted_data.m), id.vars = "Signatures")
        fitted_data.m$variable <- sapply(fitted_data.m$variable, function(x) as.numeric(toString(x)))
        if (!is.null(assigns_phylo_nodes))
        {
          fitted_data.m$tree_clusters <- tree_clusters[sapply(round(fitted_data.m$variable,3), toString)]
        }
        g <- g + ggplot2::geom_line(data=fitted_data.m, size=0.7, ggplot2::aes(x=variable, y=value, group = Signatures, color = Signatures),  alpha=alpha)
      }
    } else {
      fitted_data.m <- data.frame(signatures,fitted_data[!sigs_to_remove,])
      colnames(fitted_data.m) <- colnames(df)
      fitted_data.m <- reshape2::melt(cbind(Signatures=df[,1],fitted_data.m), id.vars = "Signatures")
      fitted_data.m$variable <- sapply(fitted_data.m$variable, function(x) as.numeric(toString(x)))
      if (!is.null(assigns_phylo_nodes))
      {
        fitted_data.m$tree_clusters <- tree_clusters[sapply(round(fitted_data.m$variable,3), toString)]
      }
      g <- g + ggplot2::geom_line(data=fitted_data.m, size=0.7, ggplot2::aes(x=variable, y=value, group = Signatures, color = Signatures))
    }
  }

  if (mark_change_points)
  {
    if (is.null(change_points))
      stop("Please provide change points to mark in the plot")

    #g <- g + geom_vline(xintercept = change_points, size = 1, show.legend = T)

      # if change points are in the list from various bootstrap runs, show all of them and adjust transparency
    if (length(change_points) > 0) {
      if (class(change_points) == "list") {
          alpha = 0.5/length(change_points)
          change_points <- unlist(change_points)
        } else {
          alpha=0.3
        }
      for (i in 1:length(change_points)) {
        g <- g +  ggplot2::annotate("rect", xmin=round(phis,3 )[change_points[i]-1], fill = "red",
                                    xmax=round(phis,3 )[change_points[i]], ymin=-Inf, ymax=Inf, alpha=alpha)
      }
    }
  }

  if (!is.null(error_bars)) {
    error.min <- dd - error_bars
    error.max <- dd + error_bars

    error.min.m <- data.frame(signatures,error.min)
    colnames(error.min.m) <- colnames(df)
    error.min.m <- reshape2::melt(cbind(Signatures=df[,1],error.min.m), id.vars = "Signatures")

    error.max.m <- data.frame(signatures,error.max)
    colnames(error.max.m) <- colnames(df)
    error.max.m <- reshape2::melt(cbind(Signatures=df[,1],error.max.m), id.vars = "Signatures")

    g <- g + geom_errorbar(ggplot2::aes(ymin=error.min.m$value, ymax=error.max.m$value), width=.2)
  }

  if (!is.null(assigns_phylo_nodes))
  {
    g <- g + facet_grid(. ~ tree_clusters, scales = "free_x", space="fixed")
  }

  if (!is.null(transition_points)) {
    g <- g + geom_vline(xintercept = round(phis,3 )[transition_points], colour="red", size=1.5)
  }

  g <- g + geom_vline(xintercept = c(phis), alpha=0.3)

  # add the below annotate function to print the signature with the maximum change
  # ggplot2::annotate("text", x=which.max(dd[maxx,]), y=max(dd[maxx,])+0.01, label=paste("S",maxx,sep=""))
  # theme(legend.position = "none")

  if (save) {
    suppressWarnings(ggplot2::ggsave(filename = plot_name, width = 14, height=4))
  }

  return(list(plot = g, data = df))




  }


























