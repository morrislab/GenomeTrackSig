makeGeneDensityPlot <- function(master, binSize, trajectory, chr_level){
  txdb <-TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  all.genes <- base::suppressMessages(GenomicFeatures::genes(txdb))
  chrom_lengths <- stats::setNames(object = all.genes@seqinfo@seqlengths[1:24], all.genes@seqinfo@seqnames[1:24])
  windows <- GenomicRanges::tileGenome(seqlengths = chrom_lengths,
                                       tilewidth = 1e6, cut.last.tile.in.chrom = TRUE)
  dens <- GenomicRanges::countOverlaps(windows, all.genes)

  start_chrom <- bin <- density <- gene_density <- mean_density <- NULL

  # assign bins to 1Mb regions of genome
  binned_master <- getBinNumber(master, binSize)

  # build dataframe of gene density values for each bin

  density_data <- data.frame(start_chrom = binned_master$start_chrom,
                             bin = binned_master$bin, density = dens)

  if (!is.null(trajectory[[4]]$genome_bin)) {
    # if plottting a shuffle trajectory, exclude Y chromosome data
    density_data <- density_data %>%
      dplyr::filter(start_chrom != 24) %>%
      dplyr::group_by(bin) %>%
      dplyr::summarize(mean_density = base::mean(density))

    change_bins <- assignChromosomeBoundsShuffle(trajectory, chr_level)
    chr_breaks <- sort(change_bins)

  }
  else {
    density_data <- density_data %>%
      dplyr::group_by(bin) %>%
      dplyr::summarize(mean_density = base::mean(density))

    change_bins <- assignChromosomeBounds(trajectory, chr_level)
    chr_breaks <- sort(change_bins)
  }

  #apply Min-Max normalization to gene density column
  density_norm <- sapply(density_data[2], min_max_norm)

  density_data <- density_data %>%
    dplyr::mutate(mean_density = density_norm)

  # make area plot of normalized gene density at each bin
  densityPlot <- ggplot2::ggplot(data = density_data, ggplot2::aes(x = bin, y = mean_density)) +
    ggplot2::geom_area(fill = "cornflowerblue") +
    ggplot2::scale_y_continuous(limits=c(0,1)) +
    ggplot2::ylab("Gene Density") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size=6))

  # # add stripes to distinguish chromosomes
  for (i in 1:length(sort(change_bins))-1) {
    if (i %% 2 != 0) {
      densityPlot <- densityPlot + ggplot2::annotate("rect", xmin=sort(change_bins)[i], xmax = sort(change_bins)[i+1],
                                                     ymin=-Inf, ymax=Inf, alpha=0.3, fill='grey')
    }
  }

  if (!is.null(trajectory[[4]]$genome_bin)) {
    densityPlot <- densityPlot + ggplot2::annotate("rect", xmin=chr_breaks[length(chr_breaks)], xmax = max(trajectory[[4]]$genome_bin),
                                                   ymin=-Inf, ymax=Inf, alpha=0.3, fill='grey')
  }

  return (densityPlot)
}


makeGCPlot <- function(master, binSize, trajectory, chr_level) {
  G <- C <- bin <- GC <- meanGC <- NULL

  # assign bins to 1Mb regions of genome
  binned_master <- getBinNumber(master, binSize)

  loci <- GenomicRanges::GRanges(GenomicRanges::seqinfo(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19))
  loci <- loci[loci@seqnames@values %in% loci@seqnames@values[1:24]]

  tiles <- GenomicRanges::tile(x = loci, width = 1000000)
  tiles <- unlist(tiles)

  seqs <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, tiles)
  counts <- Biostrings::alphabetFrequency(seqs, baseOnly=TRUE)
  freqs <- counts/rowSums(counts)

  nuc_freqs <- data.frame(freqs[,1:4]) %>%
    dplyr::mutate(bin = binned_master$bin,
                  GC = G + C) %>%
    dplyr::group_by(bin) %>%
    dplyr::summarize(meanGC = mean(GC))

  change_bins <- assignChromosomeBounds(trajectory, chr_level)
  change_bins <- sort(change_bins)

  gcPlot <- ggplot2::ggplot(data = nuc_freqs, ggplot2::aes(x = bin, y = meanGC)) +
    ggplot2::geom_area(fill = "darkgrey") +
    ggplot2::scale_y_continuous(limits=c(0,1)) +
    ggplot2::ylab("GC content") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank())


  # # add stripes to distinguish chromosomes
  for (i in 1:length(change_bins)) {
    if (i %% 2 != 0) {
      gcPlot <- gcPlot + ggplot2::annotate("rect", xmin=change_bins[i], xmax = change_bins[i+1],
                                 ymin=-Inf, ymax=Inf, alpha=0.3, fill='grey')
    }
  }

  return (gcPlot)
}

# mutation density

#define Min-Max normalization function
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

makeMutDensityPlot <- function(master, binSize, trajectory, chr_level) {
  width <- bin <- mut_density <- NULL

  binData <- trajectory[[4]]
  # calculate raw mutation density at each bin as # mutations / 1 bp
  widths <- binData %>%
    dplyr::mutate(rowsum = rowSums(binData[7:102])) %>%
    dplyr::select(width, rowsum) %>%
    dplyr::mutate(mut_density = rowsum / width)

  # assign the correct data column to bin number-- depends on GenomeTrackSig() parameters
  if (!is.null(binData$genome_bin)) {
    widths$bin <- binData$genome_bin
    change_bins <- assignChromosomeBoundsShuffle(trajectory, chr_level)
    chr_breaks <- sort(change_bins)
  }
  else if (chr_level==TRUE) {
    widths$bin <- binData$actual_bin
    change_bins <- assignChromosomeBounds(trajectory, chr_level)
    chr_breaks <- sort(change_bins)
  }
  else {
    widths$bin <- binData$bin
    change_bins <- assignChromosomeBounds(trajectory, chr_level)
    chr_breaks <- sort(change_bins)
  }

  #apply Min-Max normalization to mut_density column
  density_norm <- sapply(widths[3], min_max_norm)

  widths <- widths %>%
    dplyr::mutate(mut_density = density_norm)

  # plot sparkline of normalized mutation density at each bin
  mutDensityPlot <- ggplot2::ggplot(data = widths, ggplot2::aes(x = bin, y = mut_density)) +
    ggplot2::geom_line(color = "black") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_text(size = 8),
                   axis.text.y = ggplot2::element_text(size=6),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(0,6,0,6)) +
    ggplot2::labs(x = "Chromosome", y = "Normalized Mutation Density, Gene Density")

  # add stripes to distinguish chromosomes
  for (i in 1:length(chr_breaks)-1) {
    if (i %% 2 != 0) {
      mutDensityPlot <- mutDensityPlot + ggplot2::annotate("rect", xmin=chr_breaks[i], xmax = chr_breaks[i+1],
                                 ymin=-Inf, ymax=Inf, alpha=0.3, fill='grey')
    }
  }

  if (!is.null(binData$genome_bin)) {
    mutDensityPlot <- mutDensityPlot + ggplot2::annotate("rect", xmin=chr_breaks[length(chr_breaks)], xmax = max(binData$genome_bin),
                                                         ymin=-Inf, ymax=Inf, alpha=0.3, fill='grey')
  }


  return (mutDensityPlot)
}

binProp <- function(x, data) {
  # determine the proportion of the bin width that a chromatin state makes up
  prop <- (x / data$binwidth)*100
  return (prop)
}

binChromatinStates <- function(master, binSize, chrom_states) {
  bin <- binwidth <- NULL
  # find the distribution of chromatin states within each bin
  binned_master <- getBinNumber(master, binSize)
  chrom_states$bin <- binned_master$bin

  chrom_states <- chrom_states %>%
    dplyr::group_by(bin) %>%
    dplyr::summarize_at(dplyr::vars('0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'), sum)

  chrom_states <- chrom_states %>%
    dplyr::mutate('binwidth' = rowSums(chrom_states[2:17]))

  chrom_states <- chrom_states %>%
    dplyr::mutate_at(dplyr::vars('0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'), .funs = binProp, data = chrom_states)

  # rename variables to be more informative
  chrom_states <- chrom_states %>%
    dplyr::select(-binwidth) %>%
    dplyr::rename("Baseline" = "0", "TssA" = "1", "TssFlnk" = "2", 'TxFlnk' = "3", "Tx" = "4", "TxWk" = "5",
                  "EnhG" = "6", "Enh" = "7", 'ZNF/Rpts' = "8", 'Het' = "9", "TssBiv" = "10", "BivFlnk" = "11",
                  "EnhBiv" = "12", "ReprPC" = "13", "ReprPCWk" = "14", 'Quies' = "15")

  # convert data into long format for plotting
  chrom_states <- chrom_states %>%
    tidyr::pivot_longer(cols = c('Baseline','TssA','TssFlnk','TxFlnk','Tx','TxWk','EnhG','Enh','ZNF/Rpts','Het',
                                 'TssBiv','BivFlnk','EnhBiv','ReprPC','ReprPCWk','Quies'),
                        names_to = 'state', values_to = 'bin_prop')

  return(chrom_states)
}

makeChromStatePlot <- function(master, trajectory, binSize, chr_level,
                               chrom_states) {
  bin <- bin_prop <- state <- NULL
  # find distribution of chromatin states at each bin
  binned_states <- binChromatinStates(master, binSize, chrom_states)

  # check that bins are assigned correctly
  binData <- trajectory[[4]]
  if (!is.null(binData$genome_bin)) {
    binned_states <- binned_states %>%
      dplyr::filter(bin <= max(binData$genome_bin))
    binned_states$bin <- rep(binData$genome_bin, each=16)
  }
  else if (chr_level==TRUE) {
    binned_states$bin <- rep(base::rev(binData$actual_bin), each=16)
  }
  else {
    binned_states$bin <- rep(binData$bin, each=16)
  }


  # set color palette
  myColors <- c("red3", "red2", "chartreuse3", "green4", "darkgreen", "greenyellow", "tan1", "mediumturquoise", "slateblue1", "indianred3", "lightsalmon2", "lightgoldenrod3",
              "gray56", "gray84", "gray98", "dodgerblue")
  # order chromatin states
  binned_states$state <- factor(binned_states$state, levels = c("TssA", 'TssFlnk', "TxFlnk", "Tx", "TxWk", "EnhG", "Enh",
                                                                "ZNF/Rpts", "Het", "TssBiv", "BivFlnk", "EnhBiv", "ReprPC",
                                                                "ReprPCWk", "Quies", "Baseline"))
  # make stacked barchart
  chromatinPlot <- ggplot2::ggplot(data = binned_states, ggplot2::aes(x = bin, y = bin_prop, fill = state))
  chromatinPlot <- chromatinPlot +
    ggplot2::geom_col(position=ggplot2::position_stack()) +
    ggplot2::scale_fill_manual(values = myColors) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Chromosome", y = "Bin Composition (%)", fill = "Chromatin State") +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(0,6,0,6),
                   axis.title.y = ggplot2::element_text(size=8),
                   axis.text.y = ggplot2::element_text(size=6),
                   legend.title = ggplot2::element_text(size=10),
                   legend.text = ggplot2::element_text(size=8))

  return (chromatinPlot)
}

# AUTHOR: Caitlin Timmons
#' Plot the genomic profile of a tumor and other relevant genomic features.
#'
#'
#' @description
#' \code{plotGenomeContext} For each bin in a set of signature
#' mixtures, the mixture is plotted across the genome. Provided changepoints
#' will be highlighted. Normalized mutation density, normalized gene density, GC content,
#' and the distribution of consensus chromatin states (https://github.com/gerstung-lab/tensorsignatures) are also plotted at each bin.
#'
#' @param counts un-binned dataframe of mutation counts for the sample of interest;
#' same dataframe as `counts` parameter in @seealso \link{GenomeTrackSig}
#' @param trajectory a list containing named elements "mixtures",
#' "changepoints", and 'binData'. See @seealso \link{GenomeTrackSig}.
#' @param binSize number of mutations in each bin; same value as
#' `binSize` parameter in @seealso \link{GenomeTrackSig}.
#' @param title string containing desired plot title; default is blank
#' @param chr_level logical whether TrackSig was run on each chromosome
#' separately; default is FALSE.
#' @param cutoff minimum proportion of bootstrap samples that must agree on a
#' changepoint location in order for that changepoint to be plotted; default is
#' 0.
#' @param show logical whether to print the plot; default TRUE.
#' @param chrom_states dataframe of chromatin state information for 1Mb regions
#' throughout the genome.
#' @return ggplot object
#' @import rlang
#' @export

plotGenomeContext <- function(counts, trajectory, binSize, title='', chr_level=F,
                               cutoff=0, show=T, chrom_states = TrackSig:::roadmap_consensus_chromstates) {

  # make individual plot components
  trajPlot <- plotGenomeProfile(trajectory, show=TRUE, chr_level, cutoff)
  gcPlot <- makeGCPlot(counts, binSize, trajectory, chr_level)
  geneDensityPlot <- makeGeneDensityPlot(counts, binSize, trajectory, chr_level)
  mutDensityPlot <- makeMutDensityPlot(counts, binSize, trajectory, chr_level)
  chromatinPlot <- makeChromStatePlot(counts, trajectory, binSize,
                                      chr_level, chrom_states)

  # change default legend positions for final plot
  traj_mod <- trajPlot + ggplot2::theme(legend.position = 'bottom')
  chrom_mod <- chromatinPlot + ggplot2::theme(legend.position = 'none')
  # chrom state legend
  state_legend <- cowplot::get_legend(
    # create some space to the left of the legend
    chromatinPlot + ggplot2::theme(legend.box.margin = ggplot2::margin(0, 0, 0, 12))
  )

  # align plot components
  plots <- cowplot::align_plots(chrom_mod, geneDensityPlot, mutDensityPlot, gcPlot, traj_mod, align='v', axis='tblr')

  # make plot with all genomic features + trajectory
  fullPlot <- cowplot::plot_grid(plots[[1]],
                             plots[[2]],
                             plots[[3]],
                             plots[[4]],
                             plots[[5]],
                             ncol=1,
                             rel_heights = c(1,0.3,0.25,0.25,1.5))

  fullPlot <- cowplot::plot_grid(fullPlot, state_legend, nrow = 1, rel_widths = c(1, .05))

  # add title to plot
  plot_title <- cowplot::ggdraw() +
    cowplot::draw_label(
      title,
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    ggplot2::theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = ggplot2::margin(0, 0, 0, 1)
    )

  fullPlot <- cowplot::plot_grid(
    plot_title, fullPlot,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.05, 1)
  )

  if (show){print(fullPlot)}

  return (fullPlot)
}
