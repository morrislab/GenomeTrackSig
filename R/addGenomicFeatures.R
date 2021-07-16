# library(BiocFileCache)
# bfc <- BiocFileCache::BiocFileCache(ask=FALSE)
# K562.hmm.file <- BiocFileCache::bfcrpath(bfc, "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmK562HMM.bed.gz")
# K562.hmm <- regioneR::toGRanges(K562.hmm.file)
# K562.hmm
#
# kp <- karyoploteR::plotKaryotype(cex=2, plot.type = 4)
# kpPlotGenes(kp, data=genes.data, r0=0, r1=0.15, gene.name.cex = 2)
# karyoploteR::kpPlotRegions(kp, K562.hmm, col=K562.hmm$itemRgb, r0=0.22, r1=0.3)



makeDensityPlot <- function(master, binSize, trajectory, chr_level, cutoff){
  txdb <-TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  all.genes <- base::suppressMessages(GenomicFeatures::genes(txdb))
  kp <- karyoploteR::plotKaryotype(genome="hg19")
  kp <- karyoploteR::kpPlotDensity(kp, all.genes)
  kp_data <- kp$latest.plot$computed.values

  binned_master <- getBinNumber(master, binSize)

  density_data <- data.frame(start = kp_data[["windows"]]@ranges@start, density = kp_data[["density"]])

  density_data <- density_data %>%
    dplyr::mutate(start_chrom = binned_master$start_chrom,
                  bin = binned_master$bin) %>%
    dplyr::group_by(bin) %>%
    dplyr::summarize(mean_density = base::mean(density))

  change_bins <- assignChromosomeBounds(trajectory, chr_level)
  crPos <- assignCentromereBounds(trajectory, chr_level)
  cpPos <- assignChangepoints(trajectory, cutoff)

  densityPlot <- ggplot2::ggplot(data = density_data, ggplot2::aes(x = bin, y = mean_density)) +
    ggplot2::geom_area(fill = "cornflowerblue") +
    ggplot2::scale_y_continuous(limits=c(0,max(density_data$mean_density)+10)) +
    ggplot2::ylab("Gene Density") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank()) +
    ggplot2::geom_vline(xintercept = crPos[,1], col = 'lightblue', alpha=0.7) +
    ggplot2::geom_vline(xintercept = crPos[,2], col = "lightblue", alpha=0.7)


  # # add stripes to distinguish chromosomes
  for (i in 1:length(change_bins)) {
    if (i %% 2 != 0) {
      densityPlot <- densityPlot + ggplot2::annotate("rect", xmin=change_bins[i], xmax = change_bins[i+1],
                                           ymin=-Inf, ymax=Inf, alpha=0.3, fill='grey')
    }
  }

  # # add centromere locations
  for (i in 1:dim(crPos)[1]) {
    densityPlot <- densityPlot + ggplot2::annotate("rect", xmax=crPos[i,2], xmin=crPos[i,1],
                                         ymin=-Inf, ymax=Inf, alpha=0.5, fill = "lightblue")

  }

  # # add changepoints to plot
  for (i in 1:dim(cpPos)[1]) {
    densityPlot <- densityPlot + ggplot2::annotate("rect", xmax=cpPos$cpPos2[i], xmin=cpPos$cpPos1[i],
                                         ymin=-Inf, ymax=Inf, alpha=cpPos$prob[i]-.1, fill = "red")
  }

  return (densityPlot)
}





makeGCPlot <- function(master, binSize, trajectory, chr_level, cutoff) {
  binned_master <- getBinNumber(master, binSize)

  loci <- GenomicRanges::GRanges(GenomicRanges::seqinfo(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19))
  loci <- loci[loci@seqnames@values %in% loci@seqnames@values[1:24]]

  tiles <- GenomicRanges::tile(x = loci, width = 1000000)
  tiles <- unlist(tiles)

  seqs <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, tiles)
  counts <- Biostrings::alphabetFrequency(seqs, baseOnly=TRUE)
  freqs <- counts/rowSums(counts)

  binned_master <- getBinNumber(master, binSize)

  nuc_freqs <- data.frame(freqs[,1:4]) %>%
    dplyr::mutate(bin = binned_master$bin,
                  GC = G + C) %>%
    dplyr::group_by(bin) %>%
    dplyr::summarize(meanGC = mean(GC))

  change_bins <- assignChromosomeBounds(trajectory, chr_level)
  crPos <- assignCentromereBounds(trajectory, chr_level)
  cpPos <- assignChangepoints(trajectory, cutoff)

  gcPlot <- ggplot2::ggplot(data = nuc_freqs, ggplot2::aes(x = bin, y = meanGC)) +
    ggplot2::geom_area(fill = "darkgrey") +
    ggplot2::scale_y_continuous(limits=c(0,1)) +
    ggplot2::ylab("GC content") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank()) +
    ggplot2::geom_vline(xintercept = crPos[,1], col = 'lightblue', alpha=0.7) +
    ggplot2::geom_vline(xintercept = crPos[,2], col = "lightblue", alpha=0.7)


  # # add stripes to distinguish chromosomes
  for (i in 1:length(change_bins)) {
    if (i %% 2 != 0) {
      gcPlot <- gcPlot + ggplot2::annotate("rect", xmin=change_bins[i], xmax = change_bins[i+1],
                                 ymin=-Inf, ymax=Inf, alpha=0.3, fill='grey')
    }
  }

  # # add centromere locations
  for (i in 1:dim(crPos)[1]) {
    gcPlot <- gcPlot + ggplot2::annotate("rect", xmax=crPos[i,2], xmin=crPos[i,1],
                                         ymin=-Inf, ymax=Inf, alpha=0.5, fill = "lightblue")

  }

  # # add changepoints to plot
  for (i in 1:dim(cpPos)[1]) {
    gcPlot <- gcPlot + ggplot2::annotate("rect", xmax=cpPos$cpPos2[i], xmin=cpPos$cpPos1[i],
                               ymin=-Inf, ymax=Inf, alpha=cpPos$prob[i]-.1, fill = "red")
  }

  return (gcPlot)
}

addGenomicFeatures <- function(trajPlot, master, binSize, trajectory, chr_level=F, cutoff=0) {

  gcPlot <- makeGCPlot(master, binSize, trajectory, chr_level, cutoff)
  densityPlot <- makeDensityPlot(master, binSize, trajectory, chr_level, cutoff)

  plotHat <- cowplot::insert_xaxis_grob(trajPlot, gcPlot, position='bottom', height=ggplot2::unit(0.2, 'null'))
  plotHat <- cowplot::insert_xaxis_grob(plotHat, densityPlot, position='top', height=ggplot2::unit(0.3, 'null'))

  grid::grid.draw(plotHat)

  return (plotHat)
}
