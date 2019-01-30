#load_all_samples.R

#' \code{load_all_samples} Load previously computed mutation count data for all samples
#'
#' @rdname load_all_samples
#' @name load_all_samples
#' @export

load_samples <- function( countsDir,
                          bootstrapDir,
                          outDir = ".",
                          slidingWindow = FALSE){
  # find or create outDir
  if (!file.exists(outDir)) {
    dir.create(outDir, recursive = T)
  }

  # pick up sample names, check passed sample is present
  sel <- grep("([^/]*)\\.phi\\.txt", list.files(countsDir))
  available_samples <- gsub("([^/]*)\\.phi\\.txt","\\1", list.files(countsDir)[sel])

#  if (sampleName %in% available_samples == FALSE){
#    warning(sprintf("Sample %s not found in provided counts directory", sampleName))
#    return(NULL)
#  }

#  examples_group <- TrackSig:::get_examples_group(tumors)

  for (sample in available_samples){
    print(paste0("Loading sample ", sample, " (", which(available_samples == sample), " out of ", length(available_samples), ")"))
    sample <- load_sample(sample, countsDir)
  }
}


#' \code{load_sample} Load previously computed mutation count data for one sample
#'
#' @rdname load_all_samples
#' @export

load_sample <- function(sampleName, countsDir){

  list[tumor_id, vcfData, phis, assigns_phylo_nodes, acronym, dir_name] <- extract_data_for_example(sampleName, countsDir, tumortypes, dir_create = F)

  if (is.null(vcfData)){
    print(paste0("No data read for sample ", sampleName))
    next
  }

  if (nrow(vcfData) == 0){
    print(paste0("Zero rows for sample " , sampleName))
    next
  }

  if (nrow(vcfData) < 6){ # Plots with less than 6 lines of data are meaningless so ignored
    print(paste0("Less than 6 rows per sample ", sampleName))
    next
  }

  if (slidingWindow) {
    # Sliding window approach
    data_method <- "sliding400"
    window_size=400
    shift <- window_size/100
    gap <- 1
    vcf <- get_sliding_window_data(vcfData, shift=shift, gap = gap)
    phis_sliding_window <- get_sliding_window_data(toVerticalMatrix(phis), shift=shift)
    phis_sliding_window <- phis_sliding_window / shift
    phis_for_plot <- phis_sliding_window
  } else {
    shift <- gap <- NULL
    data_method <- "chunk100"
    vcf <- t(vcfData)
    phis_for_plot <- phis_sliding_window <- phis
  }

  purity <- get_sample_purity(sample)
  phis_for_plot <- phis_for_plot / purity

  if (sum(phis_sliding_window < 0.001) > 0.2 * length(phis_sliding_window)){
    phis_for_plot = NULL
  }

  if (!is.null(phis_for_plot)){
    colnames(vcf) <- round(phis_sliding_window, 3)
  } else {
    colnames(vcf) <- NULL
  }

  if (!is.null(assigns_phylo_nodes)){
    assigns_phylo_nodes <-  as.factor(assigns_phylo_nodes)
    assigns_phylo_nodes_levels <- levels(assigns_phylo_nodes)
    assigns_phylo_nodes <- toVerticalMatrix(as.numeric(assigns_phylo_nodes))

    assigns_phylo_nodes_sw <- assigns_phylo_nodes

    for (l in 1:length(assigns_phylo_nodes_levels)){
      assigns_phylo_nodes_sw[assigns_phylo_nodes_sw == l] <- assigns_phylo_nodes_levels[l]
    }

    stopifnot(length(phis_for_plot) == length(assigns_phylo_nodes_sw))
  } else {
    assigns_phylo_nodes_sw = assigns_phylo_nodes
  }

  list[bootstrap_vcfs, bootstrap_phis] <- lapply(extract_bootstrap_data_for_example(sample, bootstrapDir), t)

  return (list(vcfData, vcf, phis, phis_sliding_window, assigns_phylo_nodes,
               assigns_phylo_nodes_sw, acronym, window, tumor_id, phis_for_plot,
               bootstrap_vcfs, bootstrap_phis))
}
