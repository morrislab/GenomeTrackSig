#load_all_samples.R

#' \code{load_all_samples} Load previously computed mutation count data for all samples
#'
#' @rdname load_all_samples
#' @name load_all_samples

load_all_samples <- function(countsDir,
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

  allSamples <- list()

  for (sample in available_samples){
    print(paste0("Loading sample ", sample, " (", which(available_samples == sample), " out of ", length(available_samples), ")"))
    allSamples[sprintf("%s", sample)] <- list(load_sample(sample, countsDir, bootstrapDir, slidingWindow))
  }

  return (allSamples)
}


#' \code{load_sample} Load previously computed mutation count data for one sample
#'
#' @rdname load_all_samples
#' @export

load_sample <- function(sample, countsDir, bootstrapDir, slidingWindow = FALSE){

  list[tumor_id, vcfData,
       phis, quadratic_phis, assigns_phylo_nodes,
       acronym, dir_name] <- extract_data_for_example(example = sample,
                                                      dir_counts = countsDir,
                                                      tumortypes = tumortypes,
                                                      dir_results = "results_signature_trajectories/",
                                                      dir_create = F)


  if (is.null(vcfData)){
    print(paste0("No data read for sample ", sample))
    next
  }

  if (nrow(vcfData) == 0){
    print(paste0("Zero rows for sample " , sample))
    next
  }

  if (nrow(vcfData) < 6){ # Plots with less than 6 lines of data are meaningless so ignored
    print(paste0("Less than 6 rows per sample ", sample))
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

  return (list(vcfData, vcf, phis, quadratic_phis, phis_sliding_window, assigns_phylo_nodes,
               assigns_phylo_nodes_sw, acronym, window, tumor_id, phis_for_plot,
               bootstrap_vcfs, bootstrap_phis))
}

#' \code{compute_signatures_for_all_examples} Compute mutational signatures and \cr
#' find changepoints
#'
#' @rdname compute_mutational_signatures
#' @export

compute_signatures_for_all_examples <- function(countsDir, bootstrapDir){
  print("Step 2: computing signature activities")

  add_early_late_transition = TRUE

  age_signatures <- c("S1", "S5", "L1", "1", "5a", "5b")

  # pick up sample names
  sel <- grep("([^/]*)\\.phi\\.txt", list.files(countsDir))
  tumors <- gsub("([^/]*)\\.phi\\.txt","\\1", list.files(countsDir)[sel])

  examples_group <- get_examples_group(tumors)

  mutation_types <- trinucleotide #data internal to package
  mutation_types <- paste(mutation_types[,1], mutation_types[,2], mutation_types[,3], sep="_")

  #load annotation
  if (TrackSig.options()$pcawg_format == TRUE){

  }
  else{
    list[alex, tumortypes, active_signatures, active_signatures.our_samples] <- load_annotation()
  }

  for (example in examples_group)
  {
    set.seed(which(examples_group == example))
    print(paste0("Example ", example, " (", which(examples_group == example), " out of ", length(examples_group), ")"))

    list[vcfData, vcf, phis, quadratic_phis, phis_sliding_window, assigns_phylo_nodes,
         assigns_phylo_nodes_sw, acronym, window, tumor_id, phis_for_plot,
         bootstrap_vcfs, bootstrap_phis] <- load_sample(example, countsDir, bootstrapDir)

    if (TrackSig.options()$sig_amount == "onlyKnownSignatures") {
      # Fit only known signatures
      list[alex.t, matched_type, acronym] <- get_signatures_for_current_sample(tumor_id, active_signatures.our_samples, alex, TrackSig.options()$noise_sig)
    } else {
      alex.t <- alex
    }

    if (is.null(acronym) || acronym == "") {
      print(paste("ERROR: Cancer type not found for ", example))
      next
    }

    if (is.null(alex.t))
    {
      print(paste0("No active signatures for sample", example, " ...."))
      next
    }

    if (is.vector(alex.t)) {
      next
    }

    if (sum(vcf[apply(alex.t,1,sum) == 0,] != 0) != 0) {
      print(paste0("Sample ", example, ": some trinucleotides have probability 0 under the model, but their count is non-zero. Sot the count vector is impossible under the model."))
      next
    }

    dir_name <- paste0(TrackSig.options()$DIR_RESULTS, acronym, "/", tumor_id, "/")

    suppressWarnings(dir.create(dir_name, recursive = T))

    if (!is.null(phis_for_plot))
    {
      write(phis_for_plot, file=paste0(dir_name, "phis.txt"), ncolumns=length(phis_for_plot))
    }

    method_name <- "iterativeChangePoints"

    if (!file.exists(paste0(dir_name, "mixtures.csv")) || !file.exists(paste0(dir_name, "changepoints.txt")))
    {
      if (TrackSig.options()$changepoint_method == "PELT") {
        list[changepoints, mixtures] <- find_changepoints_pelt(vcf, alex.t, phis, quadratic_phis)
      } else {
        list[bics, optimal, changepoints, mixtures] <- find_changepoints_over_all_signatures_one_by_one(vcf, alex.t, n_signatures = ncol(alex.t))
      }

      write.csv(mixtures, file=paste0(dir_name, "mixtures.csv"))

      n_col <- ifelse(length(changepoints) > 0, length(changepoints), 1)
      write(changepoints, file=paste0(dir_name, "changepoints.txt"), ncolumns=n_col)
    } else {
      mixtures <- read_mixtures(paste0(dir_name, "mixtures.csv"))
      cp_file = paste0(dir_name, "changepoints.txt")
      if (file.info(cp_file)$size == 1) {
        changepoints <- c()
      } else {
        changepoints <- unlist(read.table(cp_file, header=F))
      }
    }

    if (!is.null(assigns_phylo_nodes_sw)) {
      write(assigns_phylo_nodes_sw,  file=paste0(dir_name, "assignments.txt"), ncolumns=length(assigns_phylo_nodes_sw))
    } else  {
      n_clusters = transition_points = assigns_phylo_nodes_sw = NULL
    }

    age_signatures <- intersect(rownames(mixtures), age_signatures)

    if (TrackSig.options()$simulated_data) {
      plot_name <- paste0(dir_name, "/", tumor_id,  ".pdf")
    } else {
      plot_name <- paste0(dir_name, "/", acronym, "_", tumor_id, "_", TrackSig.options()$sig_amount, TrackSig.options()$postfix, ".pdf")
    }

    if (TrackSig.options()$PLOT_FULL_NAME)
    {
      plot_name <- paste0(dir_name, "/", acronym, "_", data_method, "_multMix_fittedPerTimeSlice_", sig_amount, "_noPrior_", method_name, postfix, ".pdf")
    }

    mark_cp <- !is.null(changepoints)
    plot_signatures(mixtures*100, plot_name=plot_name, phis = phis_for_plot, mark_change_points=mark_cp, change_points=changepoints,
                    #assigns_phylo_nodes = assigns_phylo_nodes_sw,
                    transition_points = transition_points,
                    scale=1.2)

    mixtures.rescaled = NULL

    print(paste("Computed example", example))
  }

}
