#load_counts.R

#' \code{load_counts} Load data from mut_types, phi, quadratic phi files for one sample. \cr
#' Combines the functionality of extract_data_for_example and load_sample
#'
#'
#' @rdname load_counts
#' @export

load_counts <- function

(phiFile, quadphiFile,

 example, dir_counts, tumortypes, dir_results = TrackSig.options()$DIR_RESULTS, dir_create = T) {

  phiData <- tryCatch(read.table(phiFile), error=function(e) NULL) # 96 trinucleotide counts are read as input
  quadphiData <- tryCatch(read.table(quadphiFile)$V1, error=function(e) NULL)

  tumor_id <- gsub("^(.*)\\..*", "\\1", example)

  acronym <- toString(tumortypes[tumortypes$ID == tumor_id,]$tumor_type)
  dir_name <- paste0(dir_results, acronym, "/", tumor_id, "/")

  if (dir_create) {
    suppressWarnings(dir.create(dir_name, recursive = T))
  }

  if (is.null(phiData))
  {
    return(list(tumor_id, phiData, NULL, NULL, NULL, acronym, dir_name))
  }

  assigns_phylo_nodes <- NULL
  if (file.exists(paste0(TrackSig.options()$mutation_assignments, "/", example, ".tree_clusters_by100.txt")))
  {
    assigns_phylo_nodes <- tryCatch(read.delim(paste0(TrackSig.options()$TrackSig.options()$mutation_assignments, "/", example, ".tree_clusters_by100.txt"),
                                               header=F), error=function(e) NULL)
  }



  # Lydia's data
  if (exists("mutation_order"))
  {
    if (file.exists(paste0(TrackSig.options()$mutation_assignments, "/", example, "_ssms.txt")) &
        file.exists(paste0(mutation_order, "/", example, ".mut_order.txt")))
    {

      assigns_phylo_nodes <- tryCatch(read.delim(paste0(TrackSig.options()$mutation_assignments, "/", example, "_ssms.txt"),
                                                 header=T, stringsAsFactors=F), error=function(e) NULL)

      splitted <- sapply(assigns_phylo_nodes[,2], function(x) {strsplit(x, ":", fixed=FALSE)})
      chr <- gsub("chr([\\d]*)", "\\1", sapply(splitted, function(x) {x[2]}))
      pos <- sapply(splitted, function(x) {x[3]})
      cluster <- assigns_phylo_nodes[,4]

      assigns_phylo_nodes <- cbind(chr = chr, pos = pos, cluster = cluster)

      mut_order <- tryCatch(read.delim(paste0(mutation_order, "/", example, ".mut_order.txt"),
                                       header=F, stringsAsFactors=F), error=function(e) NULL)
      mut_order <- paste(mut_order[,1], mut_order[,2])

      assigns_phylo_nodes_pos <- paste(assigns_phylo_nodes[,1], assigns_phylo_nodes[,2])
      assigns_phylo_nodes <- assigns_phylo_nodes[match(mut_order, assigns_phylo_nodes_pos),3]

      n_hundreds <- length(assigns_phylo_nodes) %/% 100
      clusters_ordered_by_100 = c()
      for (i in 1:n_hundreds)
      {
        clusters_100 <- assigns_phylo_nodes[((i-1) * 100): min(length(assigns_phylo_nodes), (i * 100))]
        max_cluster <- names(which.max(table(clusters_100)))
        clusters_ordered_by_100 <- c(clusters_ordered_by_100, max_cluster)
      }
      assigns_phylo_nodes <- toVerticalMatrix(as.factor(clusters_ordered_by_100))
    }
  }

  phiData[,1] <- NULL # Filenames on the first column are deleted

  if (sum(phiData[,1] %% 1  != 0) > 0)
  {
    # second column represents phi values
    phis <- phiData[,1]
    phiData[,1] <- NULL
  }

  # Hack because of the previous bug in the code that assigned 101 mutations to each time points
  rows_keep <- ( apply(phiData, 1, sum) == (TrackSig.options()$bin_size + 1) ) | ( apply(phiData, 1, sum) == TrackSig.options()$bin_size )
  phiData <- phiData[rows_keep, ]

  if (!is.null(phis))
  {
    phis <- phis[rows_keep]
  }

  if(!is.null(assigns_phylo_nodes))
  {
    assigns_phylo_nodes <- assigns_phylo_nodes[rows_keep,1]
  }

  stopifnot(length(phis) == nrow(phiData))

  #hack!! check how hundred are computed in perl script -- there is off by one error
  if (length(assigns_phylo_nodes) == (length(phis) + 1))
  {
    assigns_phylo_nodes <- assigns_phylo_nodes[-length(assigns_phylo_nodes)]
  }
  if(length(phis) != length(assigns_phylo_nodes))
  {
    assigns_phylo_nodes = NULL
  }

  return(list(tumor_id, phiData, phis, quadphiData, assigns_phylo_nodes, acronym, dir_name))
}

#' \code{load_sample} Load previously computed mutation count data for one sample
#'
#' @rdname load_sample
#' @export

load_sample <- function(sample, countsDir, bootstrapDir, tumortypes){

  list[tumor_id, vcfData,
       phis, quadratic_phis, assigns_phylo_nodes,
       acronym, dir_name] <- extract_data_for_example(example = sample,
                                                      dir_counts = countsDir,
                                                      tumortypes = tumortypes,
                                                      dir_results = "results_signature_trajectories/",
                                                      dir_create = F)


  if (is.null(vcfData)){
    print(paste0("No data read for sample ", sample))
  }

  if (nrow(vcfData) == 0){
    print(paste0("Zero rows for sample " , sample))
  }

  if (nrow(vcfData) < 6){ # Plots with less than 6 lines of data are meaningless so ignored
    print(paste0("Less than 6 rows per sample ", sample))
  }

  shift <- gap <- NULL
  data_method <- "chunk100"
  vcf <- t(vcfData)
  phis_for_plot <- phis
  # TODO: rm phis_sliding_window
  phis_sliding_window <- NULL

  purity <- get_sample_purity(sample)
  phis_for_plot <- phis_for_plot / purity

  if (sum(phis< 0.001) > 0.2 * length(phis)){
    phis_for_plot = NULL
  }

  if (!is.null(phis_for_plot)){
    colnames(vcf) <- round(phis, 3)
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

compute_signatures_for_all_examples <- function(countsDir, bootstrapDir, samples_to_run = c()){
  print("Step 2: computing signature activities")

  # pick up sample names
  sel <- grep("([^/]*)\\.phi\\.txt", list.files(countsDir))
  tumors <- gsub("([^/]*)\\.phi\\.txt","\\1", list.files(countsDir)[sel])

  examples_group <- get_examples_group(tumors)

  if (length(samples_to_run) > 0) {
        # run only these tumour samples
        examples_group <- intersect(examples_group, samples_to_run)
  }

  mutation_types <- trinucleotide_internal #data internal to package
  mutation_types <- paste(mutation_types[,1], mutation_types[,2], mutation_types[,3], sep="_")

  #load annotation
  if (TrackSig.options()$pcawg_format == TRUE){
    list[alex, tumortypes, active_signatures, active_signatures.our_samples] <- load_annotation_pcawg()
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
         bootstrap_vcfs, bootstrap_phis] <- load_sample(example, countsDir, bootstrapDir, tumortypes)

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

    plot_name <- paste0(dir_name, "/", acronym, "_", tumor_id, "_", TrackSig.options()$sig_amount, TrackSig.options()$postfix, ".pdf")

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
