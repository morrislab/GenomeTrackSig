# AUTHOR: Yulia Rubanova
# Modified for package TrackSig by Cait Harrigan

#' Helper functions for \code{TrackSig}
#'
#' Non-exported functions called by \code{TrackSig} functions. \cr
#' Not intended for end-user use.
#'
#' @rdname helper_functions
#' @name helper_functions
NULL

#' \code{gg_color_hue} Select colour hues for plotting
#' @rdname helper_functions

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#' \code{list} Succinct object assignment
#' @rdname helper_functions
#'
list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}

#' \code{get_values_from_list} Retrieve values by name from list
#' @rdname helper_functions
get_values_from_list <- function(list, name_of_value, FUN=NULL, default = NULL, concat_function="cbind")
{
  concat_function <- match.fun(concat_function)
  res <- c()
  for (i in 1:length(list))
  {
    if (!is.null(FUN))
    {
      FUN <- match.fun(FUN)
      value <- sapply(list[[i]][names(list[[i]]) == name_of_value], FUN)

      if (is.null(value[[1]]) & !is.null(default))
      {
        value[[1]] <- default
      }
      res <- concat_function(res, value)
    } else
    {
      res <- concat_function(res, unlist(list[[i]][names(list[[i]]) == name_of_value]))
    }
  }

  return(res)
}

#' \code{toVerticalMatrix} Convert a list or vector-like object to a vertial matrix
#' @rdname helper_functions
toVerticalMatrix <- function(L)
{
  if (is.vector(L))
    return(matrix(L, ncol=1))
  else
    return(as.matrix(L))
}

#' \code{reset_parallelize} Select colour hues for plotting
#' @rdname helper_functions
reset_parallize <- function()
{
  closeAllConnections()
  cl <- makeCluster(3)
  registerDoParallel(cl)
}

IgnoreVectorOrMatrix <- function(x, FUN)
{
  if (is.vector(x))
  {
    return(x)
  } else if (is.matrix(x) | is.data.frame(x)) {
    FUN <- match.fun(FUN)
    return(FUN(x))
  } else
  {
    stop(paste("Unknown type of data:", head(x)))
  }
}


#' \code{get_sliding_window_data} <man content>
#' @rdname helper_functions
get_sliding_window_data <- function (data, shift, gap = 1) {
  if (is.null(gap)) {
    gap = 1
  }
  vcf <- c()
  last_iteration = FALSE
  for (i in 1:nrow(data)) {
    start = i * gap
    end = i*gap + shift - 1

    if ((start <= nrow(data)) && (end > nrow(data))) {
      #end = nrow(data)
      #last_iteration = TRUE
      break
    }
    vcf <- rbind(vcf, apply(toVerticalMatrix(data[start:end,]),2,sum))
  }
  vcf <- t(vcf)
  return(vcf)
}


#' \code{get_sliding_window_data} <man content>
#' @rdname helper_functions
merge_data_chunks <- function (vcfData) {
  vcf <- matrix(nrow = nrow(vcfData)/2, ncol = ncol(vcfData))
  for (i in 1:nrow(vcf)) {
    for (j in 1:ncol(vcf)) {
      vcf[i,j] <- vcfData[i*2-1,j] + vcfData[i*2,j]
    }
  }
  vcf <- t(vcf)
  return(vcf)
}

#' \code{get_sliding_window_data} <man content>
#' @rdname helper_functions
get_noise_sig <- function(noise_sig){
  if (is.null(noise_sig))
  {
    return(c())
  }

  if (noise_sig == "uniform")
  {
    noise = rep(1/96, 96)
    return(noise)
  }

  print(paste0("Noise signature type: ", noise_sig, " not found"))
  return(c())
}

#' \code{get_signatures_for_current_sample} <man content>
#' @rdname helper_functions
get_signatures_for_current_sample <- function (id, active_signatures.our_samples, alex, noise_sig) {
  if (grepl("simulation([\\d]*)", id))
  {
    num <- as.numeric(gsub("simulation([\\d]*)", "\\1", id))
    active_sigs <- sapply(unlist(read.table(paste0(DIR_RESULTS, "simulation", num, "/true_mixture_sigs.txt"), header=F)), toString)
    alex.t <- alex[,active_sigs]
    return(list(alex.t, "simulation", "simulation"))
  }

  current_type <- active_signatures.our_samples[active_signatures.our_samples$ID == id,]
  # If the known signature table does not contain this cancer type, skip this tumor sample
  if (nrow(current_type) == 0)
  {
    print(paste("Skipping file", id, "... : no such sample in active sigs table"))
    return(NULL)
  }

  if (is.na(current_type$Name))
  {
    print(paste("Skipping file", id, "with type", current_type$tumor_type, "...: no cancer type name"))
    return(NULL)
  }
  matched_type <- current_type$Name
  signature_indices <- which(as.logical(current_type[,4:ncol(current_type)]))
  if ("SNPs" %in% colnames(alex))
  {
    alex <- alex[,-which(colnames(alex) == "SNPs")]
  }
  alex.t <- toVerticalMatrix(alex[,signature_indices])
  colnames(alex.t) <- colnames(alex)[signature_indices]

  if (!is.null(noise_sig))
  {
    noise_sig_distr <- get_noise_sig(noise_sig)
    alex.t <- cbind(alex.t, noise_sig_distr)
  }
  return(list(alex.t, matched_type, current_type$tumor_type))
}

#' \code{save_data} is depricated. See load_data instead.
#' @rdname helper_functions
save_data <- function(){
  warning("save_data is depricated. See help for more")

  names_trinucleotide <- read.table(paste0(DIR, "trinucleotide.txt"), stringsAsFactors = F)
  names_trinucleotide <- apply(names_trinucleotide, 1, function(x) { do.call("paste", c(as.list(x), sep = "_"))})

  lydia_signatures = F
  pcawg_signatures = T
  cosmic_signatures = F


  # Load the tumor types for tumor IDs.
  tumortypes <- read.delim("tumortypes.PCAWG.txt", header = F, stringsAsFactors=F)
  colnames(tumortypes) <- c("ID", "tumor_type")

  stopifnot(sum(lydia_signatures + pcawg_signatures + cosmic_signatures) == 1)

  if (lydia_signatures) {
    lydia = T

    # Signatures from Lydia Liu (OICR)
    alex <- read.table(paste0(DIR, "/lydia_prostate_cancer_200/NMF_TrinucleotideSignature3wgs.tsv"), header = T)
    rownames(alex) <- names_trinucleotide
    alex <- alex[,-c(1,2)]
    colnames(alex) <- paste0("L", 1:ncol(alex))

    tumors.lydia <- gsub("(.*).phi.txt", "\\1", list.files("lydia_prostate_cancer_200//counts"))
    tumortypes.lydia <- data.frame(ID = tumors.lydia, tumor_type = "PRAD")
  } else if (pcawg_signatures)
  {
    # Signatures from PCAWG group
    alex <- read.csv(paste0(DIR, "PCAWG_signatures.csv"), header = T)
    #alex <- read.csv(paste0(DIR, "PCAWG_signature_patterns_beta.csv"), header = T)
    pcawg_trinucleotides <- paste0(gsub(">", "_", alex[,1]), "_", alex[,2])
    rownames(alex) <- pcawg_trinucleotides
    alex <- alex[,-c(1,2)]
    # to make the order of mutation types match names_trinucleotide
    alex <- alex[match(names_trinucleotide, pcawg_trinucleotides),]
    colnames(alex) <- gsub("Signature.(.*)", "\\1", colnames(alex))
    colnames(alex) <- gsub("PCAWG.(.*)", "\\1", colnames(alex))
  } else if (cosmic_signatures) {
    # ALEX DATA
    # The trinucleotide count matrix will be regressed on the 30 Alexandrov Mutational Signature frequencies
    # http://cancer.sanger.ac.uk/cosmic/signatures
    alex <- read.table(paste0(DIR, "alexSignatures.txt"))
    rownames(alex) <- names_trinucleotide
    colnames(alex) <- paste0("S", 1:ncol(alex))
  }


  if (lydia_signatures) {
    active_signatures = data.frame(Name = "Pancreatic cancer", acronym = "PRAD", toHorizontalMatrix(rep(1,ncol(alex))))
    colnames(active_signatures) <- c("Name", "acronym", paste0("L", 1:ncol(alex)))
    active_signatures.lydia_samples <- merge(tumortypes.lydia, active_signatures, by.x="tumor_type", by.y="acronym", all.x=T)
  } else if (cosmic_signatures)
  {
    # Load active signatures for each tumor type
    active_signatures <- read.delim("active_signatures_transposed.txt", stringsAsFactors=F)
    active_signatures[is.na(active_signatures)] <- 0
    # Removing column "Other.signatures"
    active_signatures <- active_signatures[, -ncol(active_signatures)]

    # Some of the signatures have several types corresponding to them.
    # The rows which correspond to several tumor types and duplicated in the table.
    indices_to_duplicate <- which(grepl(",",active_signatures$acronym))
    while (length(indices_to_duplicate) != 0)
    {
      index <- indices_to_duplicate[1]
      current_type <- active_signatures[index,]
      trim <- function (x) gsub("^\\s+|\\s+$", "", x)
      corresponding_tumor_types <- trim(unlist(strsplit(current_type$acronym, ",")))

      to_add <- current_type
      for (i in 2:length(corresponding_tumor_types))
      {
        to_add <- rbind(to_add, current_type)
      }
      to_add$acronym <- corresponding_tumor_types

      active_signatures <- rbind(active_signatures[1:(index-1),], to_add, active_signatures[(index+1):nrow(active_signatures),])
      indices_to_duplicate <- which(grepl(",",active_signatures$acronym))
    }
    active_signatures.our_samples <- merge(tumortypes, active_signatures, by.x="tumor_type", by.y="acronym", all.x=T)
    # active_signatures.lydia_samples <- merge(tumortypes.lydia, active_signatures, by.x="tumor_type", by.y="acronym", all.x=T)
  } else if (pcawg_signatures)
  {
    active_signatures <- NULL

    active_signatures.our_samples <- read.csv("PCAWG_signatures_in_samples.csv")
    #active_signatures.our_samples <- read.csv("PCAWG_signatures_in_samples_beta.csv")
    colnames(active_signatures.our_samples) <- c("tumor_type",  "ID", colnames(active_signatures.our_samples)[3:ncol(active_signatures.our_samples)])
    active_signatures.our_samples$tumor_type <- gsub("^(.*)-.*$","\\1", active_signatures.our_samples$tumor_type)

    # Making the table binary
    active_signatures.our_samples.data <- active_signatures.our_samples[,3:ncol(active_signatures.our_samples)]
    active_signatures.our_samples.data[active_signatures.our_samples.data > 0] <- 1
    colnames(active_signatures.our_samples.data) <- gsub("Signature.(.*)", "\\1", colnames(active_signatures.our_samples.data))
    colnames(active_signatures.our_samples.data) <- gsub("PCAWG.(.*)", "\\1", colnames(active_signatures.our_samples.data))

    active_signatures.our_samples <- cbind(active_signatures.our_samples[,c(1,2)], Name=NA, active_signatures.our_samples.data)
    active_signatures.our_samples$Name <- sapply(active_signatures.our_samples$tumor_type, toString)
  }

  lydia = FALSE
  if (lydia == TRUE)
  {
    tumortypes = tumortypes.lydia
    active_signatures.our_samples = active_signatures.lydia_samples
    save(alex, tumortypes, active_signatures, active_signatures.our_samples, file="saved_data/annotation_data.lydia.lydias_signatures.RData")
  } else if (pcawg_signatures) {
    save(alex, tumortypes, active_signatures, active_signatures.our_samples, file="saved_data/annotation_data.pcawg_sigs.RData")
  } else {
    save(alex, tumortypes, active_signatures, active_signatures.our_samples, file="saved_data/annotation_data.RData")
  }

}

#' \code{toHorizontalMatrix} <man content>
#' @rdname helper_functions
toHorizontalMatrix <- function(L){
  if (is.vector(L))
    return(matrix(L, nrow=1))
  else
    return(as.matrix(L))
}

#' \code{gather_statistics} <man content>
#' @rdname helper_functions
gather_statistics <- function(mixtures, changepoints, tumor_id, dir_name, tumor_type, mixtures.rescaled=NULL){
  col_names <-  c("tumor_id", "tumor_type", "n_changepoints", "n_total_sigs", "top_signature", "max_change", "sig_with_max_change", "n_sigs_greater_03")
  if (!is.null(mixtures.rescaled))
  {
    col_names <- c(col_names, "top_signature_rescaled", "max_change_rescaled", "sig_with_max_change_rescaled", "n_sigs_greater_03_rescaled")
  }
  summary <- data.frame(matrix(nrow=1, ncol=length(col_names)))
  colnames(summary) <- col_names

  summary$tumor_id <- tumor_id
  summary$tumor_type <- tumor_type
  summary$n_changepoints <- length(changepoints)
  summary$n_total_sigs <- nrow(mixtures)

  mean_sigs <- apply(mixtures, 1, mean)
  max_change_sigs <- apply(mixtures,1,max) - apply(mixtures,1,min)
  sd_sigs <- apply(mixtures, 1, sd)

  summary$top_signature <- names(which.max(mean_sigs))
  summary$max_change <- max(max_change_sigs)
  summary$sig_with_max_change <- names(which.max(max_change_sigs))
  summary$n_sigs_greater_03 <- toString(names(mean_sigs[mean_sigs > 0.3]))

  if (!is.null(mixtures.rescaled))
  {
    mean_sigs.rescaled <- apply(mixtures.rescaled, 1, mean)
    max_change_sigs.rescaled <- apply(mixtures.rescaled,1,max) - apply(mixtures.rescaled,1,min)
    sd_sigs.rescaled <- apply(mixtures.rescaled, 1, sd)

    summary$top_signature_rescaled <- names(which.max(mean_sigs.rescaled))
    summary$max_change_rescaled <- max(max_change_sigs.rescaled)
    summary$sig_with_max_change_rescaled <- names(which.max(max_change_sigs.rescaled))
    summary$n_sigs_greater_03_rescaled <- toString(names(mean_sigs.rescaled[mean_sigs.rescaled > 0.3]))
  }

  tmp <- toHorizontalMatrix(mean_sigs)
  colnames(tmp) <- paste0("mean_sigs.",names(mean_sigs))
  summary <- cbind(summary, tmp)

  tmp <- toHorizontalMatrix(max_change_sigs)
  colnames(tmp) <- paste0("max_change_sigs.",names(max_change_sigs))
  summary <- cbind(summary, tmp)

  tmp <- toHorizontalMatrix(max_change_sigs)
  colnames(tmp) <- paste0("sd_sigs.",names(max_change_sigs))
  summary <- cbind(summary, tmp)

  if (!is.null(mixtures.rescaled))
  {
    tmp <- toHorizontalMatrix(mean_sigs.rescaled)
    colnames(tmp) <- paste0("mean_sigs.rescaled.",names(mean_sigs.rescaled))
    summary <- cbind(summary, tmp)

    tmp <- toHorizontalMatrix(max_change_sigs.rescaled)
    colnames(tmp) <- paste0("max_change_sigs.rescaled.",names(max_change_sigs.rescaled))
    summary <- cbind(summary, tmp)

    tmp <- toHorizontalMatrix(max_change_sigs.rescaled)
    colnames(tmp) <- paste0("sd_sigs.rescaled.",names(max_change_sigs.rescaled))
    summary <- cbind(summary, tmp)
  }

  write.csv(summary, file=paste0(dir_name, "summary.csv"), row.names=F)
}

#' \code{extract_data_for_example} <man content>
#' @rdname helper_functions
extract_data_for_example <- function (example, dir_counts, tumortypes, dir_results = TrackSig.options()$DIR_RESULTS, dir_create = T) {
  #example <- "1c3df485-8e75-4378-87f6-c6463a520624"

  vcfFile <- paste0(dir_counts, "/", example, ".phi.txt")
  quadPhisFile <- paste0(dir_counts, "/", example, ".quadraticp.txt")

  vcfData <- tryCatch(read.table(vcfFile), error=function(e) NULL) # 96 trinucleotide counts are read as input
  quadratic_phis <- tryCatch(read.table(quadPhisFile)$V1, error=function(e) NULL)

  tumor_id <- gsub("^(.*)\\..*", "\\1", example)

  acronym <- toString(tumortypes[tumortypes$ID == tumor_id,]$tumor_type)
  dir_name <- paste0(dir_results, acronym, "/", tumor_id, "/")

  if (dir_create) {
    suppressWarnings(dir.create(dir_name, recursive = T))
  }

  if (is.null(vcfData))
  {
    return(list(tumor_id, vcfData, NULL, NULL, NULL, acronym, dir_name))
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

  vcfData[,1] <- NULL # Filenames on the first column are deleted

  if (sum(vcfData[,1] %% 1  != 0) > 0)
  {
    # second column represents phi values
    phis <- vcfData[,1]
    vcfData[,1] <- NULL
  }

  # Hack because of the previous bug in the code that assigned 101 mutations to each time points
  rows_keep <- ( apply(vcfData, 1, sum) == (TrackSig.options()$bin_size + 1) ) | ( apply(vcfData, 1, sum) == TrackSig.options()$bin_size )
  vcfData <- vcfData[rows_keep, ]

  if (!is.null(phis))
  {
    phis <- phis[rows_keep]
  }

  if(!is.null(assigns_phylo_nodes))
  {
    assigns_phylo_nodes <- assigns_phylo_nodes[rows_keep,1]
  }

  stopifnot(length(phis) == nrow(vcfData))

  #hack!! check how hundred are computed in perl script -- there is off by one error
  if (length(assigns_phylo_nodes) == (length(phis) + 1))
  {
    assigns_phylo_nodes <- assigns_phylo_nodes[-length(assigns_phylo_nodes)]
  }
  if(length(phis) != length(assigns_phylo_nodes))
  {
    assigns_phylo_nodes = NULL
  }

  return(list(tumor_id, vcfData, phis, quadratic_phis, assigns_phylo_nodes, acronym, dir_name))
}

#' \code{extract_data_for_simulation} <man content>
#' @rdname helper_functions
extract_data_for_simulation <- function (example, dir_counts, dir_results = DIR_RESULTS, dir_create = T) {
  #example <- "1c3df485-8e75-4378-87f6-c6463a520624"

  vcfFile <- paste0(dir_counts, "/", example, ".csv")

  vcfData <- tryCatch(read.csv(vcfFile, header=F), error=function(e) NULL) # 96 trinucleotide counts are read as input

  acronym <-  gsub("^(.*)\\/.*", "\\1", example)

  tumor_id <-  gsub("^.*\\/(.*)", "\\1", example)
  tumor_id <- gsub(" ", "_", tumor_id)

  dir_name <- paste0(dir_results, acronym, "/", tumor_id, "/")

  if (dir_create) {
    suppressWarnings(dir.create(dir_name, recursive = T))
  }

  if (is.null(vcfData))
  {
    return(list(tumor_id, vcfData, NULL, NULL, acronym, dir_name))
  }

  vcfData <- t(vcfData[,-c(1,2)])

  # Hack because of the previous bug in the code that assigned 101 mutations to each time points
  phis = 1:nrow(vcfData)

  stopifnot(length(phis) == nrow(vcfData))

  return(list(tumor_id, vcfData, phis, acronym, dir_name))
}

#' \code{extract_bootstrap_data_for_example} <man content>
#' @rdname helper_functions
extract_bootstrap_data_for_example <- function (example, bootstrap_counts) {

  vcfData.bootstrap <- list()
  vcfData.phis <- list()
  vcfData.bootstrap.unsorted <- list()

  for (vcfFile in (list.files(paste0(bootstrap_counts, "/", example, "/"))))
  {
    if (grepl(".*\\.(\\d+)\\.mut_counts.txt", vcfFile)) {
      index <- as.numeric(gsub(".*\\.(\\d+)\\.mut_counts.txt", "\\1", vcfFile))
      vcfFile <- paste0(bootstrap_counts, "/", example, "/", vcfFile)

      vcfData <- tryCatch(read.table(vcfFile), error=function(e) NULL) # 96 trinucleotide counts are read as input
      vcfData <- vcfData[,-c(1)] # Filenames on the first column are deleted

    if (round(vcfData[1,1]) != vcfData[1,1])
    {
      # second column represents phi values
      phis <- vcfData[,1]
      vcfData.phis[[index]] <- phis
      vcfData[,1] <- NULL
    }

       # Hack because of the previous bug in the code that assigned 101 mutations to each time points
       rows_keep <- (apply(vcfData, 1, sum) == 101) | (apply(vcfData, 1, sum) == 100)

       vcfData <- vcfData[rows_keep, ]

       vcfData.bootstrap[[index]] <- vcfData
     }
  }

  return(list(vcfData.bootstrap, vcfData.phis))
}

#' \code{gather_summaries_per_tissue} <man content>
#' @rdname helper_functions
gather_summaries_per_tissue <- function(omit_signature_information = F){
  for (dir_ in list.files(DIR_RESULTS, full.names= T))
  {
    if (!file.exists(paste0(dir_, "/")))
    {
      next
    }
    overall_summary <- c()
    for (file_ in list.files(dir_, full.names = T))
    {
      if (file.exists(paste0(file_, "/", "summary.csv")))
      {
        summary <- read.csv(paste0(file_, "/", "summary.csv"), stringsAsFactors=F)
        if (omit_signature_information) {
          summary <- summary[,1:(min(which(grepl("mean_sigs", colnames(summary))))-1)]
        }
        if (length(overall_summary) == 0) {
          overall_summary <- data.frame(summary, stringsAsFactors=F)
        } else {
          overall_summary <- rbind(overall_summary, summary)
        }
      }
    }
    write.csv(overall_summary, file=paste0(dir_, "/overall_summary.csv"), row.names=F)
  }
}

#' \code{add_noise} <man content>
#' @rdname helper_functions
add_noise <- function(vcfData, noise_rate = 0.05){
  sum_per_row <- apply(vcfData, 1, sum)
  n_noise <- round(sum_per_row * noise_rate)
  for (i in 1:nrow(vcfData))
  {
    noise_indices <- sample(ncol(vcfData), n_noise[i], replace=T)
    noise_indices <- table(noise_indices)
    noise_signs <- sample( c(1, -1), length(noise_indices), replace = T)
    vcfData[i,as.numeric(names(noise_indices))] <- vcfData[i,as.numeric(names(noise_indices))] + noise_signs * noise_indices
    vcfData[i,as.numeric(names(noise_indices))][vcfData[i,as.numeric(names(noise_indices))] < 0] <- 0
  }

  return(vcfData)
}

#' \code{convert_window_to_indices} <man content>
#' @rdname helper_functions
convert_window_to_indices <- function(total_length, window_size, gap){
  window_size = window_size/100
  gap = gap/100

  start_ends <- c()

  last_iteration = FALSE
  for (i in 1:total_length) {
    start = i * gap
    end = i*gap + window_size

    if (((i*gap + window_size) <= total_length) && ((i*gap + window_size + gap) > total_length)) {
      end = total_length
      last_iteration = TRUE
    }
    start_ends <- rbind(start_ends, data.frame(start = start, end=end))
    if (last_iteration) {
      break
    }
  }

  return(start_ends)
}

#' \code{find_indices_with_intersection} <man content>
#' @rdname helper_functions
find_indices_with_intersection <- function(covered_by_changepoints, intervals_for_change_histogram){
  converted_changepoints <- c()
  for (i in 1:nrow(covered_by_changepoints))
  {
    best_indice <- -1
    current_max <- 0
    for (j in 1:nrow(intervals_for_change_histogram))
    {
      int1 <- covered_by_changepoints[i,]
      int2 <- intervals_for_change_histogram[j,]

      intersection <- data.frame(start = max(int1$start, int2$start), end = min(int1$end, int2$end))
      intersection_max <- max(0, intersection$end - intersection$start)
      if (current_max < intersection_max)
      {
        best_indice <- j
        current_max <- intersection_max
      }
    }
    converted_changepoints <- c(converted_changepoints, best_indice)
  }

  return(converted_changepoints)
}

#' \code{read_mixtures} <man content>
#' @rdname helper_functions
read_mixtures <-  function(m_file) {
  m <- read.csv(m_file)
  rownames(m) <- m[,1]
  m <- m[,-1]
  colnames(m) <- as.numeric(gsub("X(.*)", "\\1", colnames(m)))
  return(m)
}

#' \code{get_examples_group} <man content>
#' @rdname helper_functions
get_examples_group <- function(list_, EXAMPLES_PER_GROUP = TrackSig.options()$EXAMPLES_PER_GROUP, group = TrackSig.options()$group){
  if (group == 0)
  {
    examples_group <- list_
  } else {
    start = (EXAMPLES_PER_GROUP*(group-1)) + 1
    end = min(EXAMPLES_PER_GROUP*(group), length(list_))
    if (end >= start)
    {
      indices <- start:end
      examples_group <- list_[indices]
    } else {
      print("No examples in a group")
      examples_group <- c()
    }
  }
}


#' \code{compute_mean_sd_err} <man content>
#' @rdname helper_functions
compute_mean_sd_err <- function(mixtures_bootstrap, sig_names, dir_name = NULL, descr = "", exclude_zero=F){
  mixtures.sd <- data.frame(matrix(NA,ncol=ncol(mixtures_bootstrap[[1]]), nrow=nrow(mixtures_bootstrap[[1]])))
  mixtures.err <- data.frame(matrix(NA,ncol=ncol(mixtures_bootstrap[[1]]), nrow=nrow(mixtures_bootstrap[[1]])))
  mixtures.mean <- data.frame(matrix(NA,ncol=ncol(mixtures_bootstrap[[1]]), nrow=nrow(mixtures_bootstrap[[1]])))

  for (i in 1:nrow(mixtures_bootstrap[[1]]))
    for (j in 1:ncol(mixtures_bootstrap[[1]]))
    {
      values <- sapply(mixtures_bootstrap, function(x) x[i,j])
      if (exclude_zero) {
        values <- values[values != 0]
      }
      mixtures.sd[i,j] <- sd(unlist(values), na.rm=T)
      mixtures.err[i,j] <- sd(unlist(values), na.rm=T)/sqrt(length(values))
      mixtures.mean[i,j] <- mean(unlist(values), na.rm=T)
    }
  rownames(mixtures.mean) <- rownames(mixtures.sd) <- rownames(mixtures.err) <- sig_names


  avg_phis <- c()
 for (j in 1:ncol(mixtures_bootstrap[[1]]))
  {
    values <- sapply(mixtures_bootstrap, function(x) as.numeric(colnames(x))[j])
    avg_phis <- c(avg_phis, mean(values, na.rm=T))
  }
  colnames(mixtures.mean) <- colnames(mixtures.sd) <- colnames(mixtures.err) <- avg_phis

  if (!is.null(dir_name)) {
    write.csv(mixtures.mean, file=paste0(dir_name, "mixtures.mean", descr, ".csv"))
    write.csv(mixtures.sd, file=paste0(dir_name, "mixtures.sd", descr, ".csv"))
    write.csv(mixtures.err, file=paste0(dir_name, "mixtures.err", descr, ".csv"))
  }
  return(list(mixtures.mean, mixtures.sd, mixtures.err,avg_phis))
}

#' \code{get_bootstrap_mixtures} <man content>
#' @rdname helper_functions
get_bootstrap_mixtures <- function(bootstrap_vcfs, bootstrap_phis, alex.t, dir_name, descr, verbose=TRUE){
  mixtures_bootstrap <- list()
  changepoints_bootstrap <- list()
  for (j in 1:length(bootstrap_vcfs))
  {
    if (verbose) {
      print(paste0(dir_name, "changepoints.bootstrap_", j, descr, ".txt"))
    }
    if (!file.exists(paste0(dir_name,"mixtures.bootstrap_", j, descr, ".csv")) | !file.exists(paste0(dir_name, "changepoints.bootstrap_", j, descr, ".txt")))
    {
      if (TrackSig.options()$changepoint_method == "PELT") {
        list[cp, m] <- find_changepoints_pelt(t(bootstrap_vcfs[[j]]), alex.t)
      } else {
        list[bics, optimal, cp, m] <- find_changepoints_over_all_signatures_one_by_one(bootstrap_vcfs[[j]], alex.t, n_signatures = ncol(alex.t))
      }

      write.csv(m, file=paste0(dir_name, "mixtures.bootstrap_", j, descr, ".csv"))
      n_col <- ifelse(length(cp) > 0, length(cp), 1)
      write(cp, file=paste0(dir_name, "changepoints.bootstrap_", j, descr, ".txt"), ncolumns=n_col)
    } else {
      m <- read_mixtures(paste0(dir_name,"mixtures.bootstrap_", j, descr, ".csv"))
      cp_file = paste0(dir_name, "changepoints.bootstrap_", j, descr, ".txt")
      if (file.info(cp_file)$size == 1) {
        cp <- c()
      } else {
        cp <- unlist(read.table(paste0(dir_name, "changepoints.bootstrap_", j, descr, ".txt"), header=F))
      }
    }

    colnames(m) <- bootstrap_phis[[j]][1:ncol(m)]
    mixtures_bootstrap[[j]] <- m
    changepoints_bootstrap[[j]] <- cp
  }

  return(list(mixtures_bootstrap, changepoints_bootstrap))
}

#' \code{add_list_by_element} <man content>
#' @rdname helper_functions
add_list_by_element <- function(list_total_sum, list_to_add) {
  for (name in names(list_to_add)) {
    list_total_sum[[name]] <- c(list_total_sum[[name]], list_to_add[[name]])
  }
  return(list_total_sum)
}

#' \code{not_na} <man content>
#' @rdname helper_functions
not_na <- function(x) {
  return(x[!is.na(x)])
}

#' \code{truncate_to_range} <man content>
#' @rdname helper_functions
truncate_to_range <- function(mixtures, range_) {
  min = range_[1]
  max = range_[2]

  x <- mixtures
  col_names <- as.numeric(colnames(x))
  to_leave <- which(col_names <= max+0.01 & col_names >= min-0.01)

  x2 <- x[,to_leave, drop=F]
  colnames(x2) <- col_names[to_leave]
  return(list(x2,to_leave))
}


#' \code{get_sample_purity} <man content>
#' @rdname helper_functions
get_sample_purity <- function(tumor_id) {
  purities <-read.delim(TrackSig.options()$purity_file)

  sample_purity = NULL
  if (!(tumor_id %in% purities$samplename)) {
    warning(paste("Tumor not found", tumor_id))
  } else {
    sample_purity = purities[purities$samplename == tumor_id,]$purity
  }

  return(sample_purity)
}

#' \code{merge_signatures} <man content>
#' @rdname helper_functions
merge_signatures <- function(mixtures, sigs_to_merge) {
  if (!is.null(sigs_to_merge)) {
    for (i in 1:length(sigs_to_merge)) {
      set_name <- names(sigs_to_merge)[i]
      sig_set = sigs_to_merge[[i]]
      sig_set = intersect(rownames(mixtures), sig_set)

      if (length(sig_set) == 0) {
        next
      }

      new_mixture_col <- toHorizontalMatrix(apply(mixtures[sig_set,,drop=F],2,sum))
      rownames(new_mixture_col) <- set_name
      colnames(new_mixture_col) <- colnames(mixtures)

      if (length(intersect(rownames(mixtures), sig_set)) == nrow(mixtures)) {
        mixtures <- new_mixture_col
      } else {
        mixtures <- rbind(mixtures[-which(rownames(mixtures) %in% sig_set),,drop=FALSE], new_mixture_col)
      }
    }
  }
  return(mixtures)
}

#' \code{extract_exposures_per_mutation} <man content>
#' @rdname helper_functions
#' @export

extract_exposures_per_mutation <- function(activities_dir, sorted_mutations_dir, 
  bin_size = 100, samples_to_run  = c()) {
  # activities_dir: path to the
  # sorted_mutations_dir: folder with files for each tumour sample (or simulations). Each file has a list of mutations SORTED BY CCF
  tumor_list <- list.dirs(activities_dir, recursive = F, full.names=F)

  if (length(samples_to_run) > 0) {
  	# run only these tumour samples
	tumor_list <- intersect(tumor_list, samples_to_run)
  }

  print("Extracting exposures per mutation...")
  for (tumor in tumor_list) {
    print(tumor)
    
    sorted_mut_file = paste0(sorted_mutations_dir, "/", tumor, ".mut_types.txt")
    if (!file.exists(sorted_mut_file)) {
      print(sprintf("File %s does not exist", sorted_mut_file))
      next
    }

    mut_list <- read.delim(sorted_mut_file, header=F, stringsAsFactors=F)
    colnames(mut_list) <- c("chr", "pos", "ccf", "ref", "alt", "tri")
    ccfs <- mut_list[,"ccf"]

    activity_file <- paste0(activities_dir, "/", tumor, "/mixtures.csv")
    if (!file.exists(activity_file)) {
      print(sprintf("Activity file %s does not exist", activity_file))
      next
    }

    activities <- read.csv(activity_file, header=T, stringsAsFactors=F)
    rownames(activities) <- sapply(activities[,1], toString)
    activities <- activities[,-1]
    activities <- t(activities)

    n_time_points = nrow(activities)

    #stopifnot(nrow(mut_list) <= n_time_points * bin_size)

    activities_per_mut <- c()
    for (i in 1:n_time_points) {
      mut_current_tp = mut_list[(bin_size*(i-1)+1):(bin_size * i), ]
      stopifnot(nrow(mut_current_tp) == bin_size)

      current_activities = activities[i,,drop=F]
      # repeat current_activities for every mutation
      current_activities <- current_activities[rep(1, nrow(mut_current_tp)),]

      data <- suppressWarnings(data.frame(chromosome=mut_current_tp$chr, start=mut_current_tp$pos, current_activities))
      activities_per_mut <- rbind(activities_per_mut, data)
    }

    stopifnot(nrow(activities_per_mut) == nrow(mut_list))
    write.table(activities_per_mut, file = paste0(activities_dir, "/", tumor, "/sig_exposures_per_mut.txt"), 
        sep = "\t", row.names=F, quote=F)
  }
}

#' \code{otation} <man content>
#' @rdname helper_functions
#' @export
load_annotation <- function(tumortype_file = TrackSig.options()$tumortype_file, signature_file = TrackSig.options()$signature_file,
                            active_signatures_file = TrackSig.options()$active_signatures_file){

  names_trinucleotide <- trinucleotide_internal
  names_trinucleotide <- apply(names_trinucleotide, 1, function(x) { do.call("paste", c(as.list(x), sep = "_"))})

  # Load the tumor types for tumor IDs.
  tumortypes <- read.delim(tumortype_file, header = T, stringsAsFactors=F)
  colnames(tumortypes) <- c("ID", "tumor_type")

  # ALEX DATA
  # The trinucleotide count matrix will be regressed on the 30 Alexandrov Mutational Signature frequencies
  # http://cancer.sanger.ac.uk/cosmic/signatures
  alex <- read.table(signature_file)
  rownames(alex) <- names_trinucleotide
  colnames(alex) <- paste0("S", 1:ncol(alex))

  if (TrackSig.options()$cancer_type_signatures) {
    # Load active signatures for each tumor type
    active_signatures <- read.delim(TrackSig.options()$active_signatures_file, stringsAsFactors=F)
    active_signatures[is.na(active_signatures)] <- 0
    # Removing column "Other.signatures"
    active_signatures <- active_signatures[, -ncol(active_signatures)]

    # Some of the signatures have several types corresponding to them.
    # The rows which correspond to several tumor types and duplicated in the table.
    indices_to_duplicate <- which(grepl(",",active_signatures$acronym))
    while (length(indices_to_duplicate) != 0)
    {
      index <- indices_to_duplicate[1]
      current_type <- active_signatures[index,]
      trim <- function (x) gsub("^\\s+|\\s+$", "", x)
      corresponding_tumor_types <- trim(unlist(strsplit(current_type$acronym, ",")))

      to_add <- current_type
      for (i in 2:length(corresponding_tumor_types))
      {
        to_add <- rbind(to_add, current_type)
      }
      to_add$acronym <- corresponding_tumor_types

      active_signatures <- rbind(active_signatures[1:(index-1),], to_add, active_signatures[(index+1):nrow(active_signatures),])
      indices_to_duplicate <- which(grepl(",",active_signatures$acronym))
    }
    active_signatures.our_samples <- merge(tumortypes, active_signatures, by.x="tumor_type", by.y="acronym", all.x=T)
    # active_signatures.lydia_samples <- merge(tumortypes.lydia, active_signatures, by.x="tumor_type", by.y="acronym", all.x=T)
  } else {
    active_signatures <- NULL

    active_signatures.our_samples <- read.table(TrackSig.options()$active_signatures_file, stringsAsFactors=F, header = F)

    #active_signatures.our_samples <- read.csv("PCAWG_signatures_in_samples_beta.csv")
    colnames(active_signatures.our_samples) <- c("tumor_type",  "ID", colnames(active_signatures.our_samples)[3:ncol(active_signatures.our_samples)])

    # Making the table binary
    active_signatures.our_samples.data <- active_signatures.our_samples[,3:ncol(active_signatures.our_samples)]
    active_signatures.our_samples.data[active_signatures.our_samples.data > 0] <- 1
    colnames(active_signatures.our_samples.data) <- gsub("S.(.*)", "\\1", colnames(active_signatures.our_samples.data))

    active_signatures.our_samples <- cbind(active_signatures.our_samples[,c(1,2)], Name=NA, active_signatures.our_samples.data)
    active_signatures.our_samples$Name <- sapply(active_signatures.our_samples$tumor_type, toString)
  }

  return(list(alex, tumortypes, active_signatures, active_signatures.our_samples))
}


#' \code{load_annotation_pcawg} <man content>
#' @rdname helper_functions
load_annotation_pcawg <- function(tumortype_file = TrackSig.options()$tumortype_file, signature_file =TrackSig.options()$signature_file,
                                  active_signatures_file = TrackSig.options()$active_signatures_file) {

  names_trinucleotide <- trinucleotide_internal
  names_trinucleotide <- apply(names_trinucleotide, 1, function(x) { do.call("paste", c(as.list(x), sep = "_"))})

  # Load the tumor types for tumor IDs.
  tumortypes <- read.delim(tumortype_file, header = F, stringsAsFactors=F)
  colnames(tumortypes) <- c("ID", "tumor_type")

  # Signatures from PCAWG group
  alex <- read.table(signature_file, header = T, stringsAsFactors = FALSE)
  #alex <- read.csv(paste0(DIR, "PCAWG_signature_patterns_beta.csv"), header = T)
  #pcawg_trinucleotides <- paste0(gsub(">", "_", alex[,1]), "_", alex[,2])
  pcawg_trinucleotides <- paste(substr(alex[,1], 2, 2), substr(alex[,1], 5, 5), substr(alex[,1], 1, 3), sep="_")
  rownames(alex) <- pcawg_trinucleotides
  alex <- alex[,-1]
  # to make the order of mutation types match names_trinucleotide
  alex <- alex[match(names_trinucleotide, pcawg_trinucleotides),]
  #colnames(alex) <- gsub("Signature.(.*)", "\\1", colnames(alex))
  #colnames(alex) <- gsub("PCAWG.(.*)", "\\1", colnames(alex))


  sigs_to_merge <- list()
  sigs_to_merge[["SBS7"]] <- c("SBS7a", "SBS7b", "SBS7c", "SBS7d")
  sigs_to_merge[["SBS17"]] <- c("SBS17a", "SBS17b")
  sigs_to_merge[["SBS2+13"]] <- c("SBS2", "SBS13")
  sigs_to_merge[["SBS10"]] <- c("SBS10a", "SBS10b")

  alex_merged <- t(merge_signatures(t(alex), sigs_to_merge))

  active_signatures <- active_signatures.our_samples <- NULL
  if (!is.null(active_signatures_file)) {
    active_signatures.our_samples <- read.delim(active_signatures_file, stringsAsFactors=F)
    colnames(active_signatures.our_samples) <- c("tumor_type",  "ID", colnames(active_signatures.our_samples)[3:ncol(active_signatures.our_samples)])

    active_signatures.our_samples.data <- active_signatures.our_samples[,3:ncol(active_signatures.our_samples)]
    active_signatures.our_samples.data[active_signatures.our_samples.data > 0] <- 1

    active_signatures.our_samples <- cbind(active_signatures.our_samples[,c(1,2)], Name=NA, active_signatures.our_samples.data)
    active_signatures.our_samples$Name <- sapply(active_signatures.our_samples$tumor_type, toString)

    if ("SBS2.13" %in% colnames(active_signatures.our_samples)) {
      colnames(active_signatures.our_samples)[colnames(active_signatures.our_samples) == "SBS2.13"] = "SBS2+13"
    }
  }

  return(list(alex_merged, tumortypes, active_signatures, active_signatures.our_samples))
}

#[END]
