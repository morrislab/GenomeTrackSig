# AUTHOR: Yulia Rubanova
# Modified for package TrackSig by Cait Harrigan



#' \code{generateContext} Generate a trinucleotide context from an alphabet. Note: this involves finding all three-member
#' permutations of the alphabet, which can be inconveinent for large alphabets. Nucleotides are assumed to be provided as complementary pairs,
#' where the first of each pair is used as the reference to build the context.
#'
#' @param alphabet list of pairs of characters to create combinations of as a mutation context type
#' @return data.frame containing all the possible trinucleotide contextes for a mutation in the supplied alphabet
#'
#' @examples
#' context <- TrackSig:::generateContext(c("CG", "TA"))
#' dim(context)
#' head(context)
#'
#' @rdname helper_functions
#' @name generateContext
#' @export

generateContext <- function(alphabet){

  if (any(nchar(alphabet) != 2)){
    stop("Alphabet is malformed. Please provide alphabet as a list of complementary pairs")
  }

  allpha <- unlist(strsplit(alphabet, split=NULL))
  nTypes <- (length(allpha) - 1) * length(allpha)^3 * 1/2

  context <- data.frame()

  for (i in seq(1, length(allpha), by = 2)){

    midRef <- allpha[i]
    rest <- setdiff(allpha, midRef)
    repSize <- length(allpha)^2 - length(allpha)

    midSet <- cbind(rep(midRef, length.out = repSize), rep(rest, length.out=repSize),
                    paste0(sort(rep(allpha, repSize)), rep(midRef, length.out = repSize), rep(allpha, repSize)))
    context <- rbind(context, midSet)
  }

  stopifnot( dim(context)[1] == nTypes )

  return (context)
}

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
  warning("Called a depricated function.")
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
  warning("Called a depricated function.")
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

IgnoreVectorOrMatrix <- function(x, FUN)
{
  warning("Called a depricated function.")
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

  warning("get_sliding_window_data() is depricated.")

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


#' \code{merge_data_chunks} <man content>
#' @rdname helper_functions
merge_data_chunks <- function (vcfData) {
  warning("Called a depricated function.")
  vcf <- matrix(nrow = nrow(vcfData)/2, ncol = ncol(vcfData))
  for (i in 1:nrow(vcf)) {
    for (j in 1:ncol(vcf)) {
      vcf[i,j] <- vcfData[i*2-1,j] + vcfData[i*2,j]
    }
  }
  vcf <- t(vcf)
  return(vcf)
}

#' \code{get_noise_sig} <man content>
#' @rdname helper_functions
get_noise_sig <- function(noise_sig){
  warning("Called a depricated function.")
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
  warning("Called a depricated function.")
  print("Helper function called! get_signatures_for_current_sample()")

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


#' \code{toHorizontalMatrix} <man content>
#' @rdname helper_functions
toHorizontalMatrix <- function(L){
  warning("Called a depricated function.")
  if (is.vector(L))
    return(matrix(L, nrow=1))
  else
    return(as.matrix(L))
}


#' \code{extract_data_for_example} <man content>
#' @rdname helper_functions
extract_data_for_example <- function (example, dir_counts, tumortypes, dir_results = TrackSig.options()$DIR_RESULTS, dir_create = T) {
  #example <- "1c3df485-8e75-4378-87f6-c6463a520624"
  warning("Called a depricated function.")
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






#' \code{truncate_to_range} <man content>
#' @rdname helper_functions
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


#' \code{get_sample_purity} <man content>
#' @rdname helper_functions
get_sample_purity <- function(tumor_id) {

  warning("get_sample_purity called. phi scaling should be done when corrected VAF made")
  warning("Called a depricated function.")
  purities <-read.delim(TrackSig.options()$purity_file)

  sample_purity = NULL
  if (!(tumor_id %in% purities$samplename)) {
    warning(paste("Tumor not found", tumor_id))
  } else {
    sample_purity = purities[purities$samplename == tumor_id,]$purity
  }

  return(sample_purity)
}





#[END]
