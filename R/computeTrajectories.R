# computeTrajectories.R
# Authors: Yulia Rubanova, Cait Harrigan

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


make_binary_table <- function(multinomial_vector)
{
  N = sum(multinomial_vector)
  mutation_types <- length(multinomial_vector)

  # Make a matrix of samples for fitting  of multinomials.
  # Each sample contains the one mutation.
  # data.matrix is a binary matrix with N columns and mutation_types rows.
  data.matrix <- matrix(0, ncol=N, nrow=mutation_types)
  current_index <- 1
  for (i in 1:length(multinomial_vector))
  {
    if (multinomial_vector[i] != 0)
    {
      data.matrix[i, current_index:min(current_index+multinomial_vector[i]-1, N) ] <- 1
    }
    current_index <- current_index + multinomial_vector[i]
  }

  if (current_index != N + 1)
    stop("Something went wrong during binary matrix construction: current_index is off")

  return(data.matrix)
}

# Fit mixture of multinomials for each column of the matrix
fit_mixture_of_multinomials_matrix <- function (vcf, alex.t, prior=NULL)
{
  dd <- matrix(nrow = ncol(alex.t), ncol = 0)
  for (m in 1:ncol(vcf)) {
    mixtures <- fit_mixture_of_multinomials_EM(vcf[,m], as.matrix(alex.t))
    dd <- cbind(dd, as.array(mixtures))
  }

  rownames(dd) <- colnames(alex.t)

  return(dd)
}

# fit mixture of multinomials to the vector
fit_mixture_of_multinomials_EM <- function(multinomial_vector, composing_multinomials, prior=NULL)
{
  # Number of mutations to fit
  N = sum(multinomial_vector)
  # Number of mutation types / categories of mutinomial
  mutation_types <- length(multinomial_vector)
  # Number of multinomials/signatures to fit and to make mixture of
  M <- ncol(composing_multinomials)

  if (length(multinomial_vector) != nrow(composing_multinomials))
  {
    print(length(multinomial_vector))
    print(dim(composing_multinomials))
    stop("Length of data vector is not equal to nrow of matrix to fit. Did you forget to transpose the matrix?")
  }
  data.matrix <- make_binary_table(multinomial_vector)

  # data_given_class[i,n] corresponds to class/signature i and sample/mutation n
  data_given_class <- matrix(0, nrow=M, ncol=N)
  for (i in 1:M)
  {
    data_given_class[i,] <- apply(composing_multinomials[,i]^data.matrix,2,prod)
  }

  # Mixtures of multinomials. Use uniform prior unless the prior is specified
  if (!is.null(prior))
  {
    if (length(prior) != M)
      stop(paste0("Length of prior should be equal to ", M))
  } else {
    prior <- rep(1/M, M)
  }

  pi <- prior
  pi_diff <- Inf
  iteration <- 1

  while (pi_diff > 0.001 & iteration < 1000)
  {
    # E-step: update posterior.
    p_x <- apply(data_given_class * pi, 2, sum)

    # class_given_data[i,n] corresponds to class/signature i and sample/mutation n
    class_given_data <- t(t(data_given_class * pi) / p_x)

    # S-step: update mixtures
    pi_new <- 1/N * apply(class_given_data,1,sum)

    if (sum(pi_new > 1) != 0) {
      stop("Mixture ratio is greater than 1")
    }

    if (sum(pi_new < 0) != 0)
      stop("Mixture ratio is less than 0")

    if (sum(pi_new) > 1.5)
      stop("Sum of mixture ratios is greater than 1")

    pi_diff <- sum(abs(pi_new - pi))
    pi <- pi_new
    iteration <- iteration + 1
  }

  return(pi)
}

# fit mixture of mutinomials in each time slice specified by change_points
fit_mixture_of_multinomials_in_time_slices <- function(data, change_points, alex.t, split_data_at_change_point = T)
{
  fitted_values <- matrix(NA, ncol=ncol(data), nrow=ncol(alex.t))
  rownames(fitted_values) <- colnames(alex.t)

  # Get first time slice until the first check point and get sum of it
  if (length(change_points) == 0) {
    end_of_first_slice <- ncol(data)
    slice_indices <- 1:(end_of_first_slice)
  } else {
    end_of_first_slice <- change_points[1]
    slice_indices <- 1:(end_of_first_slice-1)
  }

  d <- list()

  d[[1]] <- list()
  d[[1]]$data <-  toVerticalMatrix(data[,slice_indices])
  d[[1]]$slice_indices <- slice_indices

  if (split_data_at_change_point)
  {
    right_side <- left_side <- floor(data[,end_of_first_slice] / 2)
    leftovers <- data[,end_of_first_slice] - right_side - left_side
    leftovers_right <- leftovers
    leftovers_right[ sort(sample(which(leftovers > 0), length(which(leftovers > 0)) / 2))] <- 0
    leftovers_left <- leftovers - leftovers_right

    right_side <- right_side + leftovers_right
    left_side <- left_side + leftovers_left

    stopifnot(sum(data[,end_of_first_slice]) == sum(right_side) + sum(left_side))

    d[[1]]$data <- cbind(d[[1]]$data, left_side)
    left_side <- c()
  }

  if (length(change_points) > 1)
  {
    for (i in 2:length(change_points))
    {
      slice_indices <- (change_points[i-1]):(change_points[i]-1)
      d[[i]] <- list()
      d[[i]]$data <- toVerticalMatrix(data[,slice_indices])
      d[[i]]$slice_indices <- slice_indices

      if (split_data_at_change_point)
      {
        d[[i]]$data <- toVerticalMatrix(cbind(right_side, d[[i]]$data[,-1]))
        right_side <- c()

        right_side <- left_side <- floor(data[,slice_indices[length(slice_indices)]] / 2)
        leftovers <- data[,slice_indices[length(slice_indices)]] - right_side - left_side
        leftovers_right <- leftovers
        leftovers_right[ sort(sample(which(leftovers > 0), length(which(leftovers > 0)) / 2))] <- 0
        leftovers_left <- leftovers - leftovers_right

        right_side <- right_side + leftovers_right
        left_side <- left_side + leftovers_left

        stopifnot(sum(data[,slice_indices[length(slice_indices)]]) == sum(right_side) + sum(left_side))

        d[[i]]$data <- cbind(d[[i]]$data, left_side)
        left_side <- c()
      }
    }
  }

  if (length(change_points) != 0)
  {
    slice_indices <- (change_points[length(change_points)]):ncol(data)
    d[[length(d) + 1]] <- list()
    d[[length(d)]]$data <- toVerticalMatrix(data[,slice_indices])
    d[[length(d)]]$slice_indices <- slice_indices

    if (split_data_at_change_point)
    {
      d[[length(d)]]$data <- cbind(right_side, toVerticalMatrix(d[[length(d)]]$data)[,-1])
      right_side <- c()
    }
  }


  for (i in 1:length(d))
  {
    current_d <- d[[i]]
    current_d.sum <- IgnoreVectorOrMatrix(current_d$data, FUN=function(x) {apply(x,1, sum)})
    fitted_for_time_slice <- fit_mixture_of_multinomials_EM(current_d.sum, alex.t)

    fitted_values[,current_d$slice_indices] <- matrix(rep(fitted_for_time_slice, length(current_d$slice_indices)), nrow=nrow(fitted_values))
  }
  colnames(fitted_values) <- colnames(data)
  return(fitted_values)
}



sig_mixture_ll <- function(multinomial_vector, composing_multinomials, mixtures, ...) {
  # replaces log_likelihood_mixture_multinomials
  mutation_binary_table <-  make_binary_table(multinomial_vector)

  # mutation_probabilities_under_multinomial[i,n] corresponds to class/signature i and sample/mutation n
  mutation_probabilities_under_multinomial <- matrix(0, nrow=ncol(composing_multinomials), ncol=ncol(mutation_binary_table))
  for (sig in 1:ncol(composing_multinomials)) {
    mutation_probabilities_under_multinomial[sig,] <- apply(composing_multinomials[,sig]^mutation_binary_table,2,prod)
  }

  mutation_probabilities_under_mixture <-  log(t(mutation_probabilities_under_multinomial) %*% as.matrix(mixtures))
  stopifnot(length(mutation_probabilities_under_mixture) == sum(multinomial_vector))

  return(sum(mutation_probabilities_under_mixture))
}

# beta likelihood maximization
beta_ll <- function(qis, ...){

  #qis are the VAFs for the subproblem

  n <- length(qis)

  #assertthat::assert_that(length(qis) == length(vis), length(qis) == length(ris), msg = "problem subsetting is not good!")

  alpha <- sum(qis) + 1
  beta <- sum(1-qis) + 1

  LL <- lbeta(alpha, beta) + log(pbeta(max(qis), alpha, beta) - pbeta(min(qis), alpha, beta))

  #print(c(alpha, beta, LL))

  return(LL)

}

sum_beta_mixture_ll <- function(qis, multinomial_vector,
                                composing_multinomials, mixtures, ...){

  return( sum(sig_mixture_ll(multinomial_vector, composing_multinomials, mixtures),
              beta_ll(qis)) )

}

parseScoreMethod <- function(scoreMethod){
  # return the penalty and score function to use when computing partitions

  assertthat::assert_that(scoreMethod %in% c("TrackSig", "TrackSigFreq"), msg = "scoreMethod should be one of \"TrackSig\", \"TrackSigFreq\".
                                                                                 Please see documentation for more information on selecting a scoreMethod)")
  if(scoreMethod == "TrackSig"){
    return(list(penalty = expression((n_sigs - 1) * log(n_bins * binSize)),
                score_fxn = sig_mixture_ll))
  }

  if(scoreMethod == "TrackSigFreq"){
    return(list(penalty = expression(-log(0.1) + (n_sigs + 1) * log(n_bins * binSize)),
                score_fxn = sum_beta_mixture_ll))
  }

}

getActualMinSegLen <- function(desiredMinSegLen, binSize){
  # return the minimum segment length to use.

  # for best segment scoring, use at least 400 mutations per segment.
  if(is.null(desiredMinSegLen)){
    return (ceiling(400/binSize))
  }

  # for accurate segment scoring, reqire at least 100 mutations per segment.
  actualMinSegLen <- max(desiredMinSegLen, ceiling(100/binSize))

  if (actualMinSegLen != desiredMinSegLen){
    warning(sprintf("Could not use desiredMinSegLen, too few mutations for accurate segment scoring. minSegLen set to: %s", actualMinSegLen))
  }

  return(actualMinSegLen)
}

# Find optimal changepoint and mixtures using PELT method.
# if desiredMinSegLen is NULL, the value will be selected by default based off binSize to try to give good performance
find_changepoints_pelt <- function(countsPerBin, alex.t, vcaf, scoreMethod = "TrackSigFreq", binSize = 100, desiredMinSegLen = NULL)
{

  minSegLen <- getActualMinSegLen(desiredMinSegLen, binSize)

  score_matrix <- score_partitions_pelt(countsPerBin, alex.t, vcaf, scoreMethod, binSize, minSegLen)
  changepoints <- recover_changepoints(score_matrix)

  mixtures <- fit_mixture_of_multinomials_in_time_slices(countsPerBin, changepoints, alex.t)

  return(list(changepoints = changepoints, mixtures = mixtures))
}

# Calculate penalized BIC score for all partitions using PELT method.
score_partitions_pelt <- function(countsPerBin, alex.t, vcaf, scoreMethod, binSize, minSegLen)
{
  n_bins <- dim(countsPerBin)[2]
  n_sigs <- dim(alex.t)[2]

  list[penalty, score_fxn] <- parseScoreMethod(scoreMethod)
  penalty <- eval(penalty)

  # aggregate bin summary stats
  phis <- aggregate(vcaf$phi, by = list(vcaf$binAssignment), FUN = sum)$x
  quadratic_phis <- aggregate(vcaf$phi, by = list(vcaf$binAssignment), FUN = function(x){return(sum(x^2))})$x

  # force vaf permutation
  if(TRUE){
    phis <- aggregate(vcaf$phi2, by = list(vcaf$binAssignment), FUN = sum)$x
    quadratic_phis <- aggregate(vcaf$phi2, by = list(vcaf$binAssignment), FUN = function(x){return(sum(x^2))})$x

  }

  # Bayeisan Information Criterion penalization constant defalut parameter

  # Store score for all partitions of all sub-problems
  # Rows are length of sub-problem. Columns correspond to last changepoint
  sp_scores <- matrix(nrow=n_bins, ncol=n_bins)

  max_sp_scores <- numeric(n_bins)
  prune_set <- c()

  # Score all subproblems of length sp_len using last_cp as last changepoint
  for (sp_len in 1:n_bins)
  {
    valid_cps <- setdiff(0:(sp_len - 1), prune_set)
    print(paste0("Scoring subpartitions of length: ", sp_len, "/", n_bins))
    for (last_cp in valid_cps)
    {
      # Segments with length less than 4 cannot be accurately scored
      if (sp_len - last_cp < 4)
      {
        sp_scores[sp_len, last_cp + 1] <- -Inf
        next
      }

      sp_slice <- c((last_cp + 1), sp_len)

      r_seg_phis <- phis[sp_slice[1] : sp_slice[2]]
      r_seg_quadratic_phis <- quadratic_phis[sp_slice[1] : sp_slice[2]]

      r_seg_qis <- vcaf$phi[vcaf$binAssignment %in% (sp_slice[1] : sp_slice[2])]
      r_seg_q_cni <- vcaf$cn[vcaf$binAssignment %in% (sp_slice[1] : sp_slice[2])]

      # force vaf permutation
      if(TRUE){
        r_seg_qis <- vcaf$phi2[vcaf$binAssignment %in% (sp_slice[1] : sp_slice[2])]
      }

      r_seg_qis <- r_seg_qis / (2 + vcaf$purity[1] * (r_seg_q_cni - 2))
      r_seg_qis <- unlist(lapply(r_seg_qis, 1, FUN = min))

      r_seg_vi <- vcaf$vi[vcaf$binAssignment %in% (sp_slice[1] : sp_slice[2])]
      r_seg_ri <- vcaf$ri[vcaf$binAssignment %in% (sp_slice[1] : sp_slice[2])]

      r_seg_counts <- rowSums(countsPerBin[, sp_slice[1] : sp_slice[2], drop = FALSE])

      r_seg_mix <- fit_mixture_of_multinomials_EM(r_seg_counts, alex.t)


      r_seg_score <- 2 * score_fxn(multinomial_vector = r_seg_counts, phis = r_seg_phis, quad_phis = r_seg_quadratic_phis,
                                   composing_multinomials = alex.t, mixtures = r_seg_mix, bin_size = bin_size, qis = r_seg_qis,
                                   vis = r_seg_vi, ris = r_seg_ri)

      l_seg_score <- ifelse(last_cp == 0, penalty, max_sp_scores[last_cp])

      sp_scores[sp_len, last_cp + 1] <- l_seg_score + r_seg_score - penalty
    }

    max_sp_scores[sp_len] <- max(sp_scores[sp_len, ][!is.na(sp_scores[sp_len, ])])

    # Evaluate all changepoints for pruning condition
    for (cp in valid_cps)
    {
      if (sp_len - cp < 4)
      {
        next
      }

      if (sp_scores[sp_len, cp + 1] + penalty < max_sp_scores[sp_len])
      {
        prune_set <- c(prune_set, cp)
      }
    }
  }
  return(sp_scores)
}

# Recover optimal changepoints by from subproblem matrix
recover_changepoints <- function(sp_score_matrix)
{
  changepoints <- c()

  continue <- TRUE
  current <- dim(sp_score_matrix)[1]
  while (continue)
  {
    prev <- which.max(sp_score_matrix[current, ])
    if (prev - 1 <= 1)
    {
      continue <- FALSE
    }
    else
    {
      changepoints <- c(prev - 1, changepoints)
    }
    current <- prev - 1
  }
  return(changepoints)
}


# [END]
