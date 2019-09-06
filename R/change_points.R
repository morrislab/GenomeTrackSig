# AUTHOR: Yulia Rubanova
# modified for package TrackSig by Cait Harrigan


# Find optimal changepoint and mixtures using PELT method.
find_changepoints_pelt <- function(vcf, alex.t, vcaf)
{
  score_matrix <- score_partitions_pelt(vcf, alex.t, vcaf)
  changepoints <- recover_changepoints(score_matrix)

  mixtures <- fit_mixture_of_multinomials_in_time_slices(vcf, changepoints, alex.t)

  return(list(changepoints = changepoints, mixtures = mixtures))
}

# Calculate penalized BIC score for all partitions using PELT method.
score_partitions_pelt <- function(vcf, alex.t, vcaf,
                                  penalty = TrackSig.options()$pelt_penalty,
                                  score_fxn = TrackSig.options()$pelt_score_fxn,
                                  bin_size = TrackSig.options()$bin_size)
{
  n_bins <- ncol(vcf)
  n_sigs <- ncol(alex.t)

  penalty <- eval(penalty)

  # aggregate bin summary stats
  phis <- aggregate(vcaf$phi, by = list(vcaf$binAssignment), FUN = sum)$x
  quadratic_phis <- aggregate(vcaf$phi, by = list(vcaf$binAssignment), FUN = function(x){return(sum(x^2))})$x

  # allow vaf permutation
  if(TrackSig.options()$permute_vafs){

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

      # allow vaf permutation
      if(TrackSig.options()$permute_vafs){
        r_seg_qis <- vcaf$phi2[vcaf$binAssignment %in% (sp_slice[1] : sp_slice[2])]
      }

      r_seg_qis <- r_seg_qis / (2 + vcaf$purity[1] * (r_seg_q_cni - 2))
      r_seg_qis <- unlist(lapply(r_seg_qis, 1, FUN = min))

      r_seg_vi <- vcaf$vi[vcaf$binAssignment %in% (sp_slice[1] : sp_slice[2])]
      r_seg_ri <- vcaf$ri[vcaf$binAssignment %in% (sp_slice[1] : sp_slice[2])]

      r_seg_counts <- rowSums(vcf[, sp_slice[1] : sp_slice[2], drop = FALSE])

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
