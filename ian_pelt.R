# pelt.R
# detect (segment) change points

#from Ian's python

pelt <- function(phis, read_depth, score_function, n_sigs = 6){
  
  n_bins <- length(phis)
  
  # TODO: Figure out how to properly penalize segment.
  penalty <- n_sigs * log(n_bins * read_depth)
  
  # Matrix which contains scores for all subpartitions
  sp_scores <- matrix(NaN, ncol = n_bins, nrow = n_bins)
  
  # Array containing maximum score for subpartition of length index
  max_sp_scores <- array(0, dim = n_bins)
  
  # Set of potential changepoints to disclude
  prune_set <- c()
  
  for (sp_len in 1:(n_bins)){
    for (last_cp in 1:(sp_len)){
      
      r_phis <- phis[last_cp:sp_len] 

      if (last_cp %in% prune_set){
        sp_scores[sp_len, last_cp] <- -Inf
        next
      }
      
      # Segments of length 1 are not considered
      if (length(r_phis) <= 1){
        sp_scores[sp_len, last_cp] <- -Inf
        next
      }
      
      # Currently using left and right most phi in segment as integration bounds.
      # Changing to mean between bounds does not greatly affect result.
      r_seg_score = 2 * score_function(r_phis)
      
      # Do not penalize if not introducing changepoint
      if (last_cp == 1){
        sp_scores[sp_len, last_cp] <- r_seg_score
      }
      
      else{
        l_score = max_sp_scores[last_cp-1]
        sp_scores[sp_len, last_cp] <- l_score + r_seg_score - penalty
      }
    }
  
    # Check all changepoints for prune condition
    max_sp_scores[sp_len] <- max(sp_scores[sp_len,], na.rm = T)
    
    for (cp in 1:sp_len){
      if (is.infinite(sp_scores[sp_len, cp])){
        next
      }
      if (sp_scores[sp_len, cp] + penalty < max_sp_scores[sp_len]){
        prune_set <- c(prune_set, cp)
      }
    }
  }
  
  
  return (sp_scores)
}

score_mle <- function(phis){
 # Score a segment using maximum likelihood

  return (prod(dnorm(phis), mean(phis), sqrt(var(phis))))
}

score_const <- function(phis){
  return(60)
}

#from yulia's codebase
recover_changepoints <- function(sp_score_matrix){
  changepoints <- c()
  continue <- TRUE
  current <- dim(sp_score_matrix)[1]
  while (continue){
    #max over row
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
# ============================================================
setwd("~/Documents/BCB430/repo/")
source("simPhi.R")

sp_scores <- pelt(simPhis, read_depth, score_const)
sp_scores[1:8,1:8]

#lots of NaN, but some not
length(which(is.nan(sp_scores) != 0))
dim(sp_scores)

sp_scores[300,]

recover_changepoints(sp_scores)

# [END]