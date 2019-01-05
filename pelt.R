# pelt.R
# detect (segment) change points

pelt <- function(phis, scoreFun, read_depth = 50, no_cp_pen = 0){
  
  n <- length(phis)
  penalty <- 2 * log(n)
  
  # Initialize score matrix
  scoreMat <- matrix(NaN, ncol = n, nrow = n)
  
  # For each subproblem size > 1
  for (size in 2:n){
    
    # Score subproblem with no cp
    scoreMat[size, 1] <- scoreFun(phis[1:size], beta = no_cp_pen)
    
    
    # If there are taus in the subproblem (minimum scorable tau interval size is 2)
    if (size >= 4){
      
      # Score subproblem for each cp at position 2 <= tau <= size-2 within subproblem
      for (tau in 2:(size-2)){
        scoreMat[size, tau] <- scoreFun(phis[tau:size], beta = penalty)
      }
    }
    
    # Prune
    
  }
  return(scoreMat)
}
  
score_mle <- function(phis, beta = 0){
  # Score a segment using likelihood under normal
  n <- length(phis)
  LL <- sum( (n/2)*log(2*pi), (n/2)*log(var(phis)), 
             1/(2*var(phis)) * sum( (phis - mean(phis))^2 ) )
              
  return ( - beta - LL)
}

score_mle2 <- function(phis, beta = 0){
  # Score a segment using likelihood under normal
  n <- length(phis)
  LL <- sum( (n/2)*log(2*pi), (n/2)*log(var(phis)), 
             1/(2*var(phis)) * sum( (phis - mean(phis))^2 ) )
  
  return ( - beta - 2*LL)
}


score_const <- function(phis, beta){
  return(60)
}

#from yulia's codebase
recover_changepoints <- function(sp_score_matrix){
  changepoints <- c()
  continue <- TRUE
  current <- dim(sp_score_matrix)[1]
  while (continue){
    #max entry over last row
    prev <- which.max(sp_score_matrix[current, ])
    
    #stop if reach col 1 (cp pos 1)
    if (prev - 1 <= 1)
    {
      continue <- FALSE
    } 
    
    #add the best cp pos to cps
    else
    {
      changepoints <- c(prev - 1, changepoints)
    }
    
    #move up 1 row
    current <- prev - 1
  }
  return(changepoints)
}
# ============================================================
setwd("~/Documents/BCB430/repo/")
source("simPhi.R")
library(changepoint)

mutsPerPop <- 100
readDepth <- 15
phiIs = c(0.64, 0.35, 0.15)

simDat <- generate_mutations_binom(mutsPerPop, phiIs, readDepth)
simPhis <- generate_phis(simDat$ref_counts, simDat$read_depths)

scores_0 <- pelt(simPhis, score_mle)
scores_pen <- pelt(simPhis, score_mle, no_cp_pen = 3 * log(length(simPhis)))

#sumstat method - mu set to var when passed on later?
mu <- mean(simPhis)
sumstat <- cbind(c(0,cumsum(simPhis)),c(0,cumsum(simPhis^2)),cumsum(c(0,(simPhis-mu)^2)))
scores_sumstat <- pelt(sumstat, score_mle)

#scores[300,]

cp1 <- recover_changepoints(scores_0)
cp2 <- recover_changepoints(scores_pen)
cp3 <- cpt.meanvar(simPhis, penalty = "BIC", method = "PELT", class = F)
cp4 <- recover_changepoints(scores_sumstat)

plot(simPhis)

abline(v = cp1, col = 2)
abline(v = cp2, col = 3)
abline(v = cp3, col = 4)
abline(v = cp4, col = 5)



# [END]