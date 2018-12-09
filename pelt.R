# pelt.R
# detect (segment) change points

pelt <- function(phis, read_depth, scoreFun){
  
  n <- length(phis)
  penalty <- 3 * log(n)
  
  # Initialize score matrix
  scoreMat <- matrix(NaN, ncol = n, nrow = n)
  
  # For each subproblem size > 1
  for (size in 2:n){
    
    # Score subproblem with no cp
    scoreMat[size, 1] <- scoreFun(phis[1:size], beta = 0)
    
    
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


score_const <- function(phis){
  return(60)
}

#from yulia's codebase
recover_changepoints <- function(sp_score_matrix){
  changepoints <- c()
  continue <- TRUE
  current <- dim(sp_score_matrix)[1]
  while (continue){
    #max entry over row
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
    
    #move left 1 col
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