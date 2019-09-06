# VCF_cost.R
# VCF cost functions for likelihood
# Author: Cait Harrigan



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


# [END]
