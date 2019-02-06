# VCF_cost.R
# VCF cost functions for likelihood
# Author: Cait Harrigan

#Gaussian likelihood maximization
gaussian_ll <- function(phis, quadratic_phis, ...){
  # phis read from counts file
  # quadratic_phis read from quadratic phis file
  # ... allows for various function signatures
  # Score a segment using likelihood under normal
  n <- length(phis)
  LL <- sum( (n/2)*log(2*pi), (n/2)*log(quadratic_phis),
             1/(2*quadratic_phis) * sum( (phis - mean(phis))^2 ) )
  return(LL)
}

# log_likelihood_mixture_multinomials defined in mixture_of_multinomials.R
sum_gaussian_mixture_multinomials_ll <- function(phis, quadratic_phis, multinomial_vector,
                                                 composing_multinomials, mixtures, ...){
  return( sum(log_likelihood_mixture_multinomials(multinomial_vector, composing_multinomials, mixtures)
              , gaussian_ll(phis, quadratic_phis)) )
}

# [END]
