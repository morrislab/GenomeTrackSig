# VCF_cost.R
# VCF cost functions for likelihood
# Author: Cait Harrigan

#Gaussian likelihood maximization
gaussian_ll <- function(phis, quad_phis, ...){
  # phis read from counts file
  # quadratic_phis read from quadratic phis file
  # ... allows for various function signatures
  # Score a segment using likelihood under normal

  n <- length(phis)

  phi = sum(phis)
  quad_phi = sum(quad_phis)

  mu = phi / n
  sigmaSq = (quad_phi / n) - (mu^2)

  LL <- sum( (-1/2) * n * log(2*pi*sigmaSq), -sum((phis - mu)^2) / (2*sigmaSq) )
  return(LL)
}

# log_likelihood_mixture_multinomials defined in mixture_of_multinomials.R
sum_gaussian_mixture_multinomials_ll <- function(phis, quadratic_phis, multinomial_vector,
                                                 composing_multinomials, mixtures, ...){
  return( sum(log_likelihood_mixture_multinomials(multinomial_vector, composing_multinomials, mixtures)
              , gaussian_ll(phis, quadratic_phis)) )
}

# [END]
