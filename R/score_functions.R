# VCF_cost.R
# VCF cost functions for likelihood
# Author: Cait Harrigan

#Gaussian likelihood maximization
gaussian_ll <- function(phis, quad_phis, bin_size, ...){
  # phis read from counts file
  # quadratic_phis read from quadratic phis file
  # ... allows for various function signatures
  # Score a segment using likelihood under normal

  B <- length(phis)
  S <- bin_size
  sigmasq <- ( mean(quad_phis) - mean(phis)^2 )

  LL <- (S*B/2) * (log(2 * pi * sigmasq) - 1)

  # multiply likelihood by bin_size for scaling
  return(LL)
}

# log_likelihood_mixture_multinomials defined in mixture_of_multinomials.R
sum_gaussian_mixture_multinomials_ll <- function(phis, quad_phis, multinomial_vector,
                                                 composing_multinomials, mixtures, bin_size, ...){

  return( sum(log_likelihood_mixture_multinomials(multinomial_vector, composing_multinomials, mixtures),
              gaussian_ll(phis, quad_phis, bin_size)) )
}

# [END]
