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
  sigmasq <- mean(quad_phis) - mean(phis)^2

  if (mean(quad_phis) >= mean(phis)^2){
    warning(sprintf("mean quad_phis is %s, mean phis^2 is %s", mean(quad_phis), mean(phis)^2))
  }

  scale <- max(sqrt(2 * pi * sigmasq), 10^-8, na.rm = T)

  LL <- - (S * B * log(scale)) - (S * B / 2)

  return(LL)
}

# log_likelihood_mixture_multinomials defined in mixture_of_multinomials.R
sum_gaussian_mixture_multinomials_ll <- function(phis, quad_phis, multinomial_vector,
                                                 composing_multinomials, mixtures, bin_size, ...){

  return( sum(log_likelihood_mixture_multinomials(multinomial_vector, composing_multinomials, mixtures),
              gaussian_ll(phis, quad_phis, bin_size)) )
}

# [END]
