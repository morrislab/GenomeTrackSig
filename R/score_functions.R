# VCF_cost.R
# VCF cost functions for likelihood
# Author: Cait Harrigan

#Gaussian likelihood maximization
gaussian_ll <- function(phis, quad_phis, bin_size, ...){
  # phis read from counts file
  # quadratic_phis read from quadratic phis file
  # ... allows for various function signatures
  # Score a segment using likelihood under normal

  phi = mean(phis)
  theta = mean(quad_phis)
  sigma = sqrt(max((theta - phi^2), 0.0001, na.rm = T))

  LL <- sum(log(dnorm(phis, mean = phi, sd = sigma)))

  # multiply likelihood by bin_size for scaling
  return(LL * bin_size)
}

# log_likelihood_mixture_multinomials defined in mixture_of_multinomials.R
sum_gaussian_mixture_multinomials_ll <- function(phis, quad_phis, multinomial_vector,
                                                 composing_multinomials, mixtures, bin_size, ...){

  return( sum(log_likelihood_mixture_multinomials(multinomial_vector, composing_multinomials, mixtures),
              gaussian_ll(phis, quad_phis, bin_size)) )
}

# [END]
