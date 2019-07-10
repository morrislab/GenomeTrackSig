# VCF_cost.R
# VCF cost functions for likelihood
# Author: Cait Harrigan

#Gaussian likelihood maximization
gaussian_ll <- function(phis, quad_phis, bin_size, ...){
  # phis read from counts file
  # quadratic_phis read from quadratic phis file
  # ... allows for various function signatures
  # Score a segment using likelihood under normal

  N <- length(phis) * bin_size
  sigmasq <- mean(quad_phis) - mean(phis)^2

  if (mean(quad_phis) < mean(phis)^2){
    warning(sprintf("mean quad_phis is %s, mean phis^2 is %s", mean(quad_phis), mean(phis)^2))
    sigmasq <- max((sigmasq), 10^-8, na.rm = T)
  }

  LL <- (-N / 2) * (log(2 * pi * sigmasq) + 1)

  return(LL)
}

# log_likelihood_mixture_multinomials defined in mixture_of_multinomials.R
sum_gaussian_mixture_ll <- function(phis, quad_phis, multinomial_vector,
                                                 composing_multinomials, mixtures, bin_size, ...){

  return( sum(log_likelihood_mixture_multinomials(multinomial_vector, composing_multinomials, mixtures),
              gaussian_ll(phis, quad_phis, bin_size)) )
}

#Poisson likelihood maximization
poisson_ll <- function(phis, quad_phis, bin_size, ...){
  # phis read from counts file
  # quadratic_phis read from quadratic phis file
  # ... allows for various function signatures
  # Score a segment using likelihood under normal

  N <- length(phis) * bin_size
  t <- N * mean(quad_phis)
  s <- N * mean(phis)
  y <- (1/2) * (sqrt((t/N) + (1/4)) - 1) #y = mean = variance under poisson conditions

  warning(sprintf("y value is %s", y))

  assertthat::assert_that(y > 0 , sprintf("y value is %s", y))

  LL <- (-N / 2) * ( log(2 * pi * y) + y) - (t/(2*y)) + s

  return(LL)
}

sum_poisson_mixture_ll <- function(phis, quad_phis, multinomial_vector,
                                    composing_multinomials, mixtures, bin_size, ...){

  return( sum(log_likelihood_mixture_multinomials(multinomial_vector, composing_multinomials, mixtures),
              poisson_ll(phis, quad_phis, bin_size)) )
}

# [END]
