# VCF_cost.R
# VCF cost functions for likelihood
# Author: Cait Harrigan


beta_mom <- function(qis, vis, ris, ...){

  n <- length(qis)

  mu_hat <- mean(qis)
  v_hat <- (1/(n-1)) * sum((qis - mu_hat)^2)

  assertthat::assert_that(v_hat < mu_hat * (1-mu_hat), msg = "method of moments assumption violated")

  a_mom <- mu_hat * (((mu_hat * (1 - mu_hat)) / v_hat) -1)
  b_mom <- (1-mu_hat) * (((mu_hat * (1 - mu_hat)) / v_hat) -1)

  LL <- (a_mom - 1) * sum(log(qis)) + (b_mom - 1) * sum(log(1 - qis)) - n*lbeta(a_mom, b_mom)

  print(c(a_mom, b_mom, LL))


  return(LL)

}


# beta likelihood maximization
beta_ll <- function(qis, ...){

  #qis are the VAFs for the subproblem

  n <- length(qis)

  assertthat::assert_that(length(qis) == length(vis), length(qis) == length(ris), msg = "problem subsetting is not good!")

  alpha <- sum(qis) + 1
  beta <- sum(1-qis) + 1

  LL <- lbeta(alpha, beta) + log(pbeta(max(qis), alpha, beta) - pbeta(min(qis), alpha, beta))

  print(c(alpha, beta, LL))

  return(LL)

}

sum_beta_mixture_ll <- function(qis, multinomial_vector,
                                composing_multinomials, mixtures, ...){

  return( sum(log_likelihood_mixture_multinomials(multinomial_vector, composing_multinomials, mixtures),
              beta_ll(qis)) )

}


#Gaussian likelihood maximization
gaussian_ll_MBIC <- function(phis, quad_phis, bin_size, ...){
  # phis read from counts file
  # quadratic_phis read from quadratic phis file
  # ... allows for various function signatures
  # Score a segment using likelihood under normal

  n <- length(phis) * bin_size
  sigmasq <- (sum(quad_phis) / n) - (sum(phis)/n)^2

  assertthat::assert_that((sum(quad_phis) / n) > (sum(phis)/n)^2,
                          msg = sprintf("mean quad_phis is %s, mean phis^2 is %s",
                                        mean(quad_phis), mean(phis)^2))

  LL <- (-n / 2) * (log(2 * pi * sigmasq) + 1)  -log(n)

  return(LL)
}

#Gaussian likelihood maximization
gaussian_ll <- function(phis, quad_phis, bin_size, ...){
  # phis read from counts file
  # quadratic_phis read from quadratic phis file
  # ... allows for various function signatures
  # Score a segment using likelihood under normal

  n <- length(phis) * bin_size
  sigmasq <- (sum(quad_phis) / n) - (sum(phis)/n)^2

  assertthat::assert_that((sum(quad_phis) / n) > (sum(phis)/n)^2,
                          msg = sprintf("mean quad_phis is %s, mean phis^2 is %s",
                                        mean(quad_phis), mean(phis)^2))

  LL <- (-n / 2) * (log(2 * pi * sigmasq) + 1)

  return(LL)
}

# log_likelihood_mixture_multinomials defined in mixture_of_multinomials.R
sum_gaussian_mixture_ll <- function(phis, quad_phis, multinomial_vector,
                                                 composing_multinomials, mixtures, bin_size, ...){

  return( sum(log_likelihood_mixture_multinomials(multinomial_vector, composing_multinomials, mixtures),
              gaussian_ll(phis, quad_phis, bin_size)) )
}


# log_likelihood_mixture_multinomials defined in mixture_of_multinomials.R
sum_gaussian_MBIC_mixture_ll <- function(phis, quad_phis, multinomial_vector,
                                    composing_multinomials, mixtures, bin_size, ...){

  return( sum(log_likelihood_mixture_multinomials(multinomial_vector, composing_multinomials, mixtures),
              gaussian_ll_MBIC(phis, quad_phis, bin_size)) )
}

#Poisson likelihood maximization
poisson_ll <- function(phis, quad_phis, bin_size, ...){
  # phis read from counts file
  # quadratic_phis read from quadratic phis file
  # ... allows for various function signatures
  # Score a segment using likelihood under normal

  N <- length(phis) * bin_size
  t <- sum(quad_phis) / N
  s <- sum(phis) / N
  y <- (1/2) * (sqrt((4*t/N) + 1) - 1) #y = mean = variance under poisson conditions

  assertthat::assert_that(y > 0 , msg = sprintf("y value is %s", y))

  LL <- (-N / 2) * ( log(2 * pi * y) + y) - (t/(2*y)) + s

  return(LL)
}

sum_poisson_mixture_ll <- function(phis, quad_phis, multinomial_vector,
                                    composing_multinomials, mixtures, bin_size, ...){

  return( sum(log_likelihood_mixture_multinomials(multinomial_vector, composing_multinomials, mixtures),
              poisson_ll(phis, quad_phis, bin_size)) )
}


cpt.gamma <- function(phis, quad_phis, bin_size,...){


  n <- length(phis) * bin_size
  x <- sum(phis)

  return(2*n*(log(x)-log(n)) - log(n))
}

# [END]
