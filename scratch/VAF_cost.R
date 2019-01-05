# VAF_cost.R
# VAF cost functions for likelihood
# Author: Cait Harrigan

#from package changepoint mll_meanvar()
#raw phis must be extracted in make_hundreds.py
cpt_mll_meanvar <- function(mean_phis, mean_square_phis){ 
  # x is sumstat[j]
  # x2 is sumstat[total size n + j + 1]
  # j is length of current segment
  # uses cumsum
  sigsq <- (mean_square_phis - ((phis^2)/n))/n 

  #ensure non-negative
  if(sigsq < 0 ){
  	sigsq <- 1+e-10
  }
  
  return(n * ( log(2*pi) + log(sigsq) + 1 ))
}



#Gaussian likelihood maximization 
gaussian_mll <- function(phis, mean_square_phis){
	# mean_phis read from counts file
	# mean_square_phis read from quadratic phis file
  	# Score a segment using likelihood under normal
 	n <- length(phis)
  LL <- sum( (n/2)*log(2*pi), (n/2)*log(mean_square_phis), 
             1/(2*mean_square_phis) * sum( (phis - mean(phis))^2 ) )
  return(LL)
}

# Raw VCF PELT
raw_mll <- function(phis){
  # Score a segment using likelihood under normal
  n <- length(phis)
  LL <- sum( (n/2)*log(2*pi), (n/2)*log(var(phis)), 
             1/(2*var(phis)) * sum( (phis - mean(phis))^2 ) )
}


# [END]