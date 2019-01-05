# VAFpelt.R

#written with reference to Ian's segement.py and
#yulia's change_points.R, and the consensus supplement

library(changepoint)

setwd("~/Documents/BCB430/repo/")
source("simPhi.R")

cps <- cpt.meanvar(simPhis, penalty = "BIC", method = "PELT")
cps@param.est
plot(cps)

#==========================================================

source("simPhi.R")
data <- simPhis

costfunc = "meanvar.norm"
diffparam <- 2
n=length(data)
minseglen <- 2
penalty <- "BIC"
mul.method <- "PELT"

if(n<(2*minseglen)){stop('Minimum segment legnth is too large to include a change in this data')}
pen.value = penalty_decision(penalty, pen.value, n, diffparam, asymcheck = costfunc, method=mul.method)
out=list()
out = data_input(data=data,method=mul.method,pen.value=pen.value,costfunc=costfunc,minseglen=minseglen,Q=Q)
out[[2]]

  
a <- cpt.meanvar(data, penalty = "BIC",  method = "PELT", class=F)
a


VAFsole <- function(phis){
  
  #penalty from changepoint 
  #meanvar.norm(), penalty_decision(), data_input()
  
  costfunc = "meanvar.norm"
  n=length(data)
  minseglen <- 2
  penalty <- "BIC"
  mul.method <- "PELT"
  
  if(n<(2*minseglen)){stop('Minimum segment legnth is too large to include a change in this data')}
  pen.value =  3 * log(n)
  out=list()
  out = data_input(data=data,method=mul.method,pen.value=pen.value,costfunc=costfunc,minseglen=minseglen,Q=Q)
  out[[2]]
  
}
  
# [END]