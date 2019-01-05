#simBench.R
#benchmark tracksig edits

# ===============
# Set-up 
# ===============
setwd("~/Documents/BCB430/repo/")
source("simPhi.R")
source("pelt.R")
library(changepoint)

# ===============
# Parameters
# ===============
nSim <- 10

#number of subpopulations to simulate
poisRep <- rpois(nSim, 0.7)
poisRep[poisRep==0] <- rpois(length(poisRep[poisRep==0]),0.7)
#add one for single populations (no subclone boundary)
nSubClones <- poisRep + 1

#CCF's list with frequencies as entries
ccfReps <- list()
for (i in 1:nSim){
    ccfReps[[i]] <- runif(nSubClones[i])
}

mutsPerPop <- 100
readDepth <- 50

# ===============
# Run sims
# ===============

#for each simulation, generate phis, guess with various methods,
#store stats on best method

#package changepoint methods
cps_cpt.mean <- vector("list", nSim)
cps_cpt.meanvar <- vector("list", nSim)

#cait implementations
cps_cait1 <- vector("list", nSim)
cps_cait2 <- vector("list", nSim)
cps_cait3 <- vector("list", nSim)

for (i in 1:nSim){
  simDat <- generate_mutations_binom(mutsPerPop, ccfReps[[i]], readDepth)
  simPhis <- generate_phis(simDat$ref_counts, simDat$read_depths)
  
  cps_cpt.mean[[i]] <- cpt.mean(simPhis, penalty = "BIC",  method = "PELT", class=F)
  cps_cpt.meanvar[[i]] <- cpt.meanvar(simPhis, penalty = "BIC",  method = "PELT", class=F)
  cps_cait1[[i]] <- pelt(simPhis, score_mle)
  cps_cait1[[i]] <- pelt(simPhis, score_mle2)
  cps_cait1[[i]] <- pelt(simPhis, score_const)
}

# ===============
# Compare
# ===============






#[END]