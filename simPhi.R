#simPhi.R

#from Ian's python

generate_mutations_binom <- function(num_ssm, clonal_freqs, read_depth){

#Generate num_ssm simulated reference read counts for each clonal population.
#
#Args:
#num_ssms (int): number of ssms per population
#clonal_freqs (array of float): the population frequency of each subpopulation
#read_depth (int): simulation read depth
#Returns: 
#Array of sampled reference reads, and sampled 

  sim_ref_counts <- vector()
  sim_read_depth <- vector()
  for (clonal_freq in clonal_freqs){
    
    if (clonal_freq < 0 | clonal_freq > 1){
    stop("clonal frequency must be between 0 and 1")
    }
  
    for (i in 1:num_ssm){
      num_reads <- rpois(1, lambda = read_depth)
      ref_reads <- rbinom(1, num_reads, prob = (1 - clonal_freq + 0.5 * clonal_freq)) 
      
      if (num_reads == ref_reads){
        ref_reads = num_reads * (1 - clonal_freq + 0.5 * clonal_freq)}
      
      sim_ref_counts <- c(sim_ref_counts, ref_reads)
      sim_read_depth <- c(sim_read_depth, num_reads)

    }
  }

simBinom <- list()
simBinom$ref_counts <- sim_ref_counts
simBinom$read_depths <- sim_read_depth
  
return (simBinom)
}


generate_phis <- function(ref_counts, read_depths){

#Generate phi values from reference read counts sampled from binomial distribution
#
#Args:
#ref_counts (int): number of sampled reference reads
#read_depths (int): total number of sampled reads
#Returns:
#Sorted array of phi values

phis <- vector()

  for (i in 1:(length(ref_counts))){
    r <- ref_counts[i]
    v <- read_depths[i] - r
    phi <-  2 * rbeta(1, v + 1, r + 1)
    phis <- c(phis, phi)
  }


return (phis)
}

plotPhis <- function(phis, clonal_freqs, num_ssms, read_depth){
  
  estVars <- vector()
  trueVars <- vector()
  plot(0:1, 0:1, type = "n")
  
  for (i in 1:length(clonal_freqs)){
    print(sprintf("Cluster at mean: %f", clonal_freqs[i]))
    cluster <- phis[(1:num_ssms) + (num_ssms * (i-1))]
    freq_i <- clonal_freqs[i]
    #iEstVar <- ((2*(2 - freq_i)* freq_i) / read_depth)
    iEstVar <- ((2 - freq_i)* freq_i / read_depth) + ((2 - freq_i) * freq_i / (read_depth + 1))
    print(sprintf("Estimated variance: %f", iEstVar))
    print(sprintf("Actual variance: %f", var(cluster)))
    
    points(cluster, cluster, col = i)
    
    estVars <- c(estVars, iEstVar)
    trueVars <- c(trueVars, var(cluster))
  }
  
  return(list(estVars, trueVars))
}
    
num_ssms <- 100
clonal_freqs <- c(0.64, 0.36, 0.16)
read_depth <- 50

simDat <- generate_mutations_binom(num_ssms, clonal_freqs, read_depth)
#hist(simDat$ref_counts)

simPhis <- generate_phis(simDat$ref_counts, simDat$read_depths)
hist(simPhis, col=4,  breaks = 50, xlim = c(0, 1),
     main = sprintf("phis simulated with clonal frequencies [%s], read depth = %i", 
                    paste(clonal_freqs, collapse = ", "), read_depth),
     ylab = "number of mutations",
     xlab = expression(phi)
)


plotPhis(simPhis, clonal_freqs, num_ssms, read_depth)

#[END]