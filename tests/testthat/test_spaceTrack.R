
myCluster <- parallel::makeCluster(2, type = "FORK")
doParallel::registerDoParallel(myCluster)
############

library(ggplot2)

# CNS-Oligo
oligo_sigs <- c("SBS1", "SBS5", "SBS8", 'SBS40')
oligo_master <- readr::read_csv("~/Desktop/CBSP2021/data/CNS-Oligo_pooled.csv")

# Bootstrapping

cervix_sigs <- c("SBS1", "SBS2.13", "SBS5", "SBS18", "SBS40")
cervix_master <- readr::read_csv('~/Desktop/CBSP2021/data/Cervix-SCC_pooled.csv')
cervix_bootstraps <- bootstrapActivities(master = cervix_master, activeInSample = cervix_sigs,
                                         binSize = 200, nSamples = 20)

plotBootstrap(cervix_bootstraps) +
  labs(title = "Cervix-SCC (n = 18) | Total mutations = 114,879 | Binsize = 200 mutations | 20 bootstrap samples")


oligo_bootstraps <- bootstrapActivities(master = oligo_master, activeInSample = oligo_sigs,
                                        binSize = 200, nSamples = 20)

plotBootstrap(oligo_bootstraps) +
  labs(title = "CNS-Oligo (n = 18) | Total mutations = 48,836 | Binsize = 200 mutations | 20 bootstrap samples")






