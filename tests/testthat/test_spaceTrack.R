
myCluster <- parallel::makeCluster(2, type = "FORK")
doParallel::registerDoParallel(myCluster)
############

library(ggplot2)

plotTrajectory(cervix_traj, linearX = T) + labs(title = "Cervix-SCC (n=18), Binsize = 200 mutations")
plotTrajectory(oligo_traj, linearX = T) + labs(title = "CNS-Oligo (n=18), Binsize = 200 mutations")
plotTrajectory(thyroid_traj, linearX = T) + labs(title = "Thyroid-AdenoCA (n = 48), Binsize = 200 mutations")

# Lymph-CLL

lymph_sigs <- c("SBS1", "SBS5", "SBS9", "SBS40")
lymph_master <- poolSamples(archivePath = "~/Desktop/CBSP2021/archive", typesPath = "~/Desktop/CBSP2021/pcawg_cancer_types.csv",
                              cancerType = "Lymph-CLL")

lymph_traj <- trackParallel(lymph_master, 1, lymph_sigs, 200)
lymph_traj <- foreach(i = c(1:23), .combine = 'c') %dopar% trackParallel(lymph_master, i, lymph_sigs, 200)
lymph_traj <- combineTraj(lymph_traj)


plotTrajectory(lymph_final, linearX = T) + labs(title = "Lymph-CLL (n = 94), Total mutations = 228,153, Binsize = 200 mutations")


traj <- cleanTraj(lymph_copy)
lymph_final <- combineTraj(traj)

# CNS-Oligo
oligo_sigs <- c("SBS1", "SBS5", "SBS8", 'SBS40')
lymph_bootstraps <- bootstrapActivities(lymph_master, lymph_sigs, 200, nSamples = 5)

# Bootstrapping

cervix_sigs <- c("SBS1", "SBS2.13", "SBS5", "SBS18", "SBS40")
cervix_master <- poolSamples("~/Desktop/CBSP2021/archive", "~/Desktop/CBSP2021/pcawg_cancer_types.csv",
                             "Cervix-SCC")
cervix_bootstraps <- bootstrapActivities(master = cervix_master, activeInSample = cervix_sigs,
                                         binSize = 200, nSamples = 5)

plotBootstrap(lymph_bootstraps) +
  labs(title = "Lymph-CLL (n = 94) | Total mutations = 228,153 | Binsize = 200 mutations | 5 bootstrap samples")


