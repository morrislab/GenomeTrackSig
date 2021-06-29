library(foreach)
library(doParallel)
myCluster <- makeCluster(2, type = "FORK")
registerDoParallel(myCluster)
############

library(ggplot2)

# Cervix-SCC
path <- "~/Desktop/CBSP2021/Cervix-SCC_pooled.csv"

cervical_sigs <- c("SBS1", "SBS2.13", "SBS5", "SBS18", "SBS40")
counts <- binningNmut(path, binSize = 200)

cervix_results_bin <- foreach(i = c(1:23), .combine = 'c') %dopar% trackParallel(counts, i, cervical_sigs)

cervix_traj <- combineTraj(cervix_results_bin)
cervix_traj[['changepoints']] <- sort(cervix_traj[['changepoints']])
cervix_traj[['binData']] <- cervix_traj[['binData']] %>%
  dplyr::arrange(dplyr::desc(bin))

plotTrajectory(cervix_traj, linearX = T) + labs(title = "Cervix-SCC (n=18), Binsize = 200 mutations")

# CNS-Oligo

oligo_sigs <- c("SBS1", "SBS5", "SBS8", "SBS40")

oligo_master <- poolSamples(archivePath = "~/Desktop/CBSP2021/archive", typesPath = "~/Desktop/CBSP2021/pcawg_cancer_types.csv",
                            cancerType = "CNS-Oligo")
oligo_counts <- binningNmut(path = "~/Desktop/CBSP2021/CNS-Oligo_pooled.csv", binSize = 200)

oligo_results <- foreach(i = c(1:23), .combine = 'c') %do% trackParallel(oligo_counts, i, oligo_sigs)
oligo_traj <- combineTraj(oligo_results)
oligo_traj[['changepoints']] <- sort(oligo_traj[['changepoints']])
oligo_traj[['binData']] <- oligo_traj[['binData']] %>%
  dplyr::arrange(dplyr::desc(bin))

plotTrajectory(oligo_traj, linearX = T) + labs(title = "CNS-Oligo (n=18), Binsize = 200 mutations")

# Prostate Adeno-CA
prostate_sigs <- c("SBS1", "SBS2.13", "SBS5", "SBS8", "SBS18", "SBS33", "SBS37", "SBS40", 'SBS41')

prostate_master <- poolSamples(archivePath = "~/Desktop/CBSP2021/archive", typesPath = "~/Desktop/CBSP2021/pcawg_cancer_types.csv",
                             cancerType = "Prost-AdenoCA")

prostate_counts <- binningNmut(path = "Prost-AdenoCA_pooled.csv", binSize = 200)

prostate_results <- foreach(i = c(1:23), .combine = 'c') %dopar% trackParallel(prostate_counts, i, prostate_sigs)

prostate_traj <- combineTraj(prostate_results)
prostate_traj[['changepoints']] <- sort(prostate_traj[['changepoints']])



pcawg_cancer_types <- read_csv("~/Desktop/CBSP2021/pcawg_cancer_types.csv") %>%
  dplyr::group_by(X1) %>%
  dplyr::summarise(n = dplyr::n()) %>%
  dplyr::arrange(n)
