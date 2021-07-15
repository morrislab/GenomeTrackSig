library(ggplot2)
# Make masters

summary <- pcawg_cancer_types %>%
  dplyr::group_by(X1) %>%
  dplyr::summarize(n = dplyr::n())

# pilo_master <- poolSamples("~/Desktop/Caitlin/archive",
#                               "~/Desktop/Caitlin/data/pcawg_cancer_types.csv",
#                               "CNS-PiloAstro")
# kidney_master <- poolSamples("~/Desktop/Caitlin/archive",
#                               "~/Desktop/Caitlin/data/pcawg_cancer_types.csv",
#                               "Kidney-ChRCC")
# myeloid_master <- poolSamples("~/Desktop/Caitlin/archive",
#                               "~/Desktop/Caitlin/data/pcawg_cancer_types.csv",
#                               "Myeloid-MPN")
# thyroid_master <- poolSamples("~/Desktop/Caitlin/archive",
#                               "~/Desktop/Caitlin/data/pcawg_cancer_types.csv",
#                               "Thy-AdenoCA")
# medullo_master <- poolSamples("~/Desktop/Caitlin/archive",
#                               "~/Desktop/Caitlin/data/pcawg_cancer_types.csv",
#                               "CNS-Medullo")
# endocrine_master <- poolSamples("~/Desktop/Caitlin/archive",
#                                 "~/Desktop/Caitlin/data/pcawg_cancer_types.csv",
#                                 "Panc-Endocrine")
# osteo_master <- poolSamples("~/Desktop/Caitlin/archive",
#                             "~/Desktop/Caitlin/data/pcawg_cancer_types.csv",
#                             "Bone-Osteosarc")
# aml_master <- poolSamples("~/Desktop/Caitlin/archive",
#                               "~/Desktop/Caitlin/data/pcawg_cancer_types.csv",
#                               "Myeloid-AML")

no_aging <- TrackSig::alex_merged %>%
  dplyr::select(-SBS1, -SBS5)
# Myeloid-AML
aml_master <- readr::read_csv('~/Desktop/Caitlin/data/Myeloid-AML_pooled.csv')
aml_sigs <- c("SBS1", "SBS5", "SBS18")
aml_active <- detectActiveSignatures(aml_master, referenceSignatures = no_aging, threshold=0.05, binSize=200)
aml_traj <- SpaceTrack('~/Desktop/Caitlin/data/Myeloid-AML_pooled.csv',
                        bootstrapSamples = 5, activeInSample = aml_sigs, binSize=200, parallelize = TRUE)
base::saveRDS(aml_traj, file = "aml_5.rds")

# CNS-PiloAstro
pilo_sigs <- c("SBS1", "SBS5", "SBS19", "SBS23", "SBS40")
pilo_master <- readr::read_csv("~/Desktop/Caitlin/data/CNS-PiloAstro_pooled.csv")
pilo_active <- detectActiveSignatures(pilo_master, referenceSignatures = no_aging, threshold=0.05, binSize=200)
pilo_traj <- SpaceTrack('~/Desktop/Caitlin/data/CNS-PiloAstro_pooled.csv',
                         bootstrapSamples = 30, activeInSample = pilo_active, binSize=200, parallelize = TRUE)
base::saveRDS(pilo_traj, file = "pilo_30.rds")

# Myeloid-MPN
myeloid_sigs <- c("SBS1", "SBS5", "SBS19", "SBS23", "SBS40")
myeloid_traj <- SpaceTrack('~/Desktop/Caitlin/data/Myeloid-MPN_pooled.csv',
                        bootstrapSamples = 0, activeInSample = myeloid_sigs, binSize=200, parallelize = TRUE, chr_level = TRUE)
base::saveRDS(myeloid_traj, file = "myeloid_30.rds")

# Thy-AdenoCA
thyroid_sigs <- c("SBS1", "SBS2.13", "SBS5", "SBS40")
thyroid_traj <- SpaceTrack('~/Desktop/Caitlin/data/Thy-AdenoCA_pooled.csv',
                           bootstrapSamples = 30, activeInSample = thyroid_sigs, binSize=200, parallelize = TRUE)
base::saveRDS(thyroid_traj, file = "thyroid_30.rds")

# CNS-Medullo
medullo_sigs <- c("SBS1", "SBS5", "SBS8", "SBS18", "SBS39", "SBS40")
medullo_traj <- SpaceTrack('~/Desktop/Caitlin/data/CNS-Medullo_pooled.csv',
                           bootstrapSamples = 30, activeInSample = medullo_sigs, binSize=200, parallelize = TRUE)
base::saveRDS(medullo_traj, file = "medullo_30.rds")

# Kidney-ChRCC
kidney_sigs <- c("SBS1", "SBS2.13", "SBS5", "SBS17", "SBS29", "SBS40")
kidney_traj <- SpaceTrack('~/Desktop/Caitlin/data/Kidney-ChRCC_pooled.csv',
                           bootstrapSamples = 30, activeInSample = kidney_sigs, binSize=200, parallelize = TRUE)
base::saveRDS(kidney_traj, file = "kidney_30.rds")

# CNS-Oligo
oligo_sigs <- c("SBS1", "SBS5", "SBS8", 'SBS40')
oligo_traj <- SpaceTrack('~/Desktop/Caitlin/data/CNS-Oligo_pooled.csv',
                         bootstrapSamples = 30, activeInSample = oligo_sigs, binSize=200, parallelize = TRUE)
base::saveRDS(oligo_traj, file = "oligo_30.rds")
# plotSpaceTrajectory(oligo_traj) +
#   labs(title = "CNS-Oligo (n = 18) | Total mutations = 48,836 | Binsize = 200 mutations | 30 bootstrap samples")

# Cervix-SCC
cervix_sigs <- c("SBS1", "SBS2.13", "SBS5", "SBS18", "SBS40")
cervix_traj <- SpaceTrack('~/Desktop/Caitlin/data/Cervix-SCC_pooled.csv',
                         bootstrapSamples = 30, activeInSample = cervix_sigs, binSize=200, parallelize = TRUE)
base::saveRDS(cervix_traj, file = "cervix_30.rds")
# plotSpaceTrajectory(cervix_traj, chr_level=F, cutoff=.5) +
#   labs(title = "Cervix-SCC (n = 18) | Total mutations = 114,879 | Binsize = 200 mutations | 20 bootstrap samples")

# Lymph-CLL
lymph_sigs <- c('SBS1', 'SBS5', 'SBS9', 'SBS40')
lymph_traj <- SpaceTrack('~/Desktop/Caitlin/data/Lymph-CLL_pooled.csv',
                         bootstrapSamples = 30, activeInSample = lymph_sigs, binSize=200, parallelize = TRUE)
base::saveRDS(lymph_traj, file = "lymph_30.rds")
# plotSpaceTrajectory(lymph_traj, cutoff=0.5) +
#   labs(title = "Lymph-CLL (n = 94) | Total mutations = 228,153 | Binsize = 300 mutations | 20 bootstrap samples")

# Panc-Endocrine
endocrine_sigs <- c("SBS1", "SBS2.13", "SBS3", "SBS5", "SBS6", "SBS8", "SBS9", "SBS11", "SBS26", "SBS30", "SBS36", "SBS39")
endocrine_traj <- SpaceTrack("~/Desktop/Caitlin/data/Panc-Endocrine_pooled.csv",
                             bootstrapSamples = 30, activeInSample = endocrine_sigs, binSize = 200, parallelize = TRUE)
base::saveRDS(endocrine_traj, file = "endocrine_30.rds")

# Bone-Osteosarc
osteo_sigs <- c("SBS1", "SBS2.13", "SBS3", "SBS5", "SBS8", "SBS17", "SBS30", "SBS40")
osteo_traj <- SpaceTrack("~/Desktop/Caitlin/data/Bone-Osteosarc_pooled.csv",
                         bootstrapSamples = 30, activeInSample = osteo_sigs, binSize = 200, parallelize = TRUE)
base::saveRDS(osteo_traj, file = "osteo30.rds")


# Plotting
trajPlot <- plotSpaceTrajectory(aml_traj)
plotSpaceTrajectory(pilo_30) +
  labs(title = "CNS-PiloAstro (n = 89) | Total mutations = 22,137 | Binsize = 200 mutations | 30 bootstrap samples")
plotSpaceTrajectory(thyroid_30) +
  labs(title = "Thyroid-AdenoCA (n = 48) | Total mutations = 65,533 | Binsize = 200 mutations | 30 bootstrap samples")
plotSpaceTrajectory(aml_traj) + labs(title = "Myeloid-AML (n = 11) | Total mutations = 13,849 | Binsize = 200 mutations | 30 bootstrap samples")



myeloid_master <- readr::read_csv('~/Desktop/Caitlin/data/Myeloid-MPN_pooled.csv')
gcPlot <- makeGCPlot(master = aml_master, binSize = 200, trajectory = aml_traj)

densityPlot <- makeDensityPlot(aml_master, 200, aml_traj)
addGenomicFeatures(trajPlot, aml_master, 200, aml_traj)

