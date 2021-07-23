
# Removing SBS1 and SBS5
no_aging <- TrackSig::alex_merged %>%
  dplyr::select(-SBS1, -SBS5)

# Read in data
aml_master <- readr::read_csv('~/Desktop/Caitlin/data/Myeloid-AML_pooled.csv')
pilo_master <- readr::read_csv('~/Desktop/Caitlin/data/CNS-PiloAstro_pooled.csv')
oligo_master <- readr::read_csv('~/Desktop/Caitlin/data/CNS-Oligo_pooled.csv')

# Detect active signatures
aml_sigs <- detectActiveSignatures(aml_master, binSize = 200)
pilo_sigs <- detectActiveSignatures(pilo_master, binSize = 200)
oligo_sigs <- detectActiveSignatures(oligo_master, binSize = 200)

# Active + predefined sigs - aging sigs
aml_all <- c("SBS7", "SBS6", "SBS29", "SBS18")
pilo_all <- c("SBS7", "SBS6", "SBS29", "SBS19", "SBS23", "SBS40")
oligo_all <- c("SBS7", "SBS6", "SBS29", "SBS8", "SBS40")

# Run TrackSig
aml_active_traj <- SpaceTrack(path = aml_master, bootstrapSamples = 15, parallelize = TRUE, activeInSample = aml_all, binSize = 200)
base::saveRDS(aml_active_traj, file = "aml_traj_all.rds")
pilo_active_traj <- SpaceTrack(path = pilo_master, bootstrapSamples = 15, parallelize = TRUE, activeInSample = pilo_all, binSize = 200)
base::saveRDS(pilo_active_traj, file = "pilo_traj_all.rds")
oligo_active_traj <- SpaceTrack(path = oligo_master, bootstrapSamples = 15, parallelize = TRUE, activeInSample = oligo_all, binSize = 200)
base::saveRDS(oligo_active_traj, file = "oligo_traj_all.rds")

# Plot trajectories
plotSpaceTrajectory(aml_active_traj) + labs(title = "Myeloid-AML | n = 11 | Binsize = 200 mutations | Total mutations = 13,849 | 15 bootstrap samples")
plotSpaceTrajectory(pilo_active_traj) + labs(title = "CNS-PiloAstro | n = 89 | Binsize = 200 mutations | Total mutations = 22,137 | 15 bootstrap samples")
plotSpaceTrajectory(oligo_active_traj) + labs(title = "CNS-Oligo | n = 18 | Binsize = 200 mutations | Total mutations = 48,836 | 15 bootstrap samples")
