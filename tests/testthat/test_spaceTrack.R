
myCluster <- parallel::makeCluster(15, type = "FORK")
doParallel::registerDoParallel(myCluster)
############

library(ggplot2)

# CNS-Oligo
oligo_sigs <- c("SBS1", "SBS5", "SBS8", 'SBS40')
oligo_master <- poolSamples(archivePath = '~/Desktop/Caitlin/archive',
                            typesPath = '~/Desktop/Caitlin/pcawg_cancer_types.csv',
                            cancerType = 'CNS-Oligo')

# Bootstrapping whole genome
oligo_test <- SpaceTrack('~/Desktop/Caitlin/data/CNS-Oligo_pooled.csv', chr_level=TRUE,
                         bootstrapSamples = 10, activeInSample = oligo_sigs, binSize=200)
oligo_traj <- trackParallelGen(oligo_master, oligo_sigs, 200, 20)
cervix_traj <- trackParallelGen(cervix_master, cervix_sigs, 200, 20)
plotSpaceTrajectory(oligo_test, chr_level=TRUE) +
  labs(title = "CNS-Oligo (n = 18) | Total mutations = 48,836 | Binsize = 200 mutations | 20 bootstrap samples")

# Bootstrapping

cervix_sigs <- c("SBS1", "SBS2.13", "SBS5", "SBS18", "SBS40")
cervix_master <- poolSamples(archivePath = '~/Desktop/Caitlin/archive',
                             typesPath = '~/Desktop/Caitlin/pcawg_cancer_types.csv',
                             cancerType = 'Cervix-SCC')


plotSpaceTrajectory(cervix_traj, chr_level=F, cutoff=.5) +
  labs(title = "Cervix-SCC (n = 18) | Total mutations = 114,879 | Binsize = 200 mutations | 20 bootstrap samples")


# Lymph

lymph_sigs <- c('SBS1', 'SBS5', 'SBS9', 'SBS40')
lymph_master <- poolSamples(archivePath = '~/Desktop/Caitlin/archive',
                            typesPath = '~/Desktop/Caitlin/pcawg_cancer_types.csv',
                            cancerType = 'Lymph-CLL')
lymph_traj <- trackParallelGen(lymph_master, lymph_sigs, 300, 20)
plotSpaceTrajectory(lymph_traj, cutoff=0.5) +
  labs(title = "Lymph-CLL (n = 94) | Total mutations = 228,153 | Binsize = 300 mutations | 20 bootstrap samples")

# Skin
skin_sigs <- c('SBS1', 'SBS2.13', 'SBS5', 'SBS7', 'SBS17', 'SBS38', 'SBS40')
skin_master <- poolSamples(archivePath = '~/Desktop/Caitlin/archive',
                           typesPath = '~/Desktop/Caitlin/pcawg_cancer_types.csv',
                           cancerType = 'Skin-Melanoma')
skin_traj <- trackParallelGen(skin_master, skin_sigs, 300, 20)
