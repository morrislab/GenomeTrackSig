# CNS-cervix

archivePath <- "~/Desktop/Caitlin/archive"
typesPath <- "~/Desktop/Caitlin/data/pcawg_cancer_types.csv"
cancerType <- "Cervix-SCC"
cervix_sigs <- c("SBS1", "SBS2.13", "SBS5", "SBS18", "SBS40")

types <- readr::read_csv(typesPath,
                         col_names = c('type', 'guid'))
types <- types[2:nrow(types), ]

# list of sample filenames for desired cancer type
cervix_files <- c(types$guid[types$type == cancerType])
set.seed(23)
cervix_sample <- base::sample(cervix_files, size = length(cervix_files), replace = FALSE)

# Find mutation counts in each sample
for (i in 1:length(cervix_sample)) {
  temp <- readFormat(path = paste("~/Desktop/Caitlin/archive/", cervix_sample[i], ".MBcounts.csv", sep=''))
  print(paste(i, sum(colSums(temp[,6:101])), sep = ": "))
}

small_samples <- cervix_sample[c(11,4,9,8,7,16)]
medium_samples <- cervix_sample[c(15,13,12,10,17,14)]
large_samples <- cervix_sample[c(5,3,18,6,1,2)]

# Pool samples
cervix_medium <- readr::read_csv(as.character(paste(archivePath, "/", medium_samples[1], ".MBcounts.csv", sep = "")),
                                 col_types = readr::cols("seqnames" = "c",
                                                         "strand" = "c",
                                                         .default = "d"))

# add counts from remaining samples into master file and delete remaining files from memory
for (i in 2:length(medium_samples)) {
  temp <- readr::read_csv(as.character(paste(archivePath, "/", medium_samples[i], ".MBcounts.csv", sep = "")),
                          col_types = readr::cols(seqnames = "c",
                                                  strand = "c",
                                                  .default = "d"))
  cervix_medium[, 6:101] <- cervix_medium[, 6:101] + temp[, 6:101]
  base::remove(temp)
}

# make data type consistent in seqnames column
# split seqnames column into start chrom and end chrom
# remove unnecessary variables
cervix_medium <- cervix_medium %>%
  dplyr::mutate(seqnames = dplyr::case_when(seqnames == "X" ~ 23,
                                            seqnames == "Y" ~ 24,
                                            TRUE ~ as.numeric(seqnames)),
                start_chrom = seqnames,
                end_chrom = seqnames) %>%
  dplyr::select(-seqnames, -strand) %>%
  base::subset(select = c(100,1,101,2:99))

# Run SpaceTrack on medium samples-- chromosomal level
cervix_medium_traj <- SpaceTrack(path = cervix_medium, activeInSample = cervix_sigs, binSize = 200, chr_level = FALSE, parallelize = TRUE, bootstrapSamples = 5)
# Plot medium samples
plotSpaceTrajectory(cervix_medium_traj) + labs(title = "Cervix-SCC | n = 6 | Binsize = 200 mutations | Total mutations = 36,765 | 5 bootstrap samples")




# Bladder-TCC

archivePath <- "~/Desktop/Caitlin/archive"
typesPath <- "~/Desktop/Caitlin/data/pcawg_cancer_types.csv"
cancerType <- "Bladder-TCC"
bladder_sigs <- c("SBS1", "SBS2.13", "SBS5", "SBS8", "SBS29", "SBS40")

types <- readr::read_csv(typesPath,
                         col_names = c('type', 'guid'))
types <- types[2:nrow(types), ]

# list of sample filenames for desired cancer type
bladder_files <- c(types$guid[types$type == cancerType])
set.seed(23)
bladder_sample <- base::sample(bladder_files, size = 5, replace = FALSE)

# Find mutation counts in each sample
for (i in 1:length(bladder_sample)) {
  temp <- readFormat(path = paste("~/Desktop/Caitlin/archive/", bladder_sample[i], ".MBcounts.csv", sep=''))
  print(paste(i, sum(colSums(temp[,6:101])), sep = ": "))
}

small_samples <- bladder_sample[c(22,4,17,2,19)]
small_samples <- bladder_sample[c(5,20,18,15,6)]

# Pool samples
bladder_small <- readr::read_csv(as.character(paste(archivePath, "/", large_samples[1], ".MBcounts.csv", sep = "")),
                          col_types = readr::cols("seqnames" = "c",
                                                  "strand" = "c",
                                                  .default = "d"))

# add counts from remaining samples into master file and delete remaining files from memory
for (i in 2:length(large_samples)) {
  temp <- readr::read_csv(as.character(paste(archivePath, "/", large_samples[i], ".MBcounts.csv", sep = "")),
                          col_types = readr::cols(seqnames = "c",
                                                  strand = "c",
                                                  .default = "d"))
  bladder_large[, 6:101] <- bladder_large[, 6:101] + temp[, 6:101]
  base::remove(temp)
}

# make data type consistent in seqnames column
# split seqnames column into start chrom and end chrom
# remove unnecessary variables
bladder_large <- bladder_large %>%
  dplyr::mutate(seqnames = dplyr::case_when(seqnames == "X" ~ 23,
                                            seqnames == "Y" ~ 24,
                                            TRUE ~ as.numeric(seqnames)),
                start_chrom = seqnames,
                end_chrom = seqnames) %>%
  dplyr::select(-seqnames, -strand) %>%
  base::subset(select = c(100,1,101,2:99))

# Run SpaceTrack on small samples-- chromosomal level
bladder_large_traj <- SpaceTrack(path = bladder_large, activeInSample = bladder_sigs, binSize = 200, chr_level = TRUE, parallelize = TRUE, bootstrapSamples = 10)
# Plot small samples
plotSpaceTrajectory(bladder_large_traj, chr_level = TRUE, cutoff=0.5) + labs(title = "Bladder-TCC | n = 5 | Binsize = 200 mutations | Total mutations = 222,761 | 10 bootstrap samples")


# Read in individual dataframes
bladder1 <- readFormat(path = paste("~/Desktop/Caitlin/archive/", bladder_sample[1], ".MBcounts.csv", sep=''))
sum(colSums(bladder1[,6:101]))

bladder2 <- readFormat(path = paste("~/Desktop/Caitlin/archive/", bladder_sample[2], ".MBcounts.csv", sep=''))
sum(colSums(bladder2[,6:101]))

bladder3 <- readFormat(path = paste("~/Desktop/Caitlin/archive/", bladder_sample[3], ".MBcounts.csv", sep=''))
sum(colSums(bladder3[,6:101]))

bladder4 <- readFormat(path = paste("~/Desktop/Caitlin/archive/", bladder_sample[4], ".MBcounts.csv", sep=''))
sum(colSums(bladder4[,6:101]))

bladder5 <- readFormat(path = paste("~/Desktop/Caitlin/archive/", bladder_sample[5], ".MBcounts.csv", sep=''))
sum(colSums(bladder5[,6:101]))

# Calculate trajectories for each sample
bladder1_traj <- SpaceTrack(bladder1, activeInSample = bladder_sigs, binSize = 300)
base::saveRDS(bladder1_traj, file="bladder1_traj.rds")

bladder2_traj <- SpaceTrack(bladder2, activeInSample = bladder_sigs, binSize = 300)
base::saveRDS(bladder2_traj, file="bladder2_traj.rds")

bladder3_traj <- SpaceTrack(bladder3, activeInSample = bladder_sigs, binSize = 300)
base::saveRDS(bladder3_traj, file="bladder3_traj.rds")

bladder4_traj <- SpaceTrack(bladder4, activeInSample = bladder_sigs, binSize = 300)
base::saveRDS(bladder4_traj, file="bladder4_traj.rds")

bladder5_traj <- SpaceTrack(bladder5, activeInSample = bladder_sigs, binSize = 300)
base::saveRDS(bladder5_traj, file="bladder5_traj.rds")

# Get chromosomal activity means for each sample
bladder1_means <- chromosomalActivities(bladder1_traj, bladder1_traj[[1]], bladder_sigs, sample = 1)

bladder2_means <- chromosomalActivities(bladder2_traj, bladder2_traj[[1]], bladder_sigs, sample = 2)

bladder3_means <- chromosomalActivities(bladder3_traj, bladder3_traj[[1]], bladder_sigs, sample = 3)

bladder4_means <- chromosomalActivities(bladder4_traj, bladder4_traj[[1]], bladder_sigs, sample = 4)

bladder5_means <- chromosomalActivities(bladder5_traj, bladder5_traj[[1]], bladder_sigs, sample = 5)

# Plot individual trajectories
plotSpaceTrajectory(trajectory = bladder1_traj) + labs(title = "Bladder-TCC | n = 1 | Binsize = 300 mutations | Total mutations = 13,459")
plotSpaceTrajectory(trajectory = bladder2_traj) + labs(title = "Bladder-TCC | n = 1 | Binsize = 300 mutations | Total mutations = 9,948")
plotSpaceTrajectory(trajectory = bladder3_traj) + labs(title = "Bladder-TCC | n = 1 | Binsize = 300 mutations | Total mutations = 19,931")
plotSpaceTrajectory(trajectory = bladder4_traj) + labs(title = "Bladder-TCC | n = 1 | Binsize = 300 mutations | Total mutations = 19,090")
plotSpaceTrajectory(trajectory = bladder5_traj) + labs(title = "Bladder-TCC | n = 1 | Binsize = 300 mutations | Total mutations = 12,673")

# Plot average chromosomal activities
activity_means <- rbind(bladder1_means, bladder2_means, bladder3_means, bladder4_means, bladder5_means)
activity_long <- tidyr::pivot_longer(activity_means, cols = all_of(bladder_sigs), names_to = "Signatures", values_to = "mean_activity")
activity_long$sample <- as.factor(activity_long$sample)
activity_means <- activity_long %>%
  dplyr::group_by(chromosome, Signatures) %>%
  dplyr::summarize(mean_activity = mean(mean_activity))

ggplot2::ggplot(data = activity_means, aes(x = chromosome, y = mean_activity*100, color = Signatures)) +
  ggplot2::geom_line(alpha = 1) +
  ggplot2::theme_bw() +
  ggplot2::scale_x_continuous(breaks = c(1:24), labels = as.character(c(c(1:22), "X", "Y"))) +
  ggplot2::labs(x = "Chromosome", y = "Average Chromosomal Signature Activity (%)", title = "Bladder-TCC | n = 5 | Total mutations = 75,101")

# Calculate trajectory on pooled samples

# initialize counts dataframe
bladder_master <- readr::read_csv(as.character(paste(archivePath, "/", bladder_sample[1], ".MBcounts.csv", sep = "")),
                          col_types = readr::cols("seqnames" = "c",
                                                  "strand" = "c",
                                                  .default = "d"))

# add counts from remaining samples into master file and delete remaining files from memory
for (i in 2:length(bladder_sample)) {
  temp <- readr::read_csv(as.character(paste(archivePath, "/", bladder_sample[i], ".MBcounts.csv", sep = "")),
                          col_types = readr::cols(seqnames = "c",
                                                  strand = "c",
                                                  .default = "d"))
  bladder_master[, 6:101] <- bladder_master[, 6:101] + temp[, 6:101]
  base::remove(temp)
}

# make data type consistent in seqnames column
# split seqnames column into start chrom and end chrom
# remove unnecessary variables
bladder_master <- bladder_master %>%
  dplyr::mutate(seqnames = dplyr::case_when(seqnames == "X" ~ 23,
                                            seqnames == "Y" ~ 24,
                                            TRUE ~ as.numeric(seqnames)),
                start_chrom = seqnames,
                end_chrom = seqnames) %>%
  dplyr::select(-seqnames, -strand) %>%
  base::subset(select = c(100,1,101,2:99))



# Plot pooled trajectory
plotSpaceTrajectory(pooled_bladder_traj) + labs(title = "Bladder-TCC | n = 5 | Binsize = 300 mutations | Total mutations = 75,101")


# ColoRect-AdenoCA

archivePath <- "~/Desktop/Caitlin/archive"
typesPath <- "~/Desktop/Caitlin/data/pcawg_cancer_types.csv"
cancerType <- "ColoRect-AdenoCA"
colorect_sigs <- c("SBS1", "SBS5", "SBS8", "SBS40")

types <- readr::read_csv(typesPath,
                         col_names = c('type', 'guid'))
types <- types[2:nrow(types), ]

# list of sample filenames for desired cancer type
colorect_files <- c(types$guid[types$type == cancerType])
set.seed(25)
colorect_sample <- base::sample(colorect_files, size = length(colorect_files), replace = FALSE)

# Find mutation counts in each sample
for (i in 1:length(colorect_sample)) {
  temp <- readFormat(path = paste("~/Desktop/Caitlin/archive/", colorect_sample[i], ".MBcounts.csv", sep=''))
  print(paste(i, sum(colSums(temp[,6:101])), sep = ": "))
}

small_samples <- colorect_sample[c(7,13,57,39,42)]
large_samples <- colorect_sample[c(30,60,53,28,11)]

# Pool samples
colorect_large <- readr::read_csv(as.character(paste(archivePath, "/", large_samples[1], ".MBcounts.csv", sep = "")),
                                 col_types = readr::cols("seqnames" = "c",
                                                         "strand" = "c",
                                                         .default = "d"))

# add counts from remaining samples into master file and delete remaining files from memory
for (i in 2:length(large_samples)) {
  temp <- readr::read_csv(as.character(paste(archivePath, "/", large_samples[i], ".MBcounts.csv", sep = "")),
                          col_types = readr::cols(seqnames = "c",
                                                  strand = "c",
                                                  .default = "d"))
  colorect_large[, 6:101] <- colorect_large[, 6:101] + temp[, 6:101]
  base::remove(temp)
}

# make data type consistent in seqnames column
# split seqnames column into start chrom and end chrom
# remove unnecessary variables
colorect_large <- colorect_large %>%
  dplyr::mutate(seqnames = dplyr::case_when(seqnames == "X" ~ 23,
                                            seqnames == "Y" ~ 24,
                                            TRUE ~ as.numeric(seqnames)),
                start_chrom = seqnames,
                end_chrom = seqnames) %>%
  dplyr::select(-seqnames, -strand) %>%
  base::subset(select = c(100,1,101,2:99))

# Run SpaceTrack on small samples-- chromosomal level
colorect_large_traj <- SpaceTrack(path = colorect_large, activeInSample = colorect_sigs, binSize = 300, chr_level=TRUE, parallelize = TRUE)
# Plot small samples
plotSpaceTrajectory(colorect_large_traj, chr_level = TRUE, cutoff=0) + labs(title = "ColoRect-AdenoCA | n = 5 | Binsize = 300 mutations | Total mutations = 3,649,731")



# Read in individual dataframes
colorect1 <- readFormat(path = paste("~/Desktop/Caitlin/archive/", colorect_sample[1], ".MBcounts.csv", sep=''))
sum(colSums(colorect1[,6:101]))

colorect2 <- readFormat(path = paste("~/Desktop/Caitlin/archive/", colorect_sample[2], ".MBcounts.csv", sep=''))
sum(colSums(colorect2[,6:101]))

colorect3 <- readFormat(path = paste("~/Desktop/Caitlin/archive/", colorect_sample[3], ".MBcounts.csv", sep=''))
sum(colSums(colorect3[,6:101]))

colorect4 <- readFormat(path = paste("~/Desktop/Caitlin/archive/", colorect_sample[4], ".MBcounts.csv", sep=''))
sum(colSums(colorect4[,6:101]))

colorect5 <- readFormat(path = paste("~/Desktop/Caitlin/archive/", colorect_sample[5], ".MBcounts.csv", sep=''))
sum(colSums(colorect5[,6:101]))

# Calculate trajectories for each sample
colorect1_traj <- SpaceTrack(colorect1, activeInSample = colorect_sigs, binSize = 300)
base::saveRDS(colorect1_traj, file="colorect1_traj.rds")

colorect2_traj <- SpaceTrack(colorect2, activeInSample = colorect_sigs, binSize = 300)
base::saveRDS(colorect2_traj, file="colorect2_traj.rds")

colorect3_traj <- SpaceTrack(colorect3, activeInSample = colorect_sigs, binSize = 300)
base::saveRDS(colorect3_traj, file="colorect3_traj.rds")

colorect4_traj <- SpaceTrack(colorect4, activeInSample = colorect_sigs, binSize = 300)
base::saveRDS(colorect4_traj, file="colorect4_traj.rds")

colorect5_traj <- SpaceTrack(colorect5, activeInSample = colorect_sigs, binSize = 300)
base::saveRDS(colorect5_traj, file="colorect5_traj.rds")

# Get chromosomal activity means for each sample
colorect1_means <- chromosomalActivities(colorect1_traj, colorect1_traj[[1]], colorect_sigs, sample = 1)

colorect2_means <- chromosomalActivities(colorect2_traj, colorect2_traj[[1]], colorect_sigs, sample = 2)

colorect3_means <- chromosomalActivities(colorect3_traj, colorect3_traj[[1]], colorect_sigs, sample = 3)

colorect4_means <- chromosomalActivities(colorect4_traj, colorect4_traj[[1]], colorect_sigs, sample = 4)

colorect5_means <- chromosomalActivities(colorect5_traj, colorect5_traj[[1]], colorect_sigs, sample = 5)

# Plot individual trajectories
plotSpaceTrajectory(trajectory = colorect1_traj) + labs(title = "ColoRect-AdenoCA | n = 1 | Binsize = 300 mutations | Total mutations = 9,435")
plotSpaceTrajectory(trajectory = colorect2_traj) + labs(title = "ColoRect-AdenoCA | n = 1 | Binsize = 300 mutations | Total mutations = 114,626")
plotSpaceTrajectory(trajectory = colorect3_traj) + labs(title = "ColoRect-AdenoCA | n = 1 | Binsize = 300 mutations | Total mutations = 19,421")
plotSpaceTrajectory(trajectory = colorect4_traj) + labs(title = "ColoRect-AdenoCA | n = 1 | Binsize = 300 mutations | Total mutations = 30,359")
plotSpaceTrajectory(trajectory = colorect5_traj) + labs(title = "ColoRect-AdenoCA | n = 1 | Binsize = 300 mutations | Total mutations = 10,904")

# Plot average chromosomal activities
activity_means <- rbind(colorect1_means, colorect2_means, colorect3_means, colorect4_means, colorect5_means)
activity_long <- tidyr::pivot_longer(activity_means, cols = colorect_sigs, names_to = "Signatures", values_to = "mean_activity")
activity_long$sample <- as.factor(activity_long$sample)
activity_means <- activity_long %>%
  dplyr::group_by(chromosome, Signatures) %>%
  dplyr::summarize(mean_activity = mean(mean_activity))

ggplot2::ggplot(data = activity_means, aes(x = chromosome, y = mean_activity*100, color = Signatures)) +
  ggplot2::geom_line(alpha = 1) +
  ggplot2::theme_bw() +
  ggplot2::scale_x_continuous(breaks = c(1:24), labels = as.character(c(c(1:22), "X", "Y"))) +
  ggplot2::labs(x = "Chromosome", y = "Average Chromosomal Signature Activity (%)", title = "ColoRect-AdenoCA | n = 5 | Total mutations = 184,745")

# Calculate trajectory on pooled samples

# initialize counts dataframe
colorect_master <- readr::read_csv(as.character(paste(archivePath, "/", colorect_sample[1], ".MBcounts.csv", sep = "")),
                                  col_types = readr::cols("seqnames" = "c",
                                                          "strand" = "c",
                                                          .default = "d"))

# add counts from remaining samples into master file and delete remaining files from memory
for (i in 2:length(colorect_sample)) {
  temp <- readr::read_csv(as.character(paste(archivePath, "/", colorect_sample[i], ".MBcounts.csv", sep = "")),
                          col_types = readr::cols(seqnames = "c",
                                                  strand = "c",
                                                  .default = "d"))
  colorect_master[, 6:101] <- colorect_master[, 6:101] + temp[, 6:101]
  base::remove(temp)
}

# make data type consistent in seqnames column
# split seqnames column into start chrom and end chrom
# remove unnecessary variables
colorect_master <- colorect_master %>%
  dplyr::mutate(seqnames = dplyr::case_when(seqnames == "X" ~ 23,
                                            seqnames == "Y" ~ 24,
                                            TRUE ~ as.numeric(seqnames)),
                start_chrom = seqnames,
                end_chrom = seqnames) %>%
  dplyr::select(-seqnames, -strand) %>%
  base::subset(select = c(100,1,101,2:99))

# Determine trajectory
pooled_bladder_traj <- SpaceTrack(path = bladder_master, activeInSample = bladder_sigs, binSize = 300)
base::saveRDS(pooled_bladder_traj, file = "pooled_bladder_traj.rds")
# Determine trajectory
pooled_colorect_traj <- SpaceTrack(path = colorect_master, activeInSample = colorect_sigs, binSize = 300)
base::saveRDS(pooled_colorect_traj, file = "pooled_colorect_traj.rds")

# Plot pooled trajectory
plotSpaceTrajectory(pooled_colorect_traj) + labs(title = "ColoRect-AdenoCA | n = 5 | Binsize = 300 mutations | Total mutations = 184,745")
