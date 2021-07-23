# Initialize data - CNS-PiloAstro cancer type
pilo_sigs <- c("SBS1", "SBS5", "SBS19", "SBS23", "SBS40")
pilo_master <- readr::read_csv('~/Desktop/Caitlin/data/CNS-PiloAstro_pooled.csv')

# CNS-Oligo
oligo_sigs <- c("SBS1", "SBS5", "SBS8", "SBS40")
oligo_master <- readr::read_csv('~/Desktop/Caitlin/data/CNS-Oligo_pooled.csv')

# Bootstrap shuffling
myCluster <- parallel::makeCluster(15, type="FORK")
doParallel::registerDoParallel(myCluster)
`%dopar%` <- foreach::`%dopar%`
parallel::stopCluster(myCluster)

count <- 0
indices <- c(1,5,9,13,17,21,25,29,33,37,41,45,49,53,57,61,65,69,73,77)
for (i in indices) {
  count <- count +1
  if (pilo_shuffle30[[i]][1,5] != pilo_shuffle30[[i]][1,6]) {
    print(paste('Shuffle bootstrap #', count, sep = ""))
    print(pilo_shuffle30[[i]][,5:6])
  }
}

# chrom_bins <- oligo_master %>%
#   dplyr::group_by(start_chrom) %>%
#   dplyr::summarize(mb_bins = dplyr::n())
# chrom_bins <- chrom_bins[1:22,]
# chrom_bins <- rbind(chrom_bins, c("MT", 3), c("X", 156), c("Y", 60))
# # Epigenetic state
# chrom_pos <- data.frame(chrom = EPI@seqnames@values, length = EPI@seqnames@lengths) %>%
#   dplyr::mutate(start = 1, end = 1, mb_bins = 1)
# chrom_pos$end[1] <- chrom_pos$length[1]
# chrom_pos$mb_bins[1] <- chrom_bins$mb_bins[chrom_bins$start_chrom==chrom_pos$chrom[1]]
# for (i in 2:nrow(chrom_pos)) {
#   chrom_pos$start[i] <- chrom_pos$end[i-1] + 1
#   chrom_pos$end[i] <- chrom_pos$end[i-1] + chrom_pos$length[i]
#   chrom_pos$mb_bins[i] <- chrom_bins$mb_bins[chrom_bins$start_chrom==chrom_pos$chrom[i]]
# }
#
# start_chrom <- rep(EPI@seqinfo@seqnames, times = EPI@seqnames@lengths)
# EPI$mb_bin <- 1
# EPI$width <- EPI@ranges@width
# for (i in 1:nrow(chrom_pos)) {
#   epi_chrom <- EPI@ranges@start[chrom_pos$start[i]:chrom_pos$end[i]]
#   prev_index <- 1
#   for (i in 1:chrom_pos$mb_bins[i]-1) {
#     epi_index <- base::which.min(abs(c(epi_chrom)-(i*1000000)))
#     EPI@elementMetadata@listData$mb_bin[prev_index:epi_index] <- i
#     prev_index <- epi_index
#   }
#   EPI@elementMetadata@listData$mb_bin[prev_index:length(epi_chrom)] <- chrom_pos$mb_bins[i]
#  }
#
# base::saveRDS(EPI, file="EPI.rds")
#
# chrom_states <- data.frame(chromosome = start_chrom, mb_bin = EPI$mb_bin, state = EPI$state, width = EPI$width)
# chrom_states <- chrom_states %>%
#   dplyr::group_by_at(vars(chromosome, mb_bin, state)) %>%
#   dplyr::summarize(width = sum(width))


oligo_shuffle30 <- foreach::foreach(j = c(1:30), .combine="c") %dopar% bootstrapShuffle(master = oligo_master, binSize = 300, activeInSample = oligo_sigs, i = j)

# Random ordering of chromosomes

chroms <- c(1:23)
order <- base::sample(chroms, size = length(chroms), replace=FALSE)

pilo_chroms <- binByChromRandom(pilo_master, 200)

pilo_counts <- data.frame()
for (i in order) {
  pilo_counts <- rbind(pilo_counts, pilo_chroms[[i]])
}
pilo_counts$bin <- rep(1:nrow(pilo_counts))

# Run TrackSig
pilo_randomtraj1 <- TrackSig(df = pilo_counts, activeInSample = pilo_sigs, binSize = 200)
base::saveRDS(pilo_randomtraj1, file = "pilo_randomtraj1.rds")

# Plot results
plotSpaceTrajectoryRandom(oligo_shuffle30) + labs(title = "CNS-Oligo | n = 18 | Binsize = 300 mutations | Total mutations = 48,836 | 30 shuffle bootstraps")
# Repeat ~2 more times


# Ordering chromosomes by length
length_order <- c(1:7, 23, 8:22)

pilo_counts <- data.frame()
for (i in length_order) {
  pilo_counts <- rbind(pilo_counts, pilo_chroms[[i]])
}
pilo_counts$bin <- rep(1:nrow(pilo_counts))

# Run TrackSig with bootstrapping
myCluster <- parallel::makeCluster(15, type="FORK")
doParallel::registerDoParallel(myCluster)
`%dopar%` <- foreach::`%dopar%`


pilo_lengthtraj <- foreach::foreach(j = c(1:5), .combine='c') %dopar% TrackSig(df = pilo_counts, activeInSample = pilo_sigs, binSize = 200)
base::saveRDS(pilo_randomtraj1, file = "pilo_randomtraj1.rds")

# Plot results
plotSpaceTrajectoryRandom(pilo_randomtraj1) + labs(title = "CNS-PiloAstro | n = 89 | Binsize = 200 mutations | Total mutations = 22,137 ")





