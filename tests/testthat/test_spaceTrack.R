
myCluster <- makeCluster(2, type = "FORK")
registerDoParallel(myCluster)
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
lymph_traj <- combineTraj(lymph_traj)
lymph_traj[['changepoints']] <- sort(lymph_traj[['changepoints']])
lymph_traj[['binData']] <- lymph_traj[['binData']] %>%
  dplyr::arrange(dplyr::desc(bin))

plotTrajectory(lymph_traj, linearX = T) + labs(title = "Lymph-CLL (n = 94), Binsize = 200 mutations")

lymph_counts <- binningNmut(lymph_master, 200)



# # iteratively collapse rows of dataframe into bins with ~binSize total mutations per row
# for (i in 1:nrow(counts)) {
#   row_sums <- base::sum(base::colSums(counts[i,6:101]))
#   if (i==nrow(counts)) {break}
#   while (row_sums < binSize) {
#     counts[i, 5:101] <- as.list(base::colSums(counts[c(i,i+1),c(5:101)]))
#     counts[i,3:4] <- counts[i+1,3:4]
#     row_sums <- base::sum(base::colSums(counts[i,6:101]))
#     counts <- counts[-c(i+1), ]
#     if (i == nrow(counts)-1) {break}
#   }
# }

