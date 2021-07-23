library(ggplot2)
# Make masters

uterine_sigs <- c("SBS1", "SBS2.13", "SBS3", "SBS5", "SBS6", 'SBS10', "SBS14", "SBS15", "SBS26", 'SBS28', "SBS40", 'SBS44')
colon_sigs <- c("SBS1", "SBS5", "SBS8", "SBS40")










# uterine10 <- chromosomalActivities(uterine_trajectories[37:40], uterine_trajectories[[37]], uterine_sigs, 10)
# activity_means <- rbind(uterine1, uterine2, uterine3, uterine4, uterine5, uterine6, uterine7, uterine8, uterine9, uterine10)
# activity_long <- tidyr::pivot_longer(activity_means, cols = uterine_sigs, names_to = "Signatures", values_to = "mean_activity")
# activity_long$sample <- as.factor(activity_long$sample)
# activity_means <- activity_long %>%
#   dplyr::group_by(chromosome, Signatures) %>%
#   dplyr::summarize(mean_activity = mean(mean_activity))
#
# ggplot2::ggplot(data = activity_long, aes(x = chromosome, y = mean_activity*100, color = Signatures, shape = sample)) +
#   ggplot2::geom_point() +
#   ggplot2::geom_line(alpha = 1) +
#   ggplot2::theme_bw() +
#   ggplot2::scale_x_continuous(breaks = c(1:24), labels = as.character(c(c(1:22), "X", "Y"))) +
#   ggplot2::labs(x = "Chromosome", y = "Average Chromosomal Signature Activity (%)", shape = "Sample ID", title = "Uterine-AdenoCA | n = 10 | Total mutations = 89,375")
#


# # Panc-Endocrine
# endocrine_sigs <- c("SBS1", "SBS2.13", "SBS3", "SBS5", "SBS6", "SBS8", "SBS9", "SBS11", "SBS26", "SBS30", "SBS36", "SBS39")
# endocrine_traj <- SpaceTrack("~/Desktop/Caitlin/data/Panc-Endocrine_pooled.csv",
#                              bootstrapSamples = 30, activeInSample = endocrine_sigs, binSize = 200, parallelize = TRUE)
# base::saveRDS(endocrine_traj, file = "endocrine_30.rds")
#
# # Bone-Osteosarc
# osteo_sigs <- c("SBS1", "SBS2.13", "SBS3", "SBS5", "SBS8", "SBS17", "SBS30", "SBS40")
# osteo_traj <- SpaceTrack("~/Desktop/Caitlin/data/Bone-Osteosarc_pooled.csv",
#                          bootstrapSamples = 30, activeInSample = osteo_sigs, binSize = 200, parallelize = TRUE)
# base::saveRDS(osteo_traj, file = "osteo30.rds")
#
