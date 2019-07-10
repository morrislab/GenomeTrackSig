compare_simulations <- function(resultsDir, dataDir, outDir){

  dir.create(outDir, showWarnings = F, recursive = T)

  # ===========================================
  # Compare TrackSig exposures to ground truth
  simulations <- list.files(dataDir)
  sel <- grep(x = simulations, "^Simulation")
  simulations <- simulations[sel]

  list[res, gt_exp_l, estim_exp_l] <- TrackSig:::compare_simulation_results(
      simulations, ground_truth_dir = dataDir,
      method_results_dir = resultsDir,
      res_file_name = paste(outDir, "TrackSig_simulation_results.txt", sep = "/"))

  res_TrackSig <- res
  TrackSig:::plot_kl_results(res, paste(outDir, "TrackSig", sep = "/"))

  print("Percentage of samples with KL larger than 0.05")
  mean(res$kl > 0.1)

#  print("Average KL per type: TrackSig")
  #avg_kl_per_type_TrackSig <- TrackSig:::get_results_per_sim_type(res, paste(outDir, "TrackSig", sep = "/"))
  #print(avg_kl_per_type_TrackSig)

  # ===========================================
  # Investigate samples with large KL

  simulations_large_kl <- res[res$kl > 0.05,"sim"]
  sigs_in_bad_sims <- TrackSig:::get_sig_names_from_list(gt_exp_l[simulations_large_kl])
  table(unlist(sigs_in_bad_sims))
  # Signature 7 is present in ALL the simulations where TrackSig results differ from ground truth

  # has_sig7 <- list()
  # for (sim in names(gt_exp_l)) {
  #     has_sig7[[sim]] <- "SBS7" %in% colnames(gt_exp_l[[sim]][3:6])
  # }
  # sim_order <- res_TrackSig[,1]
  # has_sig7 = has_sig7[sim_order]
  # has_sig7 <- as.numeric(data.frame(has_sig7))
  # has_sig7[has_sig7 == 0] <- "darkgreen"
  # has_sig7[has_sig7 == 1] <- "magenta"

  pdf(paste(outDir, "TrackSig_simulation_results_KL_vs_mean.pdf", sep = "/"), width = 5, height=5)
  plot(res_TrackSig$kl, res_TrackSig$abs_diff_mean, main="TrackSig versus Truth",
     xlab="KL", ylab="mean abs diff",  xlim=c(0, 1), ylim=c(0, 1)) #col=has_sig7)
  #legend('topleft', legend = c("No SBS7", "Has SBS7"), col = c("darkgreen", "magenta"), cex = 0.8, pch = 1)
  dev.off()

  pdf(paste(outDir, "TrackSig_simulation_results_KL_vs_max.pdf", sep = "/"), width = 5, height=5)
  plot(res_TrackSig$kl, res_TrackSig$abs_diff_max, main="TrackSig versus Truth",
     xlab="KL", ylab="max abs diff", xlim=c(0, 1), ylim=c(0, 1)) # col=has_sig7)
  #legend('topleft', legend = c("No SBS7", "Has SBS7"), col = c("darkgreen", "magenta"), cex = 0.8, pch = 1)
  dev.off()


#  # ===========================================
#  # Compare number of change-points
  # simulations_depth100 <- simulations[grepl("depth100$", simulations)]
  # simulations_depth1000 <- simulations[grepl("depth1000$", simulations)]
  # Compare change-points
  cp_comparison <- TrackSig:::compare_changepoints(simulations,
    ground_truth_dir = dataDir,
    tracksig_results_dir = paste0(resultsDir, "/SIMULATED/"),
    sciclone_results_dir = NULL,
    res_file_name = paste(outDir, "cp_comparison.txt", sep = "/"),
    change_at_cp_threshold = 0.05)

  # choose depth 30 only
  # idx <-  grepl(paste0("depth30$"), sapply(cp_comparison[,1],toString))
  # cp_comparison <- cp_comparison[idx, ]

  {
#  print("TrackSig agrees with sciclone")
#  print(mean(cp_comparison$cp_tracksig_adjusted == cp_comparison$cp_SCDS_results))

  print("TrackSig agrees with GT")
  print(mean(cp_comparison$n_gt_created_cp == cp_comparison$cp_tracksig_adjusted))
#  print("SciClone agrees with GT")
#  print(mean(cp_comparison$n_gt_created_cp == cp_comparison$cp_SCDS_results))

  print("TrackSig overestimates # change-points")
  print(mean(cp_comparison$n_gt_created_cp < cp_comparison$cp_tracksig_adjusted))
  print("TrackSig underestimates # change-points")
  print(mean(cp_comparison$n_gt_created_cp > cp_comparison$cp_tracksig_adjusted))

#  print("SciClone overestimates # change-points")
#  print(mean(cp_comparison$n_gt_created_cp < cp_comparison$cp_SCDS_results))
#  print("SciClone underestimates # change-points")
#  print(mean(cp_comparison$n_gt_created_cp > cp_comparison$cp_SCDS_results))

  print("Mean difference between ground truth and method")
  print(mean(cp_comparison$n_gt_created_cp - cp_comparison$cp_tracksig_adjusted))
#  print(mean(cp_comparison$n_gt_created_cp - cp_comparison$cp_SCDS_results))

  print("Mean abs difference between ground truth and method")
  print(mean(abs(cp_comparison$n_gt_created_cp - cp_comparison$cp_tracksig_adjusted)))
#  print(mean(abs(cp_comparison$n_gt_created_cp - cp_comparison$cp_SCDS_results)))
  }



  sim_type <- "two_clusters"
  d_type <- "depth100"
  idx <- grepl(sim_type, sapply(cp_comparison[,1],toString))
  idx <- idx &  grepl(paste0(d_type,"$"), sapply(cp_comparison[,1],toString))


  print(cp_comparison[idx,]$n_gt_created_cp)
  print(cp_comparison[idx,]$cp_tracksig_adjusted)
  print(cp_comparison[idx,]$cp_tracksig)
  print(cp_comparison[idx,]$cp_SCDS_results)

  print(mean(cp_comparison[idx,]$n_gt_created_cp > cp_comparison[idx,]$cp_tracksig_adjusted))
  print(mean(cp_comparison[idx,]$n_gt_created_cp < cp_comparison[idx,]$cp_tracksig_adjusted))
  cp_comparison[idx,][cp_comparison[idx,]$cp_SCDS_results > 1,]



  return(list(res, cp_comparison))
}


# [END]

