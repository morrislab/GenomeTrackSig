compare_simulations <- function(resultsDir, dataDir, outDir){

  list <- structure(NA,class="result")
  "[<-.result" <- function(x,...,value) {
    args <- as.list(match.call())
    args <- args[-c(1:2,length(args))]
    length(value) <- length(args)
    for(i in seq(along=args)) {
      a <- args[[i]]
      if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
    }
    x
  }


  # ===========================================
  # Compare TrackSig exposures to ground truth
  simulations <- list.files(dataDir)
  sel <- grep(x = simulations, "^Simulation")
  simulations <- simulations[sel]

  tracksig_dir = resultsDir

  list[res, gt_exp_l, estim_exp_l] <- TrackSig:::compare_simulation_results(
      simulations, ground_truth_dir = dataDir,
      method_results_dir = paste0(tracksig_dir, "/SIMULATED/"),
      res_file_name = paste(outDir, "TrackSig_simulation_results.txt", sep = "/"))

  res_TrackSig <- res
  TrackSig:::plot_kl_results(res, paste(outDir, "TrackSig", sep = "/"))

  print("Percentage of samples with KL larger than 0.05")
  mean(res$kl > 0.1)

  print("Average KL per type: TrackSig")
  avg_kl_per_type_TrackSig <- TrackSig:::get_results_per_sim_type(res, paste(outDir, "TrackSig", sep = "/"))
  print(avg_kl_per_type_TrackSig)

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




  # ===========================================
  # Compare SciClone exposures to ground truth
  simulations <- list.files(dataDir)
  sel <- grep(x = simulations, "^Simulation")
  simulations <- simulations[sel]

  sciclone_dir <- "SCDS_results/"
  list[res, gt_exp_l, estim_exp_l] <- TrackSig:::compare_simulation_results(simulations,
      ground_truth_dir = dataDir,
      method_results_dir = paste0(sciclone_dir, "/SIMULATED/"),
      res_file_name = paste(outDir, "sciclone_simulation_results.txt", sep = "/"))

  res_SciClone <- res
  estim_exp_l_sciclone <- estim_exp_l
  TrackSig:::plot_kl_results(res, paste(outDir, "SciClone", sep = "/"))

  print("sciclone: Percentage of samples with KL larger than 0.05")
  mean(res$kl > 0.05)

  print("Average KL per type: SciClone")
  avg_kl_per_type_SciClone <- TrackSig:::get_results_per_sim_type(res_SciClone, paste(outDir, "SciClone", sep = "/"))
  print(avg_kl_per_type_SciClone)




  # has_sig7 <- list()
  # for (sim in names(gt_exp_l)) {
  #     has_sig7[[sim]] <- "SBS7" %in% colnames(gt_exp_l[[sim]][3:6])
  # }

  # sim_order <- res_SciClone[,1]
  # has_sig7 = has_sig7[sim_order]

  # has_sig7 <- as.numeric(data.frame(has_sig7))
  # has_sig7[has_sig7 == 0] <- "darkgreen"
  # has_sig7[has_sig7 == 1] <- "magenta"

  pdf(paste(outDir, "SciClone_simulation_results_KL_vs_mean.pdf", sep = "/"), width = 5, height=5)
  plot(res_SciClone$kl, res_SciClone$abs_diff_mean, main="SciClone versus Truth",
     xlab="KL", ylab="mean abs diff", xlim=c(0, 1), ylim=c(0, 1)) #col=has_sig7)
  #legend('topleft', legend = c("No SBS7", "Has SBS7"), col = c("darkgreen", "magenta"), cex = 0.8, pch = 1)
  dev.off()

  pdf(paste(outDir, "SciClone_simulation_results_KL_vs_max.pdf", sep = "/"), width = 5, height=5)
  plot(res_SciClone$kl, res_SciClone$abs_diff_max, main="SciClone versus Truth",
     xlab="KL", ylab="max abs diff", xlim=c(0, 1), ylim=c(0, 1)) #col=has_sig7)
  #legend('topleft', legend = c("No SBS7", "Has SBS7"), col = c("darkgreen", "magenta"), cex = 0.8, pch = 1)
  dev.off()






  # ===========================================
  # Compare TrackSig and SciClone
  sim_order <- intersect(res_SciClone[,1], res_TrackSig[,1])
  print("TrackSig has lower KL in this % of simulations:")
  mean(res_TrackSig[sim_order,]$kl < res_SciClone[sim_order,]$kl)

  print("Mean KL between the method and ground truth")
  mean(res_TrackSig$kl)
  mean(res_SciClone$kl)

  print("mean and median KL diff between two methods")
  mean(res_SciClone[sim_order,]$kl - res_TrackSig[sim_order,]$kl)
  median(res_SciClone[sim_order,]$kl - res_TrackSig[sim_order,]$kl)

  print("TrackSig has lower mean abs diff in this % of simulations:")
  mean(res_TrackSig[sim_order,]$abs_diff_mean < res_SciClone[sim_order,]$abs_diff_mean)

  print("TrackSig has lower max abs diff in this % of simulations:")
  mean(res_TrackSig[sim_order,]$abs_diff_max < res_SciClone[sim_order,]$abs_diff_max)

  print("percentage of correct reconstructions")
  mean(res_TrackSig$kl < 0.05)
  mean(res_SciClone$kl < 0.05)

  mean(res_TrackSig$kl < 0.1)
  mean(res_SciClone$kl < 0.1)


  correct_tracksig = res_TrackSig[sim_order,]$kl < 0.05

  correct_tracksig <- as.numeric(correct_tracksig)
  correct_tracksig[correct_tracksig == 0] <- "red"
  correct_tracksig[correct_tracksig == 1] <- "darkgreen"

  pdf(paste(outDir, "sciclone_results_coloured_TrackSig_mean_diff.pdf", sep = "/"), width = 5, height=5)
  plot(res_SciClone[sim_order,]$kl, res_SciClone[sim_order,]$abs_diff_mean,
    main="Sciclone versus Truth",
     xlab="KL", ylab="mean abs diff", col=correct_tracksig,
     xlim=c(0, 1), ylim=c(0, 1))
  legend('topleft', legend = c("TrackSig: KL > 0.05", "TrackSig: KL < 0.05"), col = c("red", "darkgreen"), cex = 0.8, pch = 1)
  dev.off()

  pdf(paste(outDir, "sciclone_results_coloured_TrackSig_max_diff.pdf", sep = "/"), width = 5, height=5)
  plot(res_SciClone[sim_order,]$kl, res_SciClone[sim_order,]$abs_diff_max,
    main="Sciclone versus Truth",
     xlab="KL", ylab="max abs diff", col=correct_tracksig,
     xlim=c(0, 1), ylim=c(0, 1))
  legend('topleft', legend = c("TrackSig: KL > 0.05", "TrackSig: KL < 0.05"), col = c("red", "darkgreen"), cex = 0.8, pch = 1)
  dev.off()



  # Comparison of the KL between TrackSig and SciClone
  sim_order <- intersect(res_SciClone[,1], res_TrackSig[,1])
  pdf(paste(outDir, "TrackSig_vs_SciCLone_simulation_results_KL.pdf", sep = "/"), width = 5, height=5)
  plot(res_SciClone[sim_order,]$kl, res_TrackSig[sim_order,]$kl,
     xlab="SciClone KL", ylab="TrackSig KL", xlim=c(0, 0.5), ylim=c(0, 0.5))
  dev.off()

  # for depth 30
  sim_order <- intersect(res_SciClone[res_SciClone$depth == 30,1], res_TrackSig[res_TrackSig$depth == 30,1])
  pdf(paste(outDir, "TrackSig_vs_SciCLone_simulation_results_KL_depth_30.pdf", sep = "/"), width = 5, height=5)
  plot(res_SciClone[sim_order,]$kl, res_TrackSig[sim_order,]$kl,
     xlab="SciClone KL", ylab="TrackSig KL", xlim=c(0, 0.5), ylim=c(0, 0.5))
  dev.off()


  # ===========================================
  # Compare number of change-points
  # simulations_depth100 <- simulations[grepl("depth100$", simulations)]
  # simulations_depth1000 <- simulations[grepl("depth1000$", simulations)]
  # Compare change-points
  cp_comparison <- TrackSig:::compare_changepoints(simulations,
    ground_truth_dir = dataDir,
    tracksig_results_dir = paste0(tracksig_dir, "/SIMULATED/"),
    sciclone_results_dir = paste0(sciclone_dir, "/SIMULATED/"),
    res_file_name = paste(outDir, "cp_comparison.txt", sep = "/"),
    change_at_cp_threshold = 0.05)

  # choose depth 30 only
  # idx <-  grepl(paste0("depth30$"), sapply(cp_comparison[,1],toString))
  # cp_comparison <- cp_comparison[idx, ]

  {
  print("TrackSig agrees with sciclone")
  print(mean(cp_comparison$cp_tracksig_adjusted == cp_comparison$cp_SCDS_results))

  print("TrackSig agrees with GT")
  print(mean(cp_comparison$n_gt_created_cp == cp_comparison$cp_tracksig_adjusted))
  print("SciClone agrees with GT")
  print(mean(cp_comparison$n_gt_created_cp == cp_comparison$cp_SCDS_results))

  print("TrackSig overestimates # change-points")
  print(mean(cp_comparison$n_gt_created_cp < cp_comparison$cp_tracksig_adjusted))
  print("TrackSig underestimates # change-points")
  print(mean(cp_comparison$n_gt_created_cp > cp_comparison$cp_tracksig_adjusted))

  print("SciClone overestimates # change-points")
  print(mean(cp_comparison$n_gt_created_cp < cp_comparison$cp_SCDS_results))
  print("SciClone underestimates # change-points")
  print(mean(cp_comparison$n_gt_created_cp > cp_comparison$cp_SCDS_results))

  print("Mean difference between ground truth and method")
  print(mean(cp_comparison$n_gt_created_cp - cp_comparison$cp_tracksig_adjusted))
  print(mean(cp_comparison$n_gt_created_cp - cp_comparison$cp_SCDS_results))

  print("Mean abs difference between ground truth and method")
  print(mean(abs(cp_comparison$n_gt_created_cp - cp_comparison$cp_tracksig_adjusted)))
  print(mean(abs(cp_comparison$n_gt_created_cp - cp_comparison$cp_SCDS_results)))
  }

  #print("Examples where TrackSig makes mistakes but SciClone doesn't")
  #print(cp_comparison[(cp_comparison$n_gt_created_cp != cp_comparison$cp_tracksig) & (cp_comparison$n_gt_exposure_cp == cp_comparison$cp_sciclone),])

  #print(cp_comparison[(cp_comparison$n_gt_created_cp > cp_comparison$cp_tracksig) & (cp_comparison$n_gt_exposure_cp == cp_comparison$cp_sciclone),])


  print("Comparison of number of CP per sim type")
  sim_types <- c("one_cluster", "two_clusters",
      "branching", "cna_plus",
      "inf_site_viol_plus")

  depth_types <- c("depth10", "depth30", "depth100")
  res_table <- data.frame(matrix(0, ncol = length(sim_types), nrow = length(depth_types)))
  colnames(res_table) <- sim_types
  rownames(res_table) <- depth_types

  TrackSig_cp_summary <- SciClone_cp_summary <- res_table

  # Compute the same things for over-estimating and under-estimating number of subclones
  for (sim_type in sim_types) {
    for (d_type in depth_types) {
      idx <- grepl(sim_type, sapply(cp_comparison[,1],toString))
      idx <- idx &  grepl(paste0(d_type,"$"), sapply(cp_comparison[,1],toString))

      # print(sim_type)
      # print(d_type)
      # print(sum(idx))

      # print(sim_type)
      # print(d_type)
      # print(mean(cp_comparison[idx,]$cp_tracksig ))
      # print(mean(cp_comparison[idx,]$n_gt_created_cp ))
      # print("TrackSig agrees with sciclone")
      # print(mean(cp_comparison[idx,]$cp_tracksig == cp_comparison[idx,]$cp_sciclone))

      # print("TrackSig agrees with GT")
      # print(mean(cp_comparison[idx,]$n_gt_created_cp == cp_comparison[idx,]$cp_tracksig))
      # print("SciClone agrees with GT")
      # print(mean(cp_comparison[idx,]$n_gt_created_cp == cp_comparison[idx,]$cp_sciclone))

      TrackSig_cp_summary[d_type, sim_type] <- mean((cp_comparison[idx,]$n_gt_created_cp == cp_comparison[idx,]$cp_tracksig_adjusted) )
      SciClone_cp_summary[d_type, sim_type] <- mean((cp_comparison[idx,]$n_gt_created_cp == cp_comparison[idx,]$cp_SCDS_results) )
    }
  }

  print(TrackSig_cp_summary)
  print(SciClone_cp_summary)

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

  # ===========================================
  # Save the results for each simulation type and make a bar plot
  # For this, run sciclone with three different cluster methods: clusterMethod = "bmm", "gaussian.bmm", "binomial.bmm"
  # Finally, rename the folders to SCDS_results_binomial_bmm, SCDS_results_bmm or SCDS_results_gaussian_bmm

#  cp_comparison <- TrackSig:::compare_changepoints(simulations,
#    ground_truth_dir = dataDir,
#    tracksig_results_dir = paste0(tracksig_dir, "/SIMULATED/"),
#    sciclone_results_dir = paste0(c("SCDS_results_binomial", "SCDS_results_bmm") , "/SIMULATED/"),
#    res_file_name = "cp_comparison.txt",
#    change_at_cp_threshold = 0.05)
#
#  depth_types <- c("depth10", "depth30", "depth100")
#  sim_types <- c("one_cluster", "two_clusters",
#      "branching", "cna_plus",
#      "inf_site_viol_plus")
#  method_names <- c("cp_tracksig_adjusted", "cp_SCDS_results_binomial", "cp_SCDS_results_bmm")
#
#  empty_table <- data.frame(matrix(0, ncol = length(sim_types), nrow = length(method_names)))
#  colnames(empty_table) <- sim_types
#  rownames(empty_table) <- method_names
#
#  # cp_summary <- list()
#  # cp_summary[["cp_tracksig_adjusted"]] <- empty_table
#  # cp_summary[["cp_SCDS_results_binomial"]] <- empty_table
#  # cp_summary[["cp_SCDS_results_bmm"]] <- empty_table
#  # cp_summary[["cp_SCDS_results_gaussian_bmm"]] <- empty_table
#
#  res <- list()
#
#  # Compute the same things for over-estimating and under-estimating number of subclones
#  for (d_type in depth_types) {
#    res[[d_type]] <- empty_table
#
#    for (sim_type in sim_types) {
#      idx <- grepl(sim_type, sapply(cp_comparison[,1],toString))
#      idx <- idx &  grepl(paste0(d_type,"$"), sapply(cp_comparison[,1],toString))
#
#      for (name in method_names) {
#        percentage_correct <- mean((cp_comparison[idx,]$n_gt_created_cp == cp_comparison[idx,name]) )
#        res[[d_type]][name, sim_type] <- percentage_correct
#      }
#    }
#  }
#
#  print(res)
#
#  COLORS <- c("#FFFF00", "#FF4A46","#63FFAC","#B79762")
#   "#008941"
#
#  COLORS <- c("#F3766E", "#1FBFC3", "#C280F5")  #"#7CAA1F")
#
#  for (d_type in depth_types) {
#    pdf(paste0("cp_comparison_barplot_", d_type, ".pdf"), width = 7, height=5)
#    par(mar=c(6.1, 4.1, 6.1, 2.1))
#    barplot(as.matrix(res[[d_type]]), beside=T,
#                   col=COLORS[1:nrow(res[[d_type]])],
#                   names.arg= gsub("_", " ", colnames(res[[d_type]])) ) #las=2)
#    legend("topright", inset=c(0,-0.5), xpd=TRUE,  bty="n",
#         legend = c("TrackSig", "SciClone binomial.bmm", "SciClone bmm (default)"),
#         fill = COLORS[(1:nrow(res[[d_type]]))] )
#    dev.off()
#  }
#


  # ===========================================
  # Compare results for different bin sizes
  #dataDir <- "data_bin_simulations/"
  #
  #simulations <- list.files(dataDir)
  #sel <- grep(x = simulations, "^Simulation")
  #simulations <- simulations[sel]
  #
  #tracksig_dir = "TS_results_signature_trajectories/"
  #
  #bin_sizes <- c(25, 50, 75, 100, 150, 200, 300, 500)
  #
  #bin_size_results_kl <- c()
  #bin_size_results_mean_diff <- c()
  #bin_size_results_max_diff <- c()
  #
  #for (bin_size in bin_sizes){
  #  sim_bin_size_idx = as.integer(gsub(".*_bin(.*)_depth.*", "\\1", simulations)) == bin_size
  #  simulations_w_bin_size = simulations[sim_bin_size_idx]
  #
  #  list[res, gt_exp_l, estim_exp_l] <- compare_simulation_results(
  #      simulations_w_bin_size,
  #      ground_truth_dir = dataDir,
  #      method_results_dir = paste0(tracksig_dir, "/SIMULATED/"),
  #      res_file_name = sprintf("TrackSig_simulation_results_post%d.txt", bin_size))
  #
  #  bin_size_results_kl <- rbind(bin_size_results_kl, data.frame(bin_size=bin_size, metric=mean(res$kl)))
  #  bin_size_results_mean_diff <- rbind(bin_size_results_mean_diff, data.frame(bin_size=bin_size, metric=mean(res$abs_diff_mean)))
  #  bin_size_results_max_diff <- rbind(bin_size_results_max_diff, data.frame(bin_size=bin_size, metric=mean(res$abs_diff_max)))
  #
  #  # pdf(sprintf("TrackSig_KL_post%d.pdf", bin_size), width = 5, height=5)
  #
  #  # # ploting arguments
  #  # res$pch <- rep_len(c(20, 5), length.out = dim(res)[1])
  #  # res$col <- rep_len(c(2,2,3,3,4,4,5,5,6,6,7,7,8,8), length.out = dim(res)[1])
  #
  #  # plot(res$kl, res$abs_diff_max, main=sprintf("TrackSig KL post_bin_size = %d", bin_size),
  #  #    xlab="KL", ylab="max abs diff", pch = res$pch, col = res$col)
  #
  #  # pre_bin_sizes <- c(25, 50, 100, 125, 250)
  #
  #  # legend("bottomright", pch = c(20, 5, rep(15, 8)),
  #  #        col = c(1, 1, 2:8), legend = c("depth 100", "depth 1000", as.character(pre_bin_sizes)))
  #
  #  # res$labels <- strsplit(as.character(res$sim), "^Simulation_")
  #  # res$labels <- apply(res["labels"], MARGIN = 1, FUN=unlist)[2,]
  #  # res$labels <- strsplit(res$labels, "[[:digit:]]_[[:alnum:]]*$")
  #  # res$labels <- apply(res["labels"], MARGIN = 1, FUN=unlist)
  #
  #  # text(res$kl, res$abs_diff_max, labels = res$labels, cex = 0.8, pos = 4)
  #
  #  # dev.off()
  #}
  #
  #
  #print(bin_size_results_kl)
  #print(bin_size_results_mean_diff)
  #print(bin_size_results_max_diff)
  #
  #
  #pdf("TrackSig_bin_size_vs_KL.pdf", width = 5, height=5)
  #plot(bin_size_results_kl$bin_size, bin_size_results_kl$metric,
  #  main="", xlab="bin size", ylab="KL", pch=19)
  #dev.off()
  #
  #
  #pdf("TrackSig_bin_size_vs_mean_diff.pdf", width = 5, height=5)
  #plot(bin_size_results_mean_diff$bin_size, bin_size_results_mean_diff$metric,
  #  main="", xlab="bin size", ylab="mean activity diff", pch=19)
  #dev.off()


}


# [END]

