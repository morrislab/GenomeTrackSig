kldiv_multinomials <- function(multinom1, multinom2) {
  return(apply(multinom1 * log(multinom1/multinom2),1,sum))
}

compare_simulation_results  <- function(simulation_list, 
  ground_truth_dir, method_results_dir, res_file_name) {

  results_df <- c()

  for (sim in simulation_list) {
    print(sprintf("Processing simulation %s ...", sim))

    gt_exposures_file = paste0(ground_truth_dir, "/", sim, "/", sim, "_sig_exp_per_mut.txt")

    if (!file.exists(gt_exposures_file)) {
      print(sprintf("File %s not found", gt_exposures_file))
      next
    }

    gt_exposures <- read.delim(gt_exposures_file, header=T, stringsAsFactors=F)
    gt_pos = paste0(gt_exposures[,"chromosome"], "_", gt_exposures[,"start"])

    estim_exposures_file = paste0(method_results_dir, "/", sim, "/", "sig_exposures_per_mut.txt")

    if (!file.exists(estim_exposures_file)) {
      print(sprintf("File %s not found", estim_exposures_file))
      next
    }

    estim_exposures <- read.delim(estim_exposures_file, header=T, stringsAsFactors=F)
    stopifnot(dim(estim_exposures) == dim(gt_exposures))

    # Making the rows be in the same order
    estim_pos = paste0(estim_exposures[,"chromosome"], "_", estim_exposures[,"start"])
    rownames(estim_exposures) <- estim_pos
    estim_exposures = estim_exposures[gt_pos,]

    # Making the signatures be in the same order
    estim_exposures = estim_exposures[,colnames(gt_exposures)]
    stopifnot(colnames(estim_exposures) == colnames(gt_exposures))

    gt_exposures_d = gt_exposures[,3:ncol(gt_exposures)]
    estim_exposures_d = estim_exposures[,3:ncol(estim_exposures)]

    abs_diff <- abs(gt_exposures_d - estim_exposures_d)

    d <- data.frame(sim = sim, 
      abs_diff_mean = mean(unlist(abs_diff)),
      abs_diff_max = max(unlist(abs_diff)),
      kl = mean(kldiv_multinomials(gt_exposures_d, estim_exposures_d)))

    results_df <- rbind(results_df, d)
  }

  results_df <- data.frame(results_df)
  write.table(results_df, file = res_file_name, sep = "\t", row.names=F, quote=F)
  
  return(results_df)
}