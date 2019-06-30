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

kldiv_multinomials <- function(multinom1, multinom2) {
  return(apply(multinom1 * log(multinom1/multinom2),1,sum))
}

compare_simulation_results  <- function(simulation_list,
  ground_truth_dir, method_results_dir, res_file_name) {

  results_df <- c()
  gt_exposures_list <- list()
  estim_exposures_list <- list()

  for (sim in simulation_list) {
    print(sprintf("Processing simulation %s ...", sim))

    gt_exposures_file = paste0(ground_truth_dir, "/", sim, "/", sim, "_sig_exp_per_mut.txt")

    if (!file.exists(gt_exposures_file)) {
      print(sprintf("File %s not found", gt_exposures_file))
      next
    }

    gt_exposures <- read.delim(gt_exposures_file, header=T, stringsAsFactors=F)
    gt_pos = paste0(gt_exposures[,"chromosome"], "_", gt_exposures[,"start"])
    rownames(gt_exposures) <- gt_pos

    estim_exposures_file = paste0(method_results_dir, "/", sim, "/", "sig_exposures_per_mut.txt")

    if (!file.exists(estim_exposures_file)) {
      print(sprintf("File %s not found", estim_exposures_file))
      next
    }

    estim_exposures <- read.delim(estim_exposures_file, header=T, stringsAsFactors=F)

    # Making the rows be in the same order
    estim_pos = paste0(estim_exposures[,"chromosome"], "_", estim_exposures[,"start"])
    rownames(estim_exposures) <- estim_pos
    gt_pos <- intersect(gt_pos, estim_pos)
    estim_exposures = estim_exposures[gt_pos,]
    gt_exposures = gt_exposures[gt_pos,]

    # Making the signatures be in the same order
    stopifnot(sum(sort(colnames(estim_exposures)) == sort(colnames(gt_exposures))) == ncol(estim_exposures))
    estim_exposures = estim_exposures[,colnames(gt_exposures)]
    stopifnot(colnames(estim_exposures) == colnames(gt_exposures))

    gt_exposures_d = gt_exposures[,3:ncol(gt_exposures)]
    estim_exposures_d = estim_exposures[,3:ncol(estim_exposures)]

    abs_diff <- abs(gt_exposures_d - estim_exposures_d)

    estim_exposures_d_eps <- estim_exposures_d
    estim_exposures_d_eps[estim_exposures_d_eps == 0] = 0.01
    estim_exposures_d_eps = estim_exposures_d_eps / apply(estim_exposures_d_eps, 1, sum)

    sim_type <- gsub("Simulation_([A-z_]*)\\d\\d*_.*", "\\1", sim)
    depth <- as.integer(gsub(".*_depth(.*)", "\\1", sim))

    d <- data.frame(sim = sim,
      abs_diff_mean = mean(unlist(abs_diff), na.rm=TRUE),
      abs_diff_max = max(unlist(abs_diff), na.rm=TRUE),
      kl = mean(kldiv_multinomials(gt_exposures_d, estim_exposures_d_eps), na.rm=TRUE),
      sim_type = sim_type,
      depth = depth,
      stringsAsFactors = F)

    gt_exposures_list[[sim]] <- gt_exposures
    estim_exposures_list[[sim]] <- estim_exposures

    results_df <- rbind(results_df, d)
  }

  results_df <- data.frame(results_df)
  rownames(results_df) <- results_df[,1]
  write.table(results_df, file = res_file_name, sep = "\t", row.names=F, quote=F)

  return(list(results_df, gt_exposures_list, estim_exposures_list))
}


plot_kl_results <- function(res, method, sim_type = "") {
  if (sim_type != "") {
    sim_type = paste0(".", sim_type)
  }
  pdf(paste0(method, "_simulation_results_KL_vs_max", sim_type, ".pdf"), width = 5, height=5)
  plot(res$kl, res$abs_diff_max, main=paste0(method, " versus Truth"),
     xlab="KL", ylab="max abs diff", xlim=c(0, 0.5), ylim=c(0, 0.5))
  dev.off()

  pdf(paste0(method,"_simulation_results_KL_vs_mean", sim_type, ".pdf"), width = 5, height=5)
  plot(res$kl, res$abs_diff_mean, main=paste0(method, " versus Truth"),
     xlab="KL", ylab="mean abs diff", xlim=c(0, 0.5), ylim=c(0, 0.5))
  dev.off()

}

get_sig_names_from_list <- function(exp_per_mut_list) {
  return(lapply(exp_per_mut_list, function(x) colnames(x[3:6])))
}


get_change_points <- function(mat) {
  return(c(FALSE, apply(abs(mat[,2:ncol(mat)] - mat[,1:(ncol(mat)-1)]),2,sum) != 0))
}

get_max_change_at_cp <- function(mixtures, change_points) {
    changes_at_cp <- abs(mixtures[,unlist(change_points)] - mixtures[,unlist(change_points-1)])
    if (is.vector(changes_at_cp)) {
      max_change_at_cp = max(changes_at_cp)
    } else {
      max_change_at_cp <- apply(changes_at_cp,2,max)
    }
    return(max_change_at_cp)
 }


adjust_change_points_ <- function(change_points, mixtures, change_at_cp_threshold = 0.05) {
    # filter change points: remove those where change < 0.05
    if (length(change_points) > 0) {
     max_change_at_cp <- get_max_change_at_cp(mixtures, change_points)

      # remove change points where the mixture change is less than 5% exposure -- those negligible changes
      if (sum(max_change_at_cp > change_at_cp_threshold) > 0) {
        change_points <- change_points[max_change_at_cp > change_at_cp_threshold]
      } else {
        change_points <- list()
      }
    }
    return(change_points)
}




get_max_changes_one_tumor <- function(mixtures, signatures_to_compute, tumor_id,
  change_at_cp_threshold = 0.06) {
  new_item <- toHorizontalMatrix(rep(0, length(signatures_to_compute)))

  colnames(new_item) <- signatures_to_compute

  rownames(new_item) <- tumor_id

  new_item_mean <- new_item_direction <- new_item

  new_item[,rownames(mixtures)] <- apply(mixtures, 1, max) - apply(mixtures, 1, min)

  signature_big_change = names(which(apply(mixtures, 1, max) - apply(mixtures, 1, min) >= change_at_cp_threshold)) # 0.06
  signature_neg_change = names(which(apply(mixtures, 1, which.max) < apply(mixtures, 1, which.min)))
  signature_pos_change = names(which(apply(mixtures, 1, which.max) > apply(mixtures, 1, which.min)))

  signature_neg_change <- intersect(signature_neg_change, signature_big_change)
  signature_pos_change <- intersect(signature_pos_change, signature_big_change)

  new_item_direction[,signature_pos_change] <- 1
  new_item_direction[,signature_neg_change] <- -1

  new_item_mean[,rownames(mixtures)] <- apply(mixtures, 1, mean)

  return(list(new_item, new_item_direction, new_item_mean))
}


toHorizontalMatrix <- function(L){
  if (is.vector(L))
    return(matrix(L, nrow=1))
  else
    return(as.matrix(L))
}

toVerticalMatrix <- function(L)
{
  if (is.vector(L))
    return(matrix(L, ncol=1))
  else
    return(as.matrix(L))
}

remove_cp_with_same_direction <- function(mixtures, change_points) {
  sigs <- rownames(mixtures)
  meaningful_sigs <- setdiff(sigs, c("SBS1", "SBS5"))
  mixtures <- mixtures[meaningful_sigs, ]

  if (length(change_points) > 1) {
    direction <- mixtures[,unlist(change_points)] - mixtures[,unlist(change_points-1)] > 0

    # remove change-points that have the same direction of sig change as the previous change-points
    direction_diff <- abs(toVerticalMatrix(direction[,2:ncol(direction)]) - toVerticalMatrix(direction[,1:(ncol(direction)-1)]))
    to_remove <- which(apply(direction_diff, 2, sum) == 0)
    if (length(to_remove) > 0) {
      change_points <- change_points[-to_remove]
    }
  }
  return(change_points)
}

compare_changepoints  <- function(simulation_list, ground_truth_dir,
  tracksig_results_dir, sciclone_results_dir_list, res_file_name,
  change_at_cp_threshold = 0.05) {

  if (is.character(sciclone_results_dir_list)) {
    # input for sciclone dir just one string
    sciclone_results_dir_list = c(sciclone_results_dir_list)
  } else {
    # dirs for sciclone are specified as a list of strings
    # do nothing?
  }


  results_df <- c()

  for (sim in simulation_list) {
    print(sprintf("Processing simulation %s ...", sim))
    tracksig_cp_file = paste0(tracksig_results_dir, "/", sim, "/changepoints.txt")

    if (!file.exists(tracksig_cp_file)) {
      print(sprintf("File %s not found", tracksig_cp_file))
      next
    }

   if (file.info(tracksig_cp_file)$size == 1) {
      cp_tracksig <- 0
    } else {
      cp_tracksig <- ncol(read.table(tracksig_cp_file, header=F))
    }

    mixture_file <- paste0(tracksig_results_dir, "/", sim, "/mixtures.csv")
    mixtures_tracksig <- read.csv(mixture_file, header=T, stringsAsFactors=F)
    rownames(mixtures_tracksig) <- mixtures_tracksig[,1]
    mixtures_tracksig <- mixtures_tracksig[,-1]

    cp_tracksig_unadjusted = get_change_points(mixtures_tracksig)

    adjusted_cp_tracksig = adjust_change_points_(which(cp_tracksig_unadjusted), mixtures_tracksig,
      change_at_cp_threshold = change_at_cp_threshold)

    adjusted_cp_tracksig <- remove_cp_with_same_direction(mixtures_tracksig, adjusted_cp_tracksig)

    # read change-points from sciclone
    # If we ran sciclone using multiple clustering methods, we want to aggregate results from all of them
    cp_sciclone <- list()
    for (sciclone_results_dir in sciclone_results_dir_list) {
      sciclone_cp_file = paste0(sciclone_results_dir, "/", sim, "/phis.txt")

      if (!file.exists(sciclone_cp_file)) {
        print(sprintf("File %s not found", sciclone_cp_file))
        next
      }
      # counting change-points, not the number of clusters
      method_name = paste0("cp_", gsub("([A-z_]*)//*SIMULATED/", "\\1", sciclone_results_dir))
      cp_sciclone[[method_name]] <- nrow(read.table(sciclone_cp_file, header=F, stringsAsFactors=F)) - 1
    }

    gt_exposures_file = paste0(ground_truth_dir, "/", sim, "/", sim, "_exposures.txt")

    if (!file.exists(gt_exposures_file)) {
      print(sprintf("File %s not found", gt_exposures_file))
      next
    }

    gt_exposures <- read.csv(gt_exposures_file, header=T, stringsAsFactors=F)
    rownames(gt_exposures) <- gt_exposures[,1]
    gt_exposures <- gt_exposures[,-1]

    gt_delta_exp <- abs(gt_exposures[,2:ncol(gt_exposures)]-gt_exposures[,1:(ncol(gt_exposures)-1)])
    n_gt_exposure_cp <- sum(apply(gt_delta_exp,2,max) > 0.05)

    if (grepl("two_clusters", sim)) {
      n_gt_created_cp = 1
    } else if (grepl("one_cluster", sim)) {
      n_gt_created_cp = 0
    } else if (grepl("inf_site_viol", sim)) {
      n_gt_created_cp = 2
    } else {
      n_gt_created_cp = 2
    }

    d <- data.frame(sim = sim,
      cp_tracksig = cp_tracksig,
      cp_tracksig_adjusted = length(adjusted_cp_tracksig),
      n_gt_exposure_cp = n_gt_exposure_cp,
      n_gt_created_cp = n_gt_created_cp)

    for (i in 1:length(cp_sciclone)) {
      d[[names(cp_sciclone)[i]]] <- cp_sciclone[[i]]
    }

    results_df <- rbind(results_df, d)
  }

  rownames(results_df) <- results_df[,1]
  write.table(results_df, file = res_file_name, sep = "\t", row.names=F, quote=F)

  return(results_df)
}


get_results_per_sim_type <- function(res, method_name) {
  sim_types <- c("one_cluster", "two_clusters",
    "branching", "cna_plus",
    "inf_site_viol_plus")

  print("Per sim type")
  for (type in sim_types) {
    idx <- grepl(type, sapply(res[,1],toString))

    if (sum(idx) > 1) {
      plot_kl_results(res[idx,], method_name, sim_type = type)
    }
    print(paste0("Percentage of samples with KL larger than 0.05 in type ", type))
    print(mean(res[idx,]$kl > 0.05))
  }

  print("Per depth")
  depth_types <- c("depth10", "depth30", "depth100")
  for (d_type in depth_types) {
    idx <- grepl(paste0(d_type,"$"), sapply(res[,1],toString))

    if (sum(idx) > 1) {
      plot_kl_results(res[idx,], method_name, sim_type = d_type)
    }
    print(paste0("Percentage of samples with KL larger than 0.05 in type ", d_type))
    print(mean(res[idx,]$kl > 0.05))
  }

  avg_kl_per_type <- data.frame(matrix(0, ncol = length(sim_types), nrow = length(depth_types)))
  colnames(avg_kl_per_type) <- sim_types
  rownames(avg_kl_per_type) <- depth_types

  for (sim_type in sim_types) {
    for (d_type in depth_types) {
      idx <- grepl(sim_type, sapply(res[,1],toString))
      idx <- idx &  grepl(paste0(d_type,"$"), sapply(res[,1],toString))

      avg_kl_per_type[d_type, sim_type] <- mean(res[idx,]$kl)
    }
  }
  return(avg_kl_per_type)
}
