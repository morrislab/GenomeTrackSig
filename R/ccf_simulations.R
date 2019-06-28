#' Generate simulated data for \code{TrackSig}
#'
#' Functions required for generating simulated data for TrackSig
#'

#' \code{load_sim_signatures} <man content>
#' @rdname ccf_simulations
load_sim_signatures <- function(signature_file) {
	sigs_to_merge <- list()
	sigs_to_merge[["SBS7"]] <- c("SBS7a", "SBS7b", "SBS7c", "SBS7d")
	sigs_to_merge[["SBS17"]] <- c("SBS17a", "SBS17b")
	sigs_to_merge[["SBS2.13"]] <- c("SBS2", "SBS13")
	sigs_to_merge[["SBS10"]] <- c("SBS10a", "SBS10b")

  names_trinucleotide <- trinucleotide_internal
	names_trinucleotide <- apply(names_trinucleotide, 1, function(x) { do.call("paste", c(as.list(x), sep = "_"))})

	alex <- read.table(signature_file, stringsAsFactors = F, header=T)
    pcawg_trinucleotides <- paste(substr(alex[,1], 2, 2), substr(alex[,1], 5, 5), substr(alex[,1], 1, 3), sep="_")
    rownames(alex) <- pcawg_trinucleotides
    alex <- alex[,-1]
    alex <- alex[match(names_trinucleotide, pcawg_trinucleotides),]
	alex_merged <- t(merge_signatures(t(alex), sigs_to_merge))
  	return(alex_merged)
}

#' \code{compute_uncorrected_ccf} <man content>
#' @rdname ccf_simulations
compute_uncorrected_ccf <- function(n_alt_alleles, n_ref_alleles) {
	return(2 * n_alt_alleles / (n_alt_alleles + n_ref_alleles))
}

#' \code{save_exposures_per_time_point} <man content>
#' @rdname ccf_simulations
save_exposures_per_time_point <- function(data_all_clusters, sig_names, file_path, bin_size = 100) {
	# Make exposure matrix
	n_time_steps <- as.integer(nrow(data_all_clusters) / bin_size)
	mut_ccfs <- compute_uncorrected_ccf(data_all_clusters$n_alt_alleles, data_all_clusters$n_ref_alleles)
	data_all_clusters_sorted <- data_all_clusters[order(mut_ccfs, decreasing=T),]

	exposures <- c()
	ccfs <- c()
	for (i in 1:n_time_steps) {
		data_current_tp <- data_all_clusters_sorted[(bin_size*(i-1)+1):(bin_size * i), ]
		sigs_current_tp <- data_current_tp[, sig_names]

		mut_ccfs <- compute_uncorrected_ccf(data_current_tp$n_alt_alleles, data_current_tp$n_ref_alleles)
		time_point_ccf <- round(mean(mut_ccfs),3)

		exposures <- rbind(exposures, apply(sigs_current_tp, 2, mean))
		exposures <- round(exposures, 3)
		ccfs <- c(ccfs, time_point_ccf)
	}

	exposures <- t(exposures)
	colnames(exposures) <- ccfs
	write.csv(exposures, file = paste0(file_path, "_exposures.txt"), quote=F)
}

#' Generate simulations for mutations given clusters, their CCFs, CNAs and signatures.
##' @param n_clusters integer number of clusters of mutations to simulate
##' @param cluster_ccfs list of cluster mean values
##' @param n_mut_per_cluster list of number of mutations to generate per cluster
##' @param sig_activities list of signatures and their activities for each cluster
##' @param cluster_cna_info list of named list with elements "fractions", \cr
##' "mut_cn", "total_cn". There should be one inner list per cluster.
##' @param mean_depth average read depth to simulate> Defaults to 100
##' @param to_file logical whether to save generated simulations to a file. Defaluts to TRUE
##' @param simulation_name simulation name to appear in outputted file name. Defaults to "simulation"
##' @rdname ccf_simulations
##' @export

# Generate simulations for mutations given clusters, their CCFs, CNAs and signatures.
generate_ccf_simulation <- function(
	n_clusters, cluster_ccfs,
	n_mut_per_cluster, sig_activities,
	signature_def,
	cluster_cna_info = list(), mean_depth = 100, to_file = TRUE,
	simulation_name = "simulation", outdir = ".", bin_size){

	dir.create(outdir, showWarnings = FALSE)

	normal_mut_alleles = 1
	normal_total_CN = 2

	stopifnot(length(cluster_ccfs) == n_clusters)
	stopifnot(length(n_mut_per_cluster) == n_clusters)
	if (length(cluster_cna_info) > 0) {
		stopifnot(length(cluster_cna_info) == n_clusters)
	}

	trinucl_names <- rownames(signature_def)

	data_all_clusters <- list()
	for (cl in 1:n_clusters) {
		n_mut = n_mut_per_cluster[cl]
		cluster_ccf = cluster_ccfs[cl]

		cluster_sigs = sig_activities[[cl]]
		sig_names = names(cluster_sigs)
		stopifnot(sum(unlist(cluster_sigs)) - 1. < 0.001)

		# Sample depth from Poisson (lambda = mean depth)
		depths =  rpois(n_mut, lambda = mean_depth) # rep(mean_depth, n_mut)
		n_alt_alleles <- c()
		n_ref_alleles <- c()

		if (length(cluster_cna_info) == 0) {
			# No CNA changes
			prob = rep(cluster_ccf * (normal_mut_alleles / normal_total_CN), n_mut)
			mut_alleles = rep(normal_mut_alleles, n_mut)
			total_CN = rep(normal_total_CN, n_mut)
		} else {
			# CNA change
			current_cluster_cna = cluster_cna_info[[cl]]

			stopifnot(length(current_cluster_cna[["fractions"]]) == length(current_cluster_cna[["mut_cn"]]))
			stopifnot(length(current_cluster_cna[["fractions"]]) == length(current_cluster_cna[["total_cn"]]))
			stopifnot(sum(current_cluster_cna[["fractions"]]) == 1.0)

			prob = c()
			mut_alleles = c()
			total_CN = c()
			for (i in 1:length(current_cluster_cna[["fractions"]])) {
				current_frac = current_cluster_cna[["fractions"]][i]
				current_mut_cn = current_cluster_cna[["mut_cn"]][i]
				current_total_cn = current_cluster_cna[["total_cn"]][i]

				n_mut_with_cn_change = as.integer(n_mut * current_frac) + 1
				prob = c(prob, rep(cluster_ccf * (current_mut_cn / current_total_cn),   n_mut_with_cn_change))
				mut_alleles = c(mut_alleles, rep(current_mut_cn, n_mut_with_cn_change))
				total_CN = c(total_CN, rep(current_total_cn, n_mut_with_cn_change))
			}
			prob = prob[1:n_mut]
			mut_alleles = mut_alleles[1:n_mut]
			total_CN = total_CN[1:n_mut]
		}

		stopifnot(length(prob) == n_mut)

		# For each mutation in the cluster
		for (i in 1:n_mut) {
			# Sample number of variant alleles from a Binomial( d, CCF * (# mutant alleles / avg CN))
			alt <- rbinom(1, depths[i], prob[i])
			n_alt_alleles <- c(n_alt_alleles, alt)
			n_ref_alleles <- c(n_ref_alleles, depths[i] - alt)
		}

		# Sample mutation types for this cluster
		if (sum(sig_names %in% colnames(signature_def)) != length(sig_names)) {
			stop(paste("At least one of the signatures not found: ", toString(sig_names)))
		}

		mut_types <- c()
		# Sample mutations from each signature
		for (i in 1:length(cluster_sigs)) {
			sig_activity = cluster_sigs[[i]]
			mut_per_sig = as.integer(n_mut * sig_activity)+1

			sig_def <- signature_def[,names(cluster_sigs)[i]]

			# Get counts over mutation types
			mut_type_counts <- rmultinom(1, mut_per_sig, unlist(sig_def))
			rownames(mut_type_counts) <- trinucl_names

			# Create a list of mutation types out of counts
			mut_types_from_sig <- unlist(sapply(1:length(mut_type_counts), function(ind) {rep(trinucl_names[ind], mut_type_counts[ind]) }))
			stopifnot(length(mut_types_from_sig) == mut_per_sig)

			mut_types <- c(mut_types, mut_types_from_sig)
		}
		mut_types = mut_types[1:n_mut]

		# Shuffle list of mutation types
		mut_types <- sample(mut_types)
		stopifnot(length(mut_types) == n_mut)

		#  Report #variant, #reference, and avg CN for each mutation
		data <- data.frame(n_alt_alleles, n_ref_alleles, total_CN, mut_alleles, mut_types, data.frame(cluster_sigs), stringsAsFactors = F)
		if ("SBS2.13" %in% colnames(data)) {
      		colnames(data)[colnames(data) == "SBS2.13"] = "SBS2.13"
		}

		##################
		if(TrackSig.options()$simulation_pdf == T){

		  mut_ccfs <- compute_uncorrected_ccf(data$n_alt_alleles, data$n_ref_alleles)

		  pdf(paste0(outdir, "/", simulation_name, "_mut_ccf", cl, ".pdf"), width = 8, height=5)
		  hist(mut_ccfs, breaks=200)
		  dev.off()
		}
		##################

		data_all_clusters[[cl]] <- data
	}

	data_all_clusters_list <- data_all_clusters
	data_all_clusters_table <- do.call(rbind, data_all_clusters)
	##################
	if(TrackSig.options()$simulation_pdf == T){
	  mut_ccfs <- compute_uncorrected_ccf(data_all_clusters_table$n_alt_alleles, data_all_clusters_table$n_ref_alleles)

	  pdf(paste0(outdir, "/", simulation_name, "_mut_ccf.pdf"), width = 8, height=5)
	  hist(mut_ccfs, breaks=200)
	  dev.off()
	}
	##################

	if (to_file == TRUE){
		print(paste0("Saving simulation ",simulation_name))

		list[chrom, pos] = generate_chr_pos(n_mut = nrow(data_all_clusters_table))
		file_path = paste0(outdir, "/", simulation_name)
		save_as_vcf(data_all_clusters_table,  chrom, pos, file_path)

		sig_data <- data.frame(chrom, pos, data_all_clusters_table[,sig_names])
		colnames(sig_data) <- c("chromosome", "start", sig_names)
		write.table(sig_data, file = paste0(file_path, "_sig_exp_per_mut.txt"), sep = "\t", row.names=F, quote=F)

		save_exposures_per_time_point(data_all_clusters_table, sig_names, file_path, bin_size=bin_size)

		write_sim_summary(sim_data_all_clusters = data_all_clusters_list,
			simulation_name = simulation_name,
			sig_activities = do.call(rbind,sig_activities),
			ccfs=cluster_ccfs,
			mut_per_cluster = n_mut_per_cluster,
			file_path = file_path)
	}

	return(data_all_clusters_list)
}


generate_chr_pos <- function(n_mut) {

  # starts and ends specified with refernce hg19 - for simulation only!
  starts <- c(30028083, 149790583, 93504855, 75452280, 91686129, 60001, 282485,
              86726452, 92678797, 64340157, 1212760, 37856695, 19020001, 19000001,
              29209444, 46385802, 34725849, 18510899, 27731783, 60001, 14338130, 20609432)

  ends <- c(103863906, 234003741, 194041961, 191044276, 138787073, 58087659,
            50370631, 142766515, 133073060, 116065824, 50783853, 109373470,
            86760324, 107289540, 82829645, 88389383, 62410760, 52059136,
            59118983, 26319569, 33157035, 50364777)


  chrN <- sample(1:22, n_mut, replace = TRUE)
	chrom <- paste0("chr", chrN)
	pos <- round(runif(n_mut, starts[chrN], ends[chrN]))

	return(list(chrom, pos))
}


# save data generated by generate_ccf_simulation() to a file
save_as_vcf <- function(data, chrom, pos, filename) {
	splitted_types <- t(data.frame(strsplit(data$mut_types, "_")))
	rownames(splitted_types) <- NULL
	ref <- splitted_types[,1]
	alt <- splitted_types[,2]
	tri <- splitted_types[,3]

	info <- paste0("t_ref_count=", data$n_ref_alleles, ";t_alt_count=", data$n_alt_alleles)
	vcf_data <- data.frame(chrom, pos, ".", ref, alt, ".", ".",  info)
	colnames(vcf_data) <- c("#CHROM", "POS", "",  "REF", "ALT", "", "FILTER",  "INFO")

	cna_data <- data.frame(chrom, pos, pos+1, data$total_CN)
	colnames(cna_data) <- c("chromosome", "start", "end", "total_cn")

	tri_data <- data.frame(chrom, pos, ref, alt, tri)
	colnames(tri_data) <- c("chromosome", "start", "ref", "alt", "tri")

	scicloneDS_data <- data.frame(chrom, pos, data$n_ref_alleles, data$n_alt_alleles,
	                            data$n_alt_alleles / (data$n_ref_alleles + data$n_alt_alleles), ref, alt, tri)
	colnames(scicloneDS_data) <- c("chr", "pos", "ref_reads", "var_reads", "vaf", "ref", "alt", "tri")

	write.table(vcf_data, file = paste0(filename, ".vcf"), sep = "\t", row.names=F, quote=F)
	write.table(cna_data, file = paste0(filename, "_cna.txt"), sep = "\t", row.names=F, quote=F)
	write.table(tri_data, file = paste0(filename, "_tri.txt"), sep = "\t", row.names=F, quote=F)
	write.table(scicloneDS_data, file = paste0(filename, "_scicloneDS.txt"), sep = "\t", row.names=F, quote=F)
}


#' \code{write_sim_annotation} save active signatures generated by generate_ccf_simulation() to a file
#' Uses PCAWG_sigProfiler signature set
##' @rdname simulation_functions
write_sim_annotation <- function(simulation_name, sig_activities, sig_header,
                                 sim_activity_file, sim_purity_file, sim_tumortype_file) {
  # binarize activity
  activity <- c("SIMULATED", simulation_name, rep(0, (length(sig_header) - 2)))
  activity[sig_header %in% unique(names(unlist(sig_activities)))] <- 1

  if (sum(as.integer(activity[3:length(activity)])) != ncol(do.call(rbind,sig_activities))) {
  	print("Cannot write signature activities: ")
  	print(sig_activities)
  	stop()
  }

  # write activity
  write.table(t(activity), file = sim_activity_file, append = T,
              col.names = F, row.names = F, quote = F, sep = "\t")

  # write purity
  write.table(t(c(simulation_name, 1)), file = sim_purity_file, append = T,
              col.names = F, row.names = F, quote = F, sep = "\t")

  # write tumor type
  write.table(t(c(simulation_name, "SIMULATED")), file = sim_tumortype_file, append = T,
              col.names = F, row.names = F, quote = F, sep = "\t")
}


sample_sigs_and_activities <- function(meaningful_sig1, meaningful_sig2, sig1_range=c(0.2, 0.7)) {
	exposure_time_sig = runif(1, min=0.03, max=0.1)
	cluster_sigs = list("SBS1" = exposure_time_sig)
	cluster_sigs["SBS5"] =  runif(1, min= 0.05, max= 0.15)
	remaining_activity = 1 - exposure_time_sig - cluster_sigs[["SBS5"]]

	meaningful_sig_activity1 = runif(1, min=sig1_range[1],
										max=min(remaining_activity, sig1_range[2]))

	# meaningful_sig_activity2 = runif(1,
	# 	min= min(0.1, 1- exposure_time_sig - meaningful_sig_activity1),
	# 	max= 1- exposure_time_sig - meaningful_sig_activity1)

	# cluster_sigs["SBS5"] = 1 - exposure_time_sig - meaningful_sig_activity1 - meaningful_sig_activity2

	meaningful_sig_activity2 = 1 - exposure_time_sig - meaningful_sig_activity1 - cluster_sigs[["SBS5"]]

	cluster_sigs[meaningful_sig1] = meaningful_sig_activity1
	cluster_sigs[meaningful_sig2] = meaningful_sig_activity2

	stopifnot(cluster_sigs[["SBS5"]] >= 0)
	stopifnot(meaningful_sig_activity2 > 0)
	stopifnot(sum(unlist(cluster_sigs)) - 1 < 0.001)

	return(cluster_sigs)
}

write_sim_summary <- function(sim_data_all_clusters, simulation_name, sig_activities,
			ccfs, mut_per_cluster, file_path) {

	file_name = paste0(file_path, "_info.txt")
	unlink(file_name)

	write(simulation_name, file_name, append=T)

	write("Sig activities", file_name, append=T)
	suppressWarnings(write.table(sig_activities, file_name, append=T, quote=F))

	write("CCFs per cluster", file_name, append=T)
	write(toString(ccfs), file_name, append=T)

	write("Mutation counts per cluster", file_name, append=T)
	write(toString(mut_per_cluster), file_name, append=T)
}


##' \code{create_simulation_set} generates simulated data files
##' @rdname create_simulation_set
##' @param mut_per_sim int. Number of mutations per simulation. Default: 5000
##' @export
create_simulation_set <- function(outdir = "simulations", mut_per_sim = 5000,
                                  sim_activity_file = "annotation/sim_active_in_sample.txt",
                                  sim_purity_file = "annotation/sim_purity.txt",
                                  sim_tumortype_file = "annotation/sim_tumortypes.txt",
                                  signature_file = "annotation/sigProfiler_SBS_signatures.txt",
                                  rewrite_annotations=T) {

	dir.create(outdir, showWarnings = FALSE)
	set.seed(2019)

	signature_def = load_sim_signatures(signature_file)

	if (rewrite_annotations) {
	  # remove simulation annotations (rebuilt upon simulation) and create new ones
	  unlink(sim_activity_file)
	  unlink(sim_purity_file)
	  unlink(sim_tumortype_file)

	  # write headers for sim annotation files
	  write.table(t(c("samplename", "purity")), file = sim_purity_file,
		            col.names = F, row.names = F, quote = F, sep = "\t")

	  write.table(t(c("ID", "tumortype")), file = sim_tumortype_file,
		            col.names = F, row.names = F, quote = F, sep = "\t")

	  sig_header <- colnames(signature_def)
	  sig_header <- c("Cancer_Type", "Sample_Name", sig_header)

	  write.table(t(sig_header), file = sim_activity_file,
		            col.names = F, row.names = F, quote = F, sep = "\t")
	}

  sim_list = c()

  # Variable signatures to sample from
  # meaningful_sig_list = c("SBS2.13", "SBS3", "SBS4", "SBS6", "SBS7", "SBS9")
  # Signature SBS7 is excluded from the list because both TrackSig and SciClone don't perform well with it.
  meaningful_sig_list<-  c("SBS3", "SBS4", "SBS6", "SBS8", "SBS9",
		"SBS11", "SBS12", "SBS14", "SBS15", "SBS16", "SBS18", "SBS19",
		"SBS20", "SBS21", "SBS22", "SBS23", "SBS24", "SBS25", "SBS26",
		"SBS27", "SBS28", "SBS29", "SBS30", "SBS31", "SBS32", "SBS33",
		"SBS34", "SBS35", "SBS36", "SBS37", "SBS38", "SBS39", "SBS40",
		"SBS17", "SBS2.13", "SBS10")


  print("Simulation type 0b: two clusters")
	# signature does not change, but CCFs do

	n_simulations = 225    # MUST be a square number

	dists <- seq(0.85, 0.2, length.out = sqrt(n_simulations))       # clusters of increasing distance apart
	sigAdds <- seq(0.05, 0.5, length.out = sqrt(n_simulations))     # signatures of increaing change
	bin_sizes = c(100)
	depth_list <- c(100)

	# combo indices to use
	dist_i <- rep(1:sqrt(n_simulations), each = sqrt(n_simulations))
	sigAdd_i <- rep(1:sqrt(n_simulations), times = sqrt(n_simulations))


	for (bin_size in bin_sizes){

	  for (sim_i in 1:n_simulations) {

	  	sig_activities = list()

	  	# indexing sigAdds and dists
	  	dist <- dists[dist_i[sim_i]]
	  	sigAdd <- sigAdds[sigAdd_i[sim_i]]

	  	# Sample signatures with variable presence
	  	list[meaningful_sig1, meaningful_sig2] = sample(meaningful_sig_list, size = 2)

	  	# Signatures change in cluster 2, but not in cluster 1
	  	clonal_sigs <- sample_sigs_and_activities(meaningful_sig1, meaningful_sig2, sig1_range=c(0.6, 0.6))

	  	sig_activities[[1]] <- clonal_sigs
	  	#sig_activities[[2]] <- sample_sigs_and_activities(meaningful_sig1, meaningful_sig2, sig1_range=c(0.5 + sigAdd, 0.5 + sigAdd))
	  	sig_activities[[2]] <- clonal_sigs
	  	sig_activities[[2]][[3]] <- sig_activities[[2]][[3]] - sigAdd
	  	sig_activities[[2]][[4]] <- sig_activities[[2]][[4]] + sigAdd


	  	print(do.call(rbind,sig_activities))

	  	#subclone1_ccf = runif(1, min=0.2, max=0.6)
	  	subclone1_ccf <- dist

	  	print("CCFs per cluster")
	  	print(c(1.0, subclone1_ccf))

	  	# cluster 1 and cluster 2 and two separate branches
	  	n_mut_subclone1 = as.integer(mut_per_sim * subclone1_ccf)
	  	n_mut_clonal = mut_per_sim - n_mut_subclone1

	  	print("Mutation counts per cluster")
	  	print(c(n_mut_clonal, n_mut_subclone1))

	  	for (depth in depth_list) {

	  		simulation_name = paste0("Simulation_two_clusters",
	  			sim_i, "_depth", depth, "_bin", bin_size,
	  			"_dist", dist, "_sigChange", sigAdd)

	  		print(paste0("Generating simulation ",simulation_name))

	  		sim_data_all_clusters = generate_ccf_simulation(
	  			n_clusters = 2,
	  			cluster_ccfs = c(1.0, subclone1_ccf),
	  			n_mut_per_cluster = c(n_mut_clonal, n_mut_subclone1),
	  			sig_activities = sig_activities,
	  			signature_def = signature_def,
	  			simulation_name = simulation_name,
	  			mean_depth = depth,
	  			outdir = paste0(outdir, "/", simulation_name),
	  			bin_size = bin_size)

	  		write_sim_annotation(simulation_name, sig_activities, sig_header,
	  			sim_activity_file, sim_purity_file, sim_tumortype_file)

	  		sim_list <- c(sim_list, simulation_name)
	  	}
	  }
	}

	print("Created simulations:")
	print(sim_list)
	return(sim_list)
}


