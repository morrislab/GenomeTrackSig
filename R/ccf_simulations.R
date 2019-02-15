#' Generate simulated data for \code{TrackSig}
#'
#' Functions required for generating simulated data for TrackSig
#'


toHorizontalMatrix <- function(L)
{
  if (is.vector(L))
    return(matrix(L, nrow=1))
  else 
    return(as.matrix(L))
}


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


# When I will upload this on the server -- remove merge_signatures function -- use src/helper_functions.R
merge_signatures <- function(mixtures, sigs_to_merge) {
  if (!is.null(sigs_to_merge)) {
    for (i in 1:length(sigs_to_merge)) {
        set_name <- names(sigs_to_merge)[i]
        sig_set = sigs_to_merge[[i]]
        sig_set = intersect(rownames(mixtures), sig_set)

        if (length(sig_set) == 0) {
          next
        }

        new_mixture_col <- toHorizontalMatrix(apply(mixtures[sig_set,,drop=F],2,sum))
        rownames(new_mixture_col) <- set_name
        colnames(new_mixture_col) <- colnames(mixtures)

        if (length(intersect(rownames(mixtures), sig_set)) == nrow(mixtures)) {
          mixtures <- new_mixture_col
        } else {
          mixtures <- rbind(mixtures[-which(rownames(mixtures) %in% sig_set),,drop=FALSE], new_mixture_col)
        }
    }
  }
  return(mixtures)
}


load_sim_signatures <- function(signature_file, trinucleotide_file) {
	sigs_to_merge <- list()
	sigs_to_merge[["SBS7"]] <- c("SBS7a", "SBS7b", "SBS7c", "SBS7d")
	sigs_to_merge[["SBS17"]] <- c("SBS17a", "SBS17b")
	sigs_to_merge[["SBS2+13"]] <- c("SBS2", "SBS13")
	sigs_to_merge[["SBS10"]] <- c("SBS10a", "SBS10b")

  	names_trinucleotide <- read.table(trinucleotide_file, stringsAsFactors = F)
	names_trinucleotide <- apply(names_trinucleotide, 1, function(x) { do.call("paste", c(as.list(x), sep = "_"))})
  
	alex <- read.table(signature_file, stringsAsFactors = F, header=T)
    pcawg_trinucleotides <- paste(substr(alex[,1], 2, 2), substr(alex[,1], 5, 5), substr(alex[,1], 1, 3), sep="_")
    rownames(alex) <- pcawg_trinucleotides
    alex <- alex[,-1]
    alex <- alex[match(names_trinucleotide, pcawg_trinucleotides),]
	alex_merged <- t(merge_signatures(t(alex), sigs_to_merge))
  	return(alex_merged)
}

compute_uncorrected_ccf <- function(n_alt_alleles, n_ref_alleles) {
	return(2 * n_alt_alleles / (n_alt_alleles + n_ref_alleles))
}


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
##' @rdname simulation_functions
 
# Generate simulations for mutations given clusters, their CCFs, CNAs and signatures.
generate_ccf_simulation <- function(
	n_clusters, cluster_ccfs, 
	n_mut_per_cluster, sig_activities,
	signature_def,
	cluster_cna_info = list(), mean_depth = 100, to_file = TRUE,
	simulation_name = "simulation", outdir = "."){

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
      		colnames(data)[colnames(data) == "SBS2.13"] = "SBS2+13" 
		}

		##################
		mut_ccfs <- compute_uncorrected_ccf(data$n_alt_alleles, data$n_ref_alleles)

		pdf(paste0(outdir, "/", simulation_name, "_mut_ccf", cl, ".pdf"), width = 8, height=5)
		hist(mut_ccfs, breaks=200)
		dev.off()
		##################

		data_all_clusters[[cl]] <- data
	}

	data_all_clusters_list <- data_all_clusters
	data_all_clusters_table <- do.call(rbind, data_all_clusters)
	##################
	mut_ccfs <- compute_uncorrected_ccf(data_all_clusters_table$n_alt_alleles, data_all_clusters_table$n_ref_alleles)

	pdf(paste0(outdir, "/", simulation_name, "_mut_ccf.pdf"), width = 8, height=5)
	hist(mut_ccfs, breaks=200)
	dev.off()
	##################

	if (to_file == TRUE){
		print(paste0("Saving simulation ",simulation_name))

		list[chrom, pos] = generate_chr_pos(n_mut = nrow(data_all_clusters_table))
		file_path = paste0(outdir, "/", simulation_name)
		save_as_vcf(data_all_clusters_table,  chrom, pos, file_path)
		
		sig_data <- data.frame(chrom, pos, data_all_clusters_table[,sig_names])
		colnames(sig_data) <- c("chromosome", "start", sig_names)
		write.table(sig_data, file = paste0(file_path, "_sig_exp_per_mut.txt"), sep = "\t", row.names=F, quote=F)

		save_exposures_per_time_point(data_all_clusters_table, sig_names, file_path, bin_size=100)
			
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
	chrom <- paste0("chr", sample(1:22, n_mut, replace = TRUE))
	pos <- sample(10**5, n_mut)
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

	write.table(vcf_data, file = paste0(filename, ".vcf"), sep = "\t", row.names=F, quote=F)
	write.table(cna_data, file = paste0(filename, "_cna.txt"), sep = "\t", row.names=F, quote=F)
	write.table(tri_data, file = paste0(filename, "_tri.txt"), sep = "\t", row.names=F, quote=F)
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

	meaningful_sig_activity1 = runif(1, min=sig1_range[1], max=sig1_range[2])
	
	# meaningful_sig_activity2 = runif(1, 
	# 	min= min(0.1, 1- exposure_time_sig - meaningful_sig_activity1), 
	# 	max= 1- exposure_time_sig - meaningful_sig_activity1)

	# cluster_sigs["SBS5"] = 1 - exposure_time_sig - meaningful_sig_activity1 - meaningful_sig_activity2

	cluster_sigs["SBS5"] =  runif(1, min= 0.05, max= 0.15)

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
                                  trinucleotide_file = "annotation/trinucleotide.txt",
                                  rewrite_annotations=T) {

	dir.create(outdir, showWarnings = FALSE)
	set.seed(2019)

	signature_def = load_sim_signatures(signature_file, trinucleotide_file)

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
    print("Simulation type 0a: one cluster")
	# signatures change in one cluster but not in the other"
	n_simulations = 5
	for (sim_id in 1:n_simulations) {
		sig_activities = list()

		meaningful_sig_list = c("SBS2+13", "SBS3", "SBS4", "SBS6", "SBS7", "SBS9")
		list[meaningful_sig1, meaningful_sig2] = sample(meaningful_sig_list, size = 2)
		
		# Signatures change in cluster 2, but not in cluster 1
		sig_activities[[1]] <- sample_sigs_and_activities(meaningful_sig1, meaningful_sig2, sig1_range=c(0.4, 0.8))

		print("Sig activities")
		print(do.call(rbind,sig_activities))
		
		print("Mutation counts per cluster")
		print(c(mut_per_sim))

		for (depth in c(100, 1000)) {
			simulation_name = paste0("Simulation_one_cluster", sim_id, "_depth", depth)
			print(paste0("Generating simulation ",simulation_name))

			sim_data_all_clusters = generate_ccf_simulation(
				n_clusters = 1, 
				cluster_ccfs = c(1.0), 
				n_mut_per_cluster = c(mut_per_sim), 
				sig_activities = sig_activities,
				signature_def = signature_def,
				simulation_name = simulation_name,
				mean_depth = depth,
				outdir = paste0(outdir, "/", simulation_name))

			write_sim_annotation(simulation_name, sig_activities, sig_header, 
				sim_activity_file, sim_purity_file, sim_tumortype_file)

			sim_list <- c(sim_list, simulation_name)
		}
	}


    print("Simulation type 0b: two clusters")
	# signatures change in one cluster but not in the other"
	n_simulations = 5
	for (sim_id in 1:n_simulations) {
		sig_activities = list()

		meaningful_sig_list = c("SBS2+13", "SBS3", "SBS4", "SBS6", "SBS9", "SBS7")
		list[meaningful_sig1, meaningful_sig2] = sample(meaningful_sig_list, size = 2)
		
		# Signatures change in cluster 2, but not in cluster 1
		clonal_sigs <- sample_sigs_and_activities(meaningful_sig1, meaningful_sig2, sig1_range=c(0.4, 0.8))

		sig_activities[[1]] <- clonal_sigs
		sig_activities[[2]] <- sample_sigs_and_activities(meaningful_sig1, meaningful_sig2, sig1_range=c(0.2, 0.4))

		print("Sig activities")
		print(do.call(rbind,sig_activities))

		subclone1_ccf = runif(1, min=0.1, max=0.4)

		print("CCFs per cluster")
		print(c(1.0, subclone1_ccf))

		# cluster 1 and cluster 2 and two separate branches
		n_mut_subclone1 = as.integer(mut_per_sim * subclone1_ccf)

		n_mut_clonal = mut_per_sim - n_mut_subclone1

		print("Mutation counts per cluster")
		print(c(n_mut_clonal, n_mut_subclone1))

		for (depth in c(100, 1000)) {
			simulation_name = paste0("Simulation_two_clusters", sim_id, "_depth", depth)
			print(paste0("Generating simulation ",simulation_name))

			sim_data_all_clusters = generate_ccf_simulation(
				n_clusters = 2, 
				cluster_ccfs = c(1.0, subclone1_ccf), 
				n_mut_per_cluster = c(n_mut_clonal, n_mut_subclone1), 
				sig_activities = sig_activities,
				signature_def = signature_def,
				simulation_name = simulation_name,
				mean_depth = depth,
				outdir = paste0(outdir, "/", simulation_name))

			write_sim_annotation(simulation_name, sig_activities, sig_header, 
				sim_activity_file, sim_purity_file, sim_tumortype_file)
		
			sim_list <- c(sim_list, simulation_name)
		}
	}

    sim_list = c()
    print("Simulation type 1: branching")
	# signatures change in one cluster but not in the other"
	n_simulations = 1
	for (sim_id in 1:n_simulations) {
		sig_activities = list()

		meaningful_sig_list = c("SBS2+13", "SBS3", "SBS4", "SBS6", "SBS7", "SBS9")
		list[meaningful_sig1, meaningful_sig2] = sample(meaningful_sig_list, size = 2)
		
		# Signatures change in cluster 2, but not in cluster 1
		clonal_sigs <- sample_sigs_and_activities(meaningful_sig1, meaningful_sig2, sig1_range=c(0.4, 0.8))

		sig_activities[[1]] <- clonal_sigs
		sig_activities[[2]] <- clonal_sigs
		sig_activities[[3]] <- sample_sigs_and_activities(meaningful_sig1, meaningful_sig2, sig1_range=c(0.2, 0.4))

		print("Sig activities")
		print(do.call(rbind,sig_activities))

		subclone1_ccf = runif(1, min=0.1, max=0.4)
		subclone2_ccf = runif(1, min=subclone1_ccf+0.2, max=1 - subclone1_ccf) # CCF2 > CCF1

		print("CCFs per cluster")
		print(c(1.0, subclone1_ccf, subclone2_ccf))

		stopifnot(subclone2_ccf > subclone1_ccf)

		# cluster 1 and cluster 2 and two separate branches
		n_mut_subclone1 = as.integer(mut_per_sim * subclone1_ccf)
		n_mut_subclone2 = as.integer(mut_per_sim * subclone2_ccf)

		n_mut_clonal = mut_per_sim - n_mut_subclone1 - n_mut_subclone2

		print("Mutation counts per cluster")
		print(c(n_mut_clonal, n_mut_subclone1, n_mut_subclone2))

		for (depth in c(100, 1000)) {
			simulation_name = paste0("Simulation_branching", sim_id, "_depth", depth)
			print(paste0("Generating simulation ",simulation_name))

			sim_data_all_clusters = generate_ccf_simulation(
				n_clusters = 3, 
				cluster_ccfs = c(1.0, subclone1_ccf, subclone2_ccf), 
				n_mut_per_cluster = c(n_mut_clonal, n_mut_subclone1, n_mut_subclone2), 
				sig_activities = sig_activities,
				signature_def = signature_def,
				simulation_name = simulation_name,
				mean_depth = depth,
				outdir = paste0(outdir, "/", simulation_name))

			write_sim_annotation(simulation_name, sig_activities, sig_header, 
				sim_activity_file, sim_purity_file, sim_tumortype_file)
			
			sim_list <- c(sim_list, simulation_name)
		}
	}


	print("Simulation type 2: CNA")
	# Signatures change in both cluster 1 and cluster 2

	print("Simulation 2a: 10% mutations are affected, CNA+1")
	n_simulations = 1
	for (sim_id in 1:n_simulations) {
		sig_activities = list()

		meaningful_sig_list = c("SBS2+13", "SBS3", "SBS4", "SBS6", "SBS7", "SBS9")
		list[meaningful_sig1, meaningful_sig2] = sample(meaningful_sig_list, size = 2)
		
		for (i in 1:3) {
			sig_activities[[i]] <- sample_sigs_and_activities(meaningful_sig1, meaningful_sig2)
		}

		print("Sig activities")
		print(do.call(rbind,sig_activities))

		subclone1_ccf = runif(1, min=0.1, max=0.3)
		subclone2_ccf = runif(1, min=subclone1_ccf+0.1, max=1 - subclone1_ccf - 0.05) # CCF2 > CCF1

		print("CCFs per cluster")
		print(c(1.0, subclone1_ccf, subclone2_ccf))

		stopifnot(subclone2_ccf > subclone1_ccf)

		# cluster 1 and cluster 2 and two separate branches
		n_mut_subclone1 = as.integer(mut_per_sim * subclone1_ccf)
		n_mut_subclone2 = as.integer(mut_per_sim * subclone2_ccf)

		n_mut_clonal = mut_per_sim - n_mut_subclone1 - n_mut_subclone2

		print("Mutation counts per cluster")
		print(c(n_mut_clonal, n_mut_subclone1, n_mut_subclone2))

		 cluster_cna_info = list(
		  	list("fractions"= c(0.9, 0.05, 0.05), "mut_cn"= c(1, 1, 2), "total_cn"= c(2, 3, 3)), #clonal
		  	list("fractions"= c(0.9, 0.05, 0.05), "mut_cn"= c(1, 1, 1), "total_cn"= c(2, 3, 3)), #subclone 1
		  	list("fractions"= c(0.9, 0.05, 0.05), "mut_cn"= c(1, 1, 1), "total_cn"= c(2, 3, 3))) #subclone 2

		for (depth in c(100, 1000)) {
			simulation_name = paste0("Simulation_cna_plus", sim_id, "_depth", depth)
			print(paste0("Generating simulation ",simulation_name))

			sim_data_all_clusters = generate_ccf_simulation(
				n_clusters = 3, 
				cluster_ccfs = c(1.0, subclone1_ccf, subclone2_ccf), 
				n_mut_per_cluster = c(n_mut_clonal, n_mut_subclone1, n_mut_subclone2), 
				sig_activities = sig_activities,
				signature_def = signature_def,
				cluster_cna_info = cluster_cna_info,
				simulation_name = simulation_name,
				outdir = paste0(outdir, "/", simulation_name))

			write_sim_annotation(simulation_name, sig_activities, sig_header, 
				sim_activity_file, sim_purity_file, sim_tumortype_file)

			sim_list <- c(sim_list, simulation_name)
		}
	}


	# print("Simulation 2b: 10% mutations are affected, CNA-1")
	# n_simulations = 1
	# for (sim_id in 1:n_simulations) {
	# 	sig_activities = list()

	# 	meaningful_sig_list = c("SBS2+13", "SBS3", "SBS4", "SBS6", "SBS7", "SBS9")
	# 	list[meaningful_sig1, meaningful_sig2] = sample(meaningful_sig_list, size = 2)
		
	# 	for (i in 1:3) {
	# 		sig_activities[[i]] <- sample_sigs_and_activities(meaningful_sig1, meaningful_sig2)
	# 	}

	# 	print("Sig activities")
	# 	print(do.call(rbind,sig_activities))

	# 	subclone1_ccf = runif(1, min=0.1, max=0.3)
	# 	subclone2_ccf = runif(1, min=subclone1_ccf+0.1, max=1 - subclone1_ccf - 0.05) # CCF2 > CCF1

	# 	print("CCFs per cluster")
	# 	print(c(1.0, subclone1_ccf, subclone2_ccf))

	# 	stopifnot(subclone2_ccf > subclone1_ccf)

	# 	# cluster 1 and cluster 2 and two separate branches
	# 	n_mut_subclone1 = as.integer(mut_per_sim * subclone1_ccf)
	# 	n_mut_subclone2 = as.integer(mut_per_sim * subclone2_ccf)

	# 	n_mut_clonal = mut_per_sim - n_mut_subclone1 - n_mut_subclone2

	# 	print("Mutation counts per cluster")
	# 	print(c(n_mut_clonal, n_mut_subclone1, n_mut_subclone2))

	# 	 cluster_cna_info = list(
	# 	  	list("fractions"= c(0.9, 0.05, 0.05), "mut_cn"= c(1, 0, 1), "total_cn"= c(2, 1, 1)), #clonal
	# 	  	list("fractions"= c(0.9, 0.05, 0.05), "mut_cn"= c(1, 1, 1), "total_cn"= c(2, 1, 1)), #subclone 1
	# 	  	list("fractions"= c(0.9, 0.05, 0.05), "mut_cn"= c(1, 1,1), "total_cn"= c(2, 1, 1))) #subclone 2

	#for (depth in c(100, 1000)) {
		# 	simulation_name = paste0("Simulation_cna_minus", sim_id, "_depth", depth)
		#	print(paste0("Generating simulation ",simulation_name))

		# 	sim_data_all_clusters = generate_ccf_simulation(
		# 		n_clusters = 3, 
		# 		cluster_ccfs = c(1.0, subclone1_ccf, subclone2_ccf), 
		# 		n_mut_per_cluster = c(n_mut_clonal, n_mut_subclone1, n_mut_subclone2), 
		# 		sig_activities = sig_activities,
		# 		signature_def = signature_def,
		# 		cluster_cna_info = cluster_cna_info,
		# 		simulation_name = simulation_name,
		# 		outdir = paste0(outdir, "/", simulation_name))

		# 	write_sim_annotation(simulation_name, sig_activities, sig_header,
		#                    sim_activity_file, sim_purity_file, sim_tumortype_file)

		# 	sim_list <- c(sim_list, simulation_name)
		# }
	#}



	print("Simulation 3a: Violation of infinite site assumption with CCF1+CCF2")
	n_simulations = 1
	for (sim_id in 1:n_simulations) {
		sig_activities = list()

		meaningful_sig_list = c("SBS2+13", "SBS3", "SBS4", "SBS6", "SBS7", "SBS9")
		list[meaningful_sig1, meaningful_sig2] = sample(meaningful_sig_list, size = 2)
		
		# Signatures change in cluster 2, but not in cluster 1
		clonal_sigs <- sample_sigs_and_activities(meaningful_sig1, meaningful_sig2)
		sig_activities[[1]] <- clonal_sigs
		sig_activities[[2]] <- sample_sigs_and_activities(meaningful_sig1, meaningful_sig2)
		sig_activities[[3]] <- clonal_sigs
		# Small cluster of mutations that violate infinite cite assumption
		sig_activities[[4]] <- sig_activities[[1]]

		print("Sig activities")
		print(do.call(rbind,sig_activities))

		frac_inf_site = 0.03
		subclone1_ccf = runif(1, min=0.1, max=0.3)
		subclone2_ccf = runif(1, min=subclone1_ccf+0.1, max=1 - subclone1_ccf - frac_inf_site-0.1) # CCF2 > CCF1
		infinite_site_viol_ccf = subclone1_ccf + subclone2_ccf

		stopifnot(infinite_site_viol_ccf < 1.0)

		print("CCFs per cluster")
		print(c(1.0, subclone1_ccf, subclone2_ccf, infinite_site_viol_ccf))


		stopifnot(subclone2_ccf > subclone1_ccf)

		# cluster 1 and cluster 2 and two separate branches
		n_mut_subclone1 = as.integer(mut_per_sim * subclone1_ccf)
		n_mut_subclone2 = as.integer(mut_per_sim * subclone2_ccf)
		n_mut_inf_site_viol = mut_per_sim * 0.03 # 3% of mutations violte infinite site assumption

		n_mut_clonal = mut_per_sim - n_mut_subclone1 - n_mut_subclone2 - n_mut_inf_site_viol

		print("Mutation counts per cluster")
		print(c(n_mut_clonal, n_mut_subclone1, n_mut_subclone2, n_mut_inf_site_viol))

		for (depth in c(100, 1000)) {
			simulation_name = paste0("Simulation_inf_site_viol_plus", sim_id, "_depth", depth)
			print(paste0("Generating simulation ", simulation_name))

			sim_data_all_clusters = generate_ccf_simulation(
				n_clusters = 4, 
				cluster_ccfs = c(1.0, subclone1_ccf, subclone2_ccf, infinite_site_viol_ccf), 
				n_mut_per_cluster = c(n_mut_clonal, n_mut_subclone1, n_mut_subclone2, n_mut_inf_site_viol), 
				sig_activities = sig_activities,
				signature_def = signature_def,
				simulation_name = simulation_name,
				outdir = paste0(outdir, "/", simulation_name))

			write_sim_annotation(simulation_name, sig_activities, sig_header,
				sim_activity_file, sim_purity_file, sim_tumortype_file)

			sim_list <- c(sim_list, simulation_name)
		}
	}



	# print("Simulation 3b: Violation of infinite site assumption with CCF1-CCF2")
	# n_simulations = 1
	# for (sim_id in 1:n_simulations) {
	# 	sig_activities = list()

	# 	meaningful_sig_list = c("SBS2+13", "SBS3", "SBS4", "SBS6", "SBS7", "SBS9")
	# 	list[meaningful_sig1, meaningful_sig2] = sample(meaningful_sig_list, size = 2)
		
	# 	for (i in 1:3) {
	# 		sig_activities[[i]] <- sample_sigs_and_activities(meaningful_sig1, meaningful_sig2)
	# 	}
	# 	# Small cluster of mutations that violate infinite cite assumption
	# 	sig_activities[[4]] <- sig_activities[[1]]

	# 	print("Sig activities")
	# 	print(do.call(rbind,sig_activities))

	# 	frac_inf_site = 0.03
	# 	subclone1_ccf = runif(1, min=0.1, max=0.3)
	# 	subclone2_ccf = runif(1, min=subclone1_ccf+0.1, max=1 - subclone1_ccf - frac_inf_site-0.1) # CCF2 > CCF1

	# 	infinite_site_viol_ccf = subclone2_ccf - subclone1_ccf

	# 	print("CCFs per cluster")
	# 	print(c(1.0, subclone1_ccf, subclone2_ccf, infinite_site_viol_ccf))

	# 	stopifnot(subclone2_ccf > subclone1_ccf)
	# 	stopifnot(infinite_site_viol_ccf > 0)

	# 	# cluster 1 and cluster 2 and two separate branches
	# 	n_mut_subclone1 = as.integer(mut_per_sim * subclone1_ccf)
	# 	n_mut_subclone2 = as.integer(mut_per_sim * subclone2_ccf)
	# 	n_mut_inf_site_viol = mut_per_sim * frac_inf_site # 3% of mutations violte infinite site assumption

	# 	n_mut_clonal = mut_per_sim - n_mut_subclone1 - n_mut_subclone2 - n_mut_inf_site_viol

	# 	print("Mutation counts per cluster")
	# 	print(c(n_mut_clonal, n_mut_subclone1, n_mut_subclone2, n_mut_inf_site_viol))

	#for (depth in c(100, 1000)) {
		# 	simulation_name = paste0("Simulation_inf_site_viol_minus", sim_id, "_depth", depth)
	#		print(paste0("Generating simulation ",simulation_name))
		# 	sim_data_all_clusters = generate_ccf_simulation(
		# 		n_clusters = 4, 
		# 		cluster_ccfs = c(1.0, subclone1_ccf, subclone2_ccf, infinite_site_viol_ccf), 
		# 		n_mut_per_cluster = c(n_mut_clonal, n_mut_subclone1, n_mut_subclone2, n_mut_inf_site_viol), 
		# 		sig_activities = sig_activities,
		# 		signature_def = signature_def,
		# 		simulation_name = simulation_name,
		# 		outdir = paste0(outdir, "/", simulation_name))

		# 	write_sim_annotation(simulation_name, sig_activities, sig_header,
		#                    sim_activity_file, sim_purity_file, sim_tumortype_file)

		# 	sim_list <- c(sim_list, simulation_name)
	#}
	# }



	print("Created simulations:")
	print(sim_list)
	return(sim_list)
}