if (file.exists("/mnt/raisin/")) {
  DIR= "/mnt/raisin/yulia/"
  READ_DIR = WRITE_DIR = DIR
} else {
  DIR <- "~/Documents/BCB430/"
  READ_DIR = WRITE_DIR = DIR
}

setwd(DIR)

#source("src/helper_functions.R")

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


load_signatures <- function() {
	sigs_to_merge <- list()
	sigs_to_merge[["SBS7"]] <- c("SBS7a", "SBS7b", "SBS7c", "SBS7d")
	sigs_to_merge[["SBS17"]] <- c("SBS17a", "SBS17b")
	sigs_to_merge[["SBS2+13"]] <- c("SBS2", "SBS13")
	sigs_to_merge[["SBS10"]] <- c("SBS10a", "SBS10b")

	names_trinucleotide <- read.table(paste0(DIR, "data/trinucleotide.txt"), stringsAsFactors = F)
	names_trinucleotide <- apply(names_trinucleotide, 1, function(x) { do.call("paste", c(as.list(x), sep = "_"))})

	alex <- read.csv(paste0(DIR, "data/sigProfiler_SBS_signatures.csv"), header = T, stringsAsFactors = FALSE)
    pcawg_trinucleotides <- paste(substr(alex[,1], 2, 2), substr(alex[,1], 5, 5), substr(alex[,1], 1, 3), sep="_")
    rownames(alex) <- pcawg_trinucleotides
    alex <- alex[,-1]
    alex <- alex[match(names_trinucleotide, pcawg_trinucleotides),]
	alex_merged <- t(merge_signatures(t(alex), sigs_to_merge))
  	return(alex)
}


# Generate simulations for mutations given clusters, their CCFs, CNAs and signatures.
generate_ccf_simulation <- function(n_clusters, cluster_ccfs, n_mut_per_cluster, sig_activities,
	cluster_cna_info = list(), mean_depth = 100,
	simulation_name = "simulation") {
	# Inputs:
	# n_clusters

	#cluster_ccfs

	#Clusters are ordered in the decreasing order of the CCF

	# n_mut_per_cluster
	# sig_activities -- list of signatures and their activities for each cluster
	#total_CN = 2, mut_alleles = 1
	# cluster_cna_info
	#fraction_cna_affected = 0
	# mean_depth

	normal_mut_alleles = 1
	normal_total_CN = 2

	stopifnot(length(cluster_ccfs) == n_clusters)
	stopifnot(length(n_mut_per_cluster) == n_clusters)

	alex <- load_signatures()
	trinucl_names <- rownames(alex)

	data_all_clusters <- c()
	for (cl in 1:n_clusters) {
		n_mut = n_mut_per_cluster[cl]
		cluster_ccf = cluster_ccfs[cl]

		cluster_sigs = sig_activities[[cl]]
		sig_names = names(cluster_sigs)
		stopifnot(sum(unlist(cluster_sigs)) == 1.)

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

				n_mut_with_cn_change = as.integer(n_mut * current_frac)
				prob = c(prob, rep(cluster_ccf * (current_mut_cn / current_total_cn),   n_mut_with_cn_change))
				mut_alleles = c(mut_alleles, rep(current_mut_cn, n_mut_with_cn_change))
				total_CN = c(total_CN, rep(current_total_cn, n_mut_with_cn_change))
			}
		}

		stopifnot(length(prob) == n_mut)

		# For each mutation in the cluster
		for (i in 1:n_mut) {
			# Sample number of variant alleles from a Binomial( d, CCF * (# mutant alleles / avg CN))

			alt <- rbinom(1, depths[i], prob[i])
			n_alt_alleles <- c(n_alt_alleles, alt)
			n_ref_alleles <- c(n_ref_alleles, depths[i] - alt)
		}


		print(mean(n_alt_alleles))
		print(sd(n_alt_alleles))

		# Sample mutation types for this cluster
		if (sum(sig_names %in% colnames(alex)) != length(sig_names)) {
			stop(paste("At least one of the signatures not found: ", toString(sig_names)))
		}

		mut_types <- c()
		# Sample mutations from each signature
		for (i in 1:length(cluster_sigs)) {
			sig_activity = cluster_sigs[[i]]
			mut_per_sig = as.integer(n_mut * sig_activity)
			sig_def <- alex[[names(cluster_sigs)[i]]]

			# Get counts over mutation types
			mut_type_counts <- rmultinom(1, mut_per_sig, unlist(sig_def))
			rownames(mut_type_counts) <- trinucl_names

			# Create a list of mutation types out of counts
			mut_types_from_sig <- unlist(sapply(1:length(mut_type_counts), function(ind) {rep(trinucl_names[ind], mut_type_counts[ind]) }))
			stopifnot(length(mut_types_from_sig) == mut_per_sig)

			mut_types <- c(mut_types, mut_types_from_sig)
		}

		# Shuffle list of mutation types
		mut_types <- sample(mut_types)
		stopifnot(length(mut_types) == n_mut)

		#  Report #variant, #reference, and avg CN for each mutation
		data <- data.frame(n_alt_alleles, n_ref_alleles, total_CN, mut_alleles, mut_types, stringsAsFactors = F)


		##################
		mut_ccfs <- data$n_alt_alleles / (data$n_alt_alleles + data$n_ref_alleles)
		print(mean(mut_ccfs))


		pdf(paste0("tmp/mut_ccf", cl, ".pdf"), width = 8, height=5)
		hist(mut_ccfs, breaks=200)
		dev.off()
		##################

		data_all_clusters <- rbind(data_all_clusters, data)
	}


	##################
	mut_ccfs <- data_all_clusters$n_alt_alleles / (data_all_clusters$n_alt_alleles + data_all_clusters$n_ref_alleles)

	pdf("tmp/mut_ccf.pdf", width = 8, height=5)
	hist(mut_ccfs, breaks=200)
	dev.off()
	##################

	save_as_vcf(data_all_clusters, simulation_name)
}


save_as_vcf <- function(data, filename) {
	n_mut = nrow(data)

	splitted_types <- t(data.frame(strsplit(data$mut_types, "_")))
	rownames(splitted_types) <- NULL
	ref <- splitted_types[,1]
	alt <- splitted_types[,2]
	tri <- splitted_types[,3]

	chrom <- paste0("chr", sample(1:22, n_mut, replace = TRUE))
	pos <- sample(10**5, n_mut)

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




create_simulation_set <- function(with_CNA = FALSE) {
	if (with_CNA == FLASE) {
	  # Basic simulation
	  sig_activities = list(
	  	list("SBS1" = 0.05, "SBS4" = 0.6, "SBS5" = 0.35),
	  	list("SBS1" = 0.05, "SBS4" = 0.2, "SBS5" = 0.75),
	  	list("SBS1" = 0.05, "SBS4" = 0.3, "SBS5" = 0.65))

	  generate_ccf_simulation(
	  	n_clusters = 3,
	  	cluster_ccfs = c(1.0, 0.5, 0.3),
	  	n_mut_per_cluster = c(10000, 5000, 3000),
	  	sig_activities = sig_activities,
	  	simulation_name = "basic_simulation")
	}
  else{
    # Simulation with mutations affected by CNA

	  cluster_cna_info = list(
	  	list("fractions"= c(1.0), "mut_cn"= c(1), "total_cn"= c(2)), #cluster 1
	  	list("fractions"= c(1.0), "mut_cn"= c(1), "total_cn"= c(2)), #cluster 2
	  	list("fractions"= c(0.75, 0.25), "mut_cn"= c(1, 2), "total_cn"= c(3, 3))) #cluster 3


	  generate_ccf_simulation(
	  	n_clusters = 3,
	  	cluster_ccfs = c(1.0, 0.8, 0.5),
	  	n_mut_per_cluster = c(1000, 500, 300),
	  	sig_activities = sig_activities,
	  	cluster_cna_info = cluster_cna_info,
	  	simulation_name = "cna_simulation")
  }
}
