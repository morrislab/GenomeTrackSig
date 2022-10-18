## \code{getBinNumber} Identify which bin each Mb region of the genome belongs to,
## depending on the desired bin size.
##
## @param master Un-binned data frame of mutation counts
## @param binSize Desired number of mutations in each bin
##
## @name getBinNumber

getBinNumber <- function(master, binSize) {
  counts <- data.frame()
  start_chrom <- NULL

  for (i in c(1:23)) {
    # find mutation counts in each row
    # bin sex chromosomes together
    if (i == 23) {
      master_subset <- master %>%
        dplyr::filter(start_chrom>=i)
      master_subset <- master_subset %>%
        dplyr::mutate(rowsum = rowSums(master_subset[,6:101]),
                      bin = 1) %>%
        base::subset(select = c(1:5,102,103,6:101))

      bin <- 1
      sums <- 0
      # add up counts in each row until desired bin size is reached
      for (i in 1:nrow(master_subset)) {
        sums <- sums + master_subset$rowsum[i]
        master_subset$bin[i] <- bin
        if (sums >= binSize) {
          sums <- 0
          bin <- bin + 1
        }
      }
      counts <- rbind(counts, master_subset)
    }
    else {
      # bin autosomal chromosomes individually
      master_subset <- master %>%
        dplyr::filter(start_chrom==i)
      master_subset <- master_subset %>%
        dplyr::mutate(rowsum = rowSums(master_subset[,6:101]),
                      bin = 1) %>%
        base::subset(select = c(1:5,102,103,6:101))

      bin <- 1
      sums <- 0
      # add up counts in each row until desired bin size is reached
      for (i in 1:nrow(master_subset)) {
        sums <- sums + master_subset$rowsum[i]
        master_subset$bin[i] <- bin
        if (sums >= binSize) {
          sums <- 0
          bin <- bin + 1
        }
      }
      counts <- rbind(counts, master_subset)
    }
  }
  # define bins in relation to the whole genome
  for (i in c(2:23)) {
    if (i == 23) {
      counts$bin[counts$start_chrom>=i] <- counts$bin[counts$start_chrom>=i] + max(counts$bin[counts$start_chrom==i-1])
    }
    else {
      counts$bin[counts$start_chrom==i] <- counts$bin[counts$start_chrom==i] + max(counts$bin[counts$start_chrom==i-1])
    }
  }

  return (counts)
}

## \code{binningNmut} Group mutation counts dataframe into bins with specified # of mutations per bin.
##
## @param master dataframe of mutation counts in your sample/s
## @param binSize desired number of mutations per bin (should be >= 100 for optimal results)
##
## @examples
## binCounts <- binningNmut(pooled_breastcancer_counts,
##                          binSize = 200)
##
## @name binningNmut

binningNmut <- function(master, binSize) {

  start_chrom <- end_chrom <- start <- end <- width <- C_A_ACA <- T_G_TTT <- NULL

  # find mutation counts in each row
  master <- master %>%
    dplyr::mutate(rowsum = rowSums(master[,6:101]),
                  bin = 1) %>%
    base::subset(select = c(1:5,102,103,6:101))

  # assign bin numbers to rows--bin number changes when the sum of all mutations in
  # the previous rows in that bin reaches or exceeds the desired bin size
  bin <- 1
  sums <- 0

  for (i in 1:nrow(master)) {
    sums <- sums + master$rowsum[i]
    master$bin[i] <- bin
    if (sums >= binSize) {
      sums <- 0
      bin <- bin + 1
    }
  }

  # find the total number of mutation counts for each type per bin
  dup <- data.table::copy(master)
  dup <- dup %>%
    dplyr::group_by(bin) %>%
    dplyr::summarize(across("C_A_ACA":"T_G_TTT", ~ sum(.x))) %>%
    dplyr::select(-bin)

  # find the genomic positions that demarcate the start and end of each bin
  master <- master %>%
    dplyr::group_by(bin) %>%
    dplyr::summarize(start_chrom = min(start_chrom),
                     start = min(start),
                     end_chrom = max(end_chrom),
                     end = max(end),
                     width = sum(width)) %>%
    base::cbind(dup)

  # return binned dataframe where each row is one bin
  return (master)
}

## \code{binByChrom} Perform binning on each chromosome separately.
## Binned counts for each chromosome are merged into dataframe of counts for full genome.
## This solution addresses the issue of bins that span multiple chromosomes without requiring the user to
## run TrackSig on each chromosome separately.
##
## @param master dataframe of mutation counts for your sample/s
## @param binSize desired number of mutations per bin (should be >= 100 for optimal results)
##
## @examples
## binCounts <- binByChrom(pooled_breastcancer_counts,
##                          binSize = 200)
##
## @name binByChrom

binByChrom <- function(master, binSize) {
  start_chrom <- NULL
  # initialize empty dataframe
  counts <- data.frame()
  # iterate through chromosomes
  for (i in c(1:23)) {
    if (i == 23) {
      master_subset <- master %>%
        dplyr::filter(start_chrom>=i)
      # bin the counts for each chromosome
      counts_subset <- binningNmut(master_subset, binSize)
      # append binned chromosome data to final dataframe of bin counts
      counts <- rbind(counts, counts_subset)
    }
    else {
      master_subset <- master %>%
        dplyr::filter(start_chrom==i)
      # bin the counts for each chromosome
      counts_subset <- binningNmut(master_subset, binSize)
      # append binned chromosome data to final dataframe of bin counts
      counts <- rbind(counts, counts_subset)
    }
  }
  # define bins in relation to whole genome rather than whole chromosome
  counts$bin <- c(1:nrow(counts))
  return (counts)
}
