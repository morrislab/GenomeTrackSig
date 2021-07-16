## Contains functions used to aggregate PCAWG mutation counts data and represent counts data
## as 'bins,' each bin containing a user-defined number of mutations and spanning a variable
## segment of the genome.

#' Summarize mutation counts across all samples of interest
#'
#' @description
#' \code{poolSamples} will take an input CSV file of mutation counts, and
#' determine a trajectory over the genome for the sample, based on changepoints
#' found using the PELT segmentation algorithm.
#'
#' @param archivePath file path to folder containing all cancer samples
#' @param typesPath file path to dataframe (csv) containing file IDs and
#' corresponding cancer type
#' @param cancerType string denoting cancer type for which all
#' samples will be pooled
#'
#' @return dataframe of pooled mutation counts, csv file to working directory.
#' @export

poolSamples <- function(archivePath, typesPath, cancerType) {
  seqnames <- strand <- NULL
  # read in dataframe of file IDs with corresponding cancer types
  types <- readr::read_csv(typesPath,
                           col_names = c('type', 'guid'))
  types <- types[2:nrow(types), ]

  # list of sample filenames for desired cancer type
  get_files <- c(types$guid[types$type == cancerType])

  # initialize counts dataframe
  master <- readr::read_csv(as.character(paste(archivePath, "/", get_files[1], ".MBcounts.csv", sep = "")),
                            col_types = readr::cols("seqnames" = "c",
                                                    "strand" = "c",
                                                    .default = "d"))

  # add counts from remaining samples into master file and delete remaining files from memory
  for (i in 2:4) {
    temp <- readr::read_csv(as.character(paste(archivePath, "/", get_files[i], ".MBcounts.csv", sep = "")),
                            col_types = readr::cols(seqnames = "c",
                                                    strand = "c",
                                                    .default = "d"))
    master[, 6:101] <- master[, 6:101] + temp[, 6:101]
    base::remove(temp)
  }

  # make data type consistent in seqnames column
  # split seqnames column into start chrom and end chrom
  # remove unnecessary variables
  master <- master %>%
    dplyr::mutate(seqnames = dplyr::case_when(seqnames == "X" ~ 23,
                                              seqnames == "Y" ~ 24,
                                              TRUE ~ as.numeric(seqnames)),
                  start_chrom = seqnames,
                  end_chrom = seqnames) %>%
    dplyr::select(-seqnames, -strand) %>%
    base::subset(select = c(100,1,101,2:99))

  # save output into csv and return dataframe
  readr::write_csv(master, paste(cancerType, "_pooled.csv", sep=""))
  return (master)
}

getBinNumber <- function(master, binSize) {
  counts <- data.frame()
  for (i in c(1:24)) {
    # find mutation counts in each row
    master_subset <- master %>%
      dplyr::filter(start_chrom==i)
    master_subset <- master_subset %>%
      dplyr::mutate(rowsum = rowSums(master_subset[,6:101]),
                    bin = 1) %>%
      base::subset(select = c(1:5,102,103,6:101))

    bin <- 1
    sums <- 0

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

  for (i in c(2:24)) {
    counts$bin[counts$start_chrom==i] <- counts$bin[counts$start_chrom==i] + max(counts$bin[counts$start_chrom==i-1])
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
    dplyr::summarize_at(dplyr::vars(C_A_ACA:T_G_TTT), sum) %>%
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
  for (i in c(1:24)) {
    master_subset <- master %>%
      dplyr::filter(start_chrom==i)
    # bin the counts for each chromosome
    counts_subset <- binningNmut(master_subset, binSize)
    # append binned chromosome data to final dataframe of bin counts
    counts <- rbind(counts, counts_subset)
  }
  # define bins in relation to whole genome rather than whole chromosome
  counts$bin <- c(1:nrow(counts))
  return (counts)
}
