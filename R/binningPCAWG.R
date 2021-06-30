
## \code{poolSamples} Generate dataframe summarizing mutation counts across all samples of a
## given cancer type.
##
## @param archivePath file path to folder containing all cancer samples
## @param typesPath file path to dataframe (csv) containing file IDs and corresponding cancer type
## @param cancerType string denoting cancer type for which all samples will be pooled
##
## @examples
## pooled_breastcancer_counts <- poolSamples(archivePath = '~/Desktop/archive',
##                                           typesPath = '~/Desktop/pcawg_cancer_types.csv',
##                                           cancerType = 'Breast-AdenoCA')
##
## @name poolSamples

poolSamples <- function(archivePath, typesPath, cancerType) {
  # read in dataframe of file IDs with corresponding cancer types
  types <- readr::read_csv(typesPath,
                           col_names = c('type', 'guid'))
  types <- types[2:nrow(types), ]

  # list of sample filenames for desired cancer type
  get_files <- c(types$guid[types$type == cancerType])
  # get_files <- c()
  # for (i in nrow(types)) {
  #   if (types$type[i] %in% cancerType) {
  #     get_files <- c(get_files, types$guid[i])
  #   }
  # }

  # initialize counts dataframe
  master <- readr::read_csv(as.character(paste(archivePath, "/", get_files[1], ".MBcounts.csv", sep = "")),
                            col_types = cols("seqnames" = "c",
                                             "strand" = "c",
                                             .default = "d"))

  # add counts from remaining samples into master file and delete remaining files from memory
  for (i in 2:length(get_files)) {
    temp <- readr::read_csv(as.character(paste(archivePath, "/", get_files[i], ".MBcounts.csv", sep = "")),
                            col_types = cols(seqnames = "c",
                                             strand = "c",
                                             .default = "d"))
    master[, 6:101] <- master[, 6:101] + temp[, 6:101]
    base::remove(temp)
  }
  readr::write_csv(master, paste(cancerType, "_pooled.csv", sep=""))
  return (master)
}

## \code{binningNmut} Group mutation counts dataframe into bins with specified # of mutations per bin.
##
## @param path file path to sample counts dataframe
## @param binSize desired number of mutations per bin (should be >= 100 for optimal results)
##
## @examples
## binCounts <- binningNmut(path = '~/Desktop/Breast-AdenoCA_pooled.csv',
##                          binSize = 100)
##
## @name binningNmut

binningNmut <- function(path, binSize) {
  # read in counts dataframe and organize descriptive columns
  counts <- readr::read_csv(path,
                            col_types = readr::cols("seqnames" = "c",
                                             "strand" = "c",
                                             .default = "d")) %>%
    dplyr::mutate(seqnames = dplyr::case_when(seqnames == "X" ~ 23,
                                        seqnames == "Y" ~ 24,
                                        TRUE ~ as.numeric(seqnames)),
                   start_chrom = seqnames,
                   end_chrom = seqnames) %>%
    dplyr::select(-seqnames, -strand) %>%
    base::subset(select = c(100,1,101,2:99))

  # iteratively collapse rows of dataframe into bins with ~binSize total mutations per row
  for (i in 1:nrow(counts)) {
    row_sums <- base::sum(base::colSums(counts[i,6:101]))
    if (i==nrow(counts)) {break}
    while (row_sums < binSize) {
      counts[i, 5:101] <- as.list(base::colSums(counts[c(i,i+1),c(5:101)]))
      counts[i,3:4] <- counts[i+1,3:4]
      row_sums <- base::sum(base::colSums(counts[i,6:101]))
      counts <- counts[-c(i+1), ]
      if (i == nrow(counts)-1) {break}
    }
  }

  # add in additional descriptive variables and sort dataframe in descending order (for plotting)
  counts <- counts %>%
    dplyr::mutate(bin = rep(1:nrow(counts)),
                  chr_pos = as.numeric(paste(counts$start_chrom, counts$start, sep = '.'))) %>%
    base::subset(select = c(1:5,103,102,6:101)) %>%
    dplyr::arrange(dplyr::desc(bin))

  return (counts)
}

