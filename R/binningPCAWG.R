poolSamples <- function(archivePath, typesPath, cancerType) {
  types <- readr::read_csv(typesPath,
                           col_names = c('type', 'guid'))
  types <- types[2:nrow(types), ]
  get_files <- c(types$guid[types$type == cancerType])
  master <- readr::read_csv(as.character(paste(archivePath, "/", get_files[1], ".MBcounts.csv", sep = "")),
                            col_types = cols("seqnames" = "c",
                                             "strand" = "c",
                                             .default = "d"))
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


binningNmut <- function(path, binSize) {
  counts <- readr::read_csv(path,
                            col_types = cols("seqnames" = "c",
                                             "strand" = "c",
                                             .default = "d")) %>%
    dplyr::mutate(seqnames = dplyr::case_when(seqnames == "X" ~ 23,
                                        seqnames == "Y" ~ 24,
                                        TRUE ~ as.numeric(seqnames)),
                   start_chrom = seqnames,
                   end_chrom = seqnames) %>%
    dplyr::select(-seqnames, -strand) %>%
    base::subset(select = c(100,1,101,2:99))

  for (i in 1:nrow(counts)) {
    row_sums <- base::sum(base::colSums(counts[i,6:101]))
    if (i==nrow(counts)) {break}
    while (row_sums < 100) {
      counts[i, 5:101] <- as.list(base::colSums(counts[c(i,i+1),c(5:101)]))
      counts[i,3:4] <- counts[i+1,3:4]
      row_sums <- base::sum(base::colSums(counts[i,6:101]))
      counts <- counts[-c(i+1), ]
      if (i == nrow(counts)-1) {break}
    }
  }

  counts <- counts %>%
    dplyr::mutate(bin = rep(1:nrow(counts)),
                  chr_pos = as.numeric(paste(counts$start_chrom, counts$start, sep = '.'))) %>%
    base::subset(select = c(1:5,103,102,6:101)) %>%
    dplyr::arrange(dplyr::desc(bin))

  return (counts)
}

