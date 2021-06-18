
# library(TrackSig)
# library(ggplot2)
# #
# vcfFile = system.file(package = "TrackSig", "extdata/Example.vcf")
# cnaFile = system.file(package = "TrackSig", "extdata/Example_cna.txt")
# purity = 1
# refGenome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
# binSize = 100
# context = generateContext(c("CG", "TA"))
#
# vcf <- parseVcfFile(vcfFile)
# cnaGR <- parseCnaFile(cnaFile)
# vcf <- annotateCn(vcf, cnaGR)
# vcaf <- getVcaf(vcf, purity, cnaGR, refGenome)
# vcaf <- vcafConstruction(vcf, refGenome)
# vcaf <- getTrinuc(vcaf, refGenome)
# binCounts <- getBinCounts(vcaf, binSize, context)
#
#
# counts <- vcfToCounts(vcfFile, cnaFile = NULL, purity = 1, binSize = 100,
#                         context = generateContext(c("CG", "TA")),
#                         nCutoff = 10000, verbose = F,
#                         refGenome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
#
# nMut <- dim(vcaf)[1]
# assertthat::assert_that(nMut >= binSize, msg = "number of mutations may not be less than specified binSize\n")
# assertthat::assert_that(dim(unique(vcaf[c("ref","alt","mutType")]))[1] <= dim(context)[1], msg = sprintf("too many mutation types (%s) for context (%s)\n",
#                                                                                                          dim(unique(vcaf[c("ref","alt","mutType")]))[1],  dim(context)[1]) )
#
# # populate all possible bins, up to one possible remaining partial bin
# nFullBins <- floor( (nMut / binSize) )
# sizePartialBin <- round(((nMut / binSize) %% 1) * binSize)
#
# vcaf <- vcaf %>%
#   dplyr::mutate(chr_pos = as.numeric(paste(.data$chr, .data$pos, sep = "."))) %>%
#   dplyr::arrange(chr_pos)
#
# vcaf$bin <- c(rep(1:nFullBins, each = binSize), rep( (nFullBins + 1) , times = sizePartialBin))
#
# # get count of each mutation type for each bin
# vcaf %>%
#   dplyr::mutate(cat = paste(.data$ref, .data$alt, .data$mutType, sep = "_")) %>%
#   dplyr::group_by(.data$bin, .data$cat) %>%
#   dplyr::summarize(sum = sum(.data$bin)) %>%
#   tidyr::spread(cat, sum) -> binCounts
#
# # replace NAs with 0
# binCounts[is.na(binCounts)] <- 0
#
# context <- generateContext(c("CG", "TA"))
