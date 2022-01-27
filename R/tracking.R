# Tracking.R
# Defines functions that carry out core capabilities for GenomeTrackSig
# Authors: Cait Harrigan and Caitlin Timmons


#' Get a list of signatures active in a sample
#'
#' @description
#' \code{detectActiveSignatures} determines what signatures are active above an
#' activity threshold, from a list of signature definintions. To do this, a
#' multinomial mixture model of signature activities is fit to the collection of
#' all mutations in the sample.
#'
#' @param df dataframe of mutation counts in 1 Mb segments across the genome
#' @param threshold minimum activity level that signature must have to be
#'   detected
#' @param prior prior on the likelihood of observing a given signature (must
#'   match signatures present in referenceSignatures)
#' @param binSize number of mutations per bin
#' @param referenceSignatures dataframe containing definitions of mutational
#'   signatures.
#' @param refGenome BSgenome to use as reference
#'
#' @return Names of signatures active in sample.
#' @export

detectActiveSignatures <- function(df,
                                   threshold = 0.05, prior = NULL, binSize,
                                   referenceSignatures = GenomeTrackSig:::cosmicV3,
                                   refGenome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19){

  # group mutations into bins with specified size
  binCounts <- binByChrom(df, binSize)

  # form counts matrix where each row is a mutation type and each column is a bin
  counts_subset <- binCounts[,7:102]
  counts <- data.table::transpose(counts_subset)
  colnames(counts) <- rownames(counts_subset)
  rownames(counts) <- colnames(counts_subset)

  # return list of active signatures in sample, whether by matching per-cancer-type to provided data,
  # or fitting all counts by EM. If not using this function, must provide active signatures per sample

  mixtures <- fitMixturesEM(counts, referenceSignatures, prior=prior)
  mixtures <- mixtures[mixtures >= threshold]
  mixtures <- sort(mixtures, decreasing = T)

  return(names(mixtures))
}


#' Determine a genomic signature profile.
#'
#'
#' @description
#' \code{TrackSig} will take an input CSV file of mutation counts and determine a genomic profile for the sample, based on changepoints found using the PELT segmentation algorithm.
#'
#' @param df dataframe of mutation counts from which we can fit activity profiles.
#' @param activeInSample list of signatures names to fit the exposures of. All listed signatures must be present in the referenceSignatures dataframe.
#' @param sampleID name to call sample. If none provided, name will be automatically drawn from the provided vcf file name.
#' @param referenceSignatures dataframe containing definitions of mutational signatures.
#' @param scoreMethod string to indicate what scoring method to apply when finding changepoints. Default = "Signature"
#' @param binSize number of mutations per bin
#' @param nCutoff maximum number of total mutations to consider (samples with more than nCutoff muations will be down-sampled)
#' @param desiredMinSegLen minimum number of mutations to include in a PELT segment (the desiredMinSegLen will be overridden if there are too few for accurate scoring)
#' @param refGenome BSgenome to use as reference
#'


TrackSig <- function(df,
                         binSize,
                         activeInSample,
                         sampleID = NULL,
                         referenceSignatures = GenomeTrackSig:::cosmicV3,
                         scoreMethod = "Signature",
                         nCutoff = 10000,
                         desiredMinSegLen = 1,
                         refGenome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19) {

  # input checking

  if(missing(activeInSample)){
    assertthat::assert_that(scoreMethod == "Frequency", msg = "When scoreMethod is not equal to \"Frequency\", activeInSample must be provided.")
  }else{
    assertthat::assert_that(all(activeInSample %in% colnames(referenceSignatures)))
  }


  # TODO: activeSignatures %in% colnames(referenceSignatures) must be TRUE
  # TODO: length(activeInSample) >1 should be true, else no mixture to fit
  # TODO: binSize has to make sense; positive, not larger than nMut, maybe throw warning if it's some ratio too large for low-resolution.
  # TODO: generateContext and mut types in referenceSignatures should make sense together.

  # TODO: get context from supplied referenceSignatures
  context <- generateContext(c("CG", "TA"))

  # TODO: other parameters non-default options
  countsPerBin <- binCounts <- counts_subset <- NULL

  binCounts <- df

  # form counts matrix where each row is a mutation type and each column is a bin
  counts_subset <- binCounts[,7:102]
  countsPerBin <- data.table::transpose(counts_subset)
  colnames(countsPerBin) <- rownames(counts_subset)
  rownames(countsPerBin) <- colnames(counts_subset)

  assertthat::assert_that(all(rownames(countsPerBin) %in%
                                rownames(referenceSignatures)),
                          msg = "Mutation type counts failed.")

  countsPerBin <- countsPerBin[rownames(referenceSignatures),,drop = FALSE]

  # subset referenceSignatures with activeInSample
  referenceSignatures <- referenceSignatures[activeInSample]

  if ( any(rowSums(countsPerBin)[rowSums(referenceSignatures) == 0] != 0) ) {
    print(sprintf("Error in sample %s: Some mutation types have probability 0 under the model, but their count is non-zero. This count vector is impossible under the model.", sampleID))
  }

  # compute results
  mixtures <- changepoints <- NULL
  list[mixtures, changepoints] <- getChangepointsPELT(binCounts = binCounts,
                                                      countsPerBin = countsPerBin,
                                                      referenceSignatures = referenceSignatures,
                                                      scoreMethod = scoreMethod,
                                                      binSize = binSize,
                                                      desiredMinSegLen = desiredMinSegLen)

  return (list(mixtures = mixtures, changepoints = changepoints, sampleID = sampleID, binData = binCounts))
}


# list unpacker util: used internally in package GenomeTrackSig
# source: https://stat.ethz.ch/pipermail/r-help/2004-June/053343.html

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

# TODO: function recommendBinSize()


# TODO: see http://adv-r.had.co.nz/S3.html for best practices
# constructor function for tracksig results. Returned by TrackSig when class=T
TS.trajectory <- function(sampleID = NULL, scoreMethod = NULL,
                          mixtures = NULL, changepoints = NULL,
                          binData = NULL, binSize = NULL){

  return(structure(list(), class = "TS.trajectory",
                   sampleID = sampleID, scoreMethod = scoreMethod,
                   mixtures = mixtures, changepoints = changepoints,
                   binData = binData, binSize = binSize
  )
  )
}

is.TS.trajectory <- function(x) inherits(x, "TS.trajectory")

# [END]

