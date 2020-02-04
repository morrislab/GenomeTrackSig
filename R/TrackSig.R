# TrackSig.R
# Defines main functions for user to interact with package TrackSig.
# Author: Cait Harrigan

#' get binCounts from vcfToCounts
#' @param binCounts TODO
#' @param referenceSignatures TODO
#' @param threshold TODO
#' @param prior TODO
#'
#' @export
detectActiveSignatures <- function(binCounts, referenceSignatures = alex_merged, threshold = 0.05, prior = NULL){

  # return list of active signatures in sample, whether by matching per-cancer-type to provided data,
  # or fitting all counts by EM. If not using this function, must provide active signatures per sample

  counts <-  rowSums(binCounts)
  mixtures <- fitMixturesEM(counts, referenceSignatures, prior=prior)
  mixtures <- mixtures[mixtures >= threshold]

  return(names(mixtures))
}


#' \code{TrackSig} Take an input vcf file and annotation and generate the counts data.
#' Create all plotting output that compute_signatures_for_all_examples does.
#'
#' @rdname TrackSig
#' @name TrackSig
#'
#' @param vcfFile path to variant calling format (vcf) file
#' @param activeInSample list of signatures that are active. All listed signatures
#' must be present in the referenceSIgnatures dataframe.
#' @param cnaFile path to copy number abberation (cna) file
#' @param sampleID name to call sample. If none provided, name will be automatically drawn from the provided vcf file name.
#' @param referenceSignatures TODO
#' @param purity TODO
#' @param scoreMethod TODO
#' @param binSize TODO
#' @param desiredMinSegLen TODO
#' @param refGenome TODO
#'
#' @export

TrackSig <- function(vcfFile,
                     activeInSample,
                     cnaFile = NULL,
                     purity = NULL,
                     sampleID = NULL,
                     referenceSignatures = alex_merged,
                     scoreMethod = "SigFreq",
                     binSize = 100,
                     desiredMinSegLen = NULL,
                     refGenome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19) {

  # input checking

  assertthat::assert_that(grepl(".vcf$", vcfFile) | grepl(".txt$", vcfFile), msg = "Unsupported VCF file extension. Expected file type .vcf or .txt")

  assertthat::assert_that(scoreMethod %in% c("SigFreq", "Signature", "Frequency"),
  msg = "scoreMethod should be one of \"SigFreq\", \"Signature\", \"Frequency\". \n Please see documentation for more information on selecting a scoreMethod)")

  # TODO: activeSignatures %in% colnames(referenceSignatures) must be TRUE
  # TODO: binSize must be a natureal number >0
  # TODO: length(activeInSample) >1 should be true, else no mixture to fit
  # TODO: binSize has to make sense; positive, not larger than nMut, maybe throw warning if it's some ratio too large for low-resolution.
  # TODO: generateContext and mut types in referenceSignatures should make sense together.

  # take sampleID from file name if not provided
  if (is.null(sampleID)){

    if (grepl(".txt$", vcfFile)){
      sampleID <- strsplit( unlist(strsplit(vcfFile, "/"))[ length( strsplit(vcfFile, "/")[[1]] ) ] , ".txt")[[1]]
    }

    if (grepl(".vcf$", vcfFile)){
      sampleID <- strsplit( unlist(strsplit(vcfFile, "/"))[ length( strsplit(vcfFile, "/")[[1]] ) ] , ".vcf")[[1]]
    }else{stop("Failed setting sampleID. Please check input vcf file.")}

  }

  # TODO: geet context from supplied referenceSignatures
  context <- generateContext(c("CG", "TA"))

  # TODO: other parameters non-default options
  list[vcaf, countsPerBin] <- vcfToCounts(vcfFile = vcfFile, cnaFile = cnaFile,
                                          purity = purity, binSize = binSize,
                                          context = context, refGenome = refGenome)

  assertthat::assert_that(all(rownames(countsPerBin) %in%
                                rownames(referenceSignatures)),
                          msg = "Mutation type counts failed.")

  countsPerBin <- countsPerBin[rownames(referenceSignatures),]

  # subset referenceSignatures with activeInSample
  referenceSignatures <- referenceSignatures[activeInSample]

  if ( any(rowSums(countsPerBin)[rowSums(referenceSignatures) == 0] != 0) ) {
    print(sprintf("Error in sample %s: Some mutation types have probability 0 under the model, but their count is non-zero. This count vector is impossible under the model.", sampleID))
  }

  # compute results
  list[mixtures, changepoints] <- getChangepointsPELT(vcaf = vcaf,
                                                      countsPerBin = countsPerBin,
                                                      referenceSignatures = referenceSignatures,
                                                      scoreMethod = scoreMethod,
                                                      binSize = binSize,
                                                      desiredMinSegLen = desiredMinSegLen)

  # side effect: plot
  #plot <- NULL
  #tryCatch({

  #          binned_phis <- aggregate(vcaf$phi, by = list(vcaf$bin), FUN = mean)$x

  #          plot <- ( plotTrajectory(mixtures * 100, phis = binned_phis, changepoints, linearX = T, anmac = T)
  #                    + ggtitle(paste0(sampleID, " Signature Trajectory"))
  #                  )

  #          print(plot)

  #         },
  #         warning = function(w){w},
  #         error = function(e){print("Error: failed to plot signature trajectory")}
  #        )


  return (list(mixtures = mixtures, changepoints = changepoints, binData = vcaf))
}

# list unpacker util: used internally in package TrackSig
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
