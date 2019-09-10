# TrackSig.R
# Defines main functions for user to interact with package TrackSig.
# Author: Cait Harrigan


detectActiveSignatures <- function(){

  # return list of active signatures in sample, whether by matching per-cancer-type to provided data,
  # or fitting all counts by EM. If not using this function, must provide active signatures per sample

  NULL
}

#' \code{TrackSig} Take an input vcf file and annotation and generate the counts data.
#' Create all plotting output that compute_signatures_for_all_examples does.
#'
#' @rdname TrackSig
#' @name TrackSig
#'
#' @param vcfFile path to variant calling format (vcf) file
#' @param cnaFile path to copy number abberation (cna) file
#' @param purityFile path to sample purity file
#' @param sampleID name to call sample. If none provided, name will be automatically drawn from the provided vcf file name.
#' @param saveIntermediate boolean whether to save intermediate results (mutation types)
#'
#'
#' activeInSample is list used to subset refrenceSignatures
#'
#' @export

TrackSig <- function(vcfFile,
                     cnaFile = NULL,
                     purity = NULL,
                     activeInSample = c("SBS1", "SBS5"),
                     sampleID = NULL,
                     refrenceSignatures = alex,
                     scoreMethod = "SigFreq",
                     binSize = 100,
                     desiredMinSegLen = NULL,
                     refGenome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19) {

  # input checking

  assertthat::assert_that(grepl(".vcf$", vcfFile) | grepl(".txt$", vcfFile), msg = "Unsupported VCF file extension. Expected file type .vcf or .txt")

  assertthat::assert_that(scoreMethod %in% c("SigFreq", "Signature", "Frequency"),
  msg = "scoreMethod should be one of \"SigFreq\", \"Signature\", \"Frequency\". \n Please see documentation for more information on selecting a scoreMethod)")

  # TODO: activeSignatures %in% rownames(referenceSignatures) must be TRUE
  # TODO: length(activeInSample) >1 should be true, else no mixture to fit

  # take sampleID from file name if not provided
  if (is.null(sampleID)){

    if (grepl(".txt$", vcfFile)){
      sampleID <- strsplit( unlist(strsplit(vcfFile, "/"))[ length( strsplit(vcfFile, "/")[[1]] ) ] , ".txt")[[1]]
    }

    if (grepl(".vcf$", vcfFile)){
      sampleID <- strsplit( unlist(strsplit(vcfFile, "/"))[ length( strsplit(vcfFile, "/")[[1]] ) ] , ".vcf")[[1]]
    }else{stop("Failed setting sampleID. Please check input vcf file.")}

  }

  # TODO: other parameters non-default options
  data <- vcfToCounts(vcfFile, cnaFile, purity)
  vcaf <- data[[1]]
  countsPerBin <- data[[2]]

  assertthat::assert_that(all(rownames(countsPerBin) == rownames(refrenceSignatures)), msg = "Mutation type counts failed.")

  # subset refrenceSignatures with activeInSample
  refrenceSignatures <- refrenceSignatures[activeInSample]

  if ( any(rowSums(countsPerBin)[rowSums(refrenceSignatures) == 0] != 0) ) {
    print(sprintf("Error in sample %s: Some mutation types have probability 0 under the model, but their count is non-zero. This count vector is impossible under the model.", sampleID))
  }

  # compute results
  # TODO: other parameters non-default options
  trajectory <- find_changepoints_pelt(countsPerBin, refrenceSignatures, vcaf, scoreMethod, binSize, desiredMinSegLen)
  changepoints <- trajectory[[1]]
  mixtures <- trajectory[[2]]


  # side effect: plot
  tryCatch({
            plot_name <- paste0(sampleID, " Signature Trajectory")
            binned_phis <- aggregate(vcaf$phi, by = list(vcaf$binAssignment), FUN = sum)$x / binSize
            mark_cp <- !is.null(changepoints)
            print(plot_signatures_real_scale(mixtures * 100, plot_name=plot_name, phis = binned_phis, mark_change_points=mark_cp,
                                       change_points=changepoints, transition_points = NULL, save = F)[[1]])
           },
           warning = function(w){w},
           error = function(e){print("Error: failed to plot signature trajectory")}
          )


  return (NULL)
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



# [END]
