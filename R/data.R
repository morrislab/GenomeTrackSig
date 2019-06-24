#' example.vcf data.
#'
#' Description of the data:
#'
#' @format An example VCF file with 8 columns (only 5 columns contain data) and 3181 rows:
#' \describe{
#'   \item{#CHROM}{The name of the sequence (typically a chromosome) on which the variation is being called. This sequence is usually known as 'the reference sequence', i.e. the sequence against which the given sample varies.}
#'   \item{POS}{The 1-based position of the variation on the given sequence.}
#'   \item{REF}{The reference base (or bases in the case of an indel) at the given position on the given reference sequence.}
#'   \item{ALT}{The list of alternative alleles at this position.}
#'   \item{FILTER}{A flag indicating which of a given set of filters the variation has passed.}
#'   \item{INFO}{	An extensible list of key-value pairs (fields) describing the variation. See below for some common fields. Multiple fields are separated by semicolons with optional values in the format: "<key>=[,data]".}
#' }
#' @source \url{https://en.wikipedia.org/wiki/Variant_Call_Format}
#' @examples
#' \dontrun{
#' vcfFile <- system.file("extdata", "example.vcf", package="TrackSig")
#' myExample <- read.delim(vcfFile, stringsAsFactors = FALSE)
#' }
#' @docType data
#' @name example.vcf
NULL


