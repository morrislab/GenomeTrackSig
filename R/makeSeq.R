# makeSeq.R

#' \code{makeSeq} output a random mRNA sequence.
#'
#' \code{makeSeq} outputs a string that contains a random sequence of characters (A/C/G/T by default). By default the first three letters are "ATG" and the last three letters are a stop-codon ("TAA", "TAG" or "TGA"). The number of characters in the output is three times the \code{len} argument, i.e. it condains \code{len} codons
#'
#' @param len integer number of codons in the output sequence.
#' @param useInit  (logical)  make the first codon an initiation codon. Default is \code{TRUE}.
#' @param useStop  (logical)  make the last codon a stop codon. Default is \code{TRUE}.
#' @param alphabet (character vector)  the elements that are randomly concatenated. Default is A/C/G/T.
#' @param p (numeric vector)  Probabilites for each character to occur. Default is equiprobable.
#' @param seed  (numeric)  the seed for the RNG. Default is \code{Sys.time()}.

#' @return A string.
#' @examples
#' makeSeq(7)
#' makeSeq(7, p = c(0.2, 0.4, 0.4, 0.2), seed = 112358)
#' gsub("T", "U", makeSeq(7)) # for RNA
#' @export
makeSeq <- function(len,
                    useInit = TRUE,
                    useStop = TRUE,
                    alphabet = c("A", "C", "G", "T"),
                    p = rep(0.25, 4),
                    seed = Sys.time()) {
  set.seed(seed)

  if (useInit) {
    myInit <- "ATG"
    len <- len - 1
  } else {
    myInit <- NULL
  }

  if (useStop) {
    myStop <- sample(c("TAA", "TAG", "TGA"), 1)
    len <- len - 1
  } else {
    myStop <- NULL
  }

  v <- c(myInit,
         sample(alphabet, len * 3, prob = p, replace = TRUE),
         myStop)

  return(paste(v, collapse = ""))
}

# [END]
