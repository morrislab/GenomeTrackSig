# AUTHOR: Yulia Rubanova
# Modified for package TrackSig by Cait Harrigan



#' \code{generateContext} Generate a trinucleotide context from an alphabet. Note: this involves finding all three-member
#' permutations of the alphabet, which can be inconveinent for large alphabets. Nucleotides are assumed to be provided as complementary pairs,
#' where the first of each pair is used as the reference to build the context.
#'
#' @param alphabet list of pairs of characters to create combinations of as a mutation context type
#' @return data.frame containing all the possible trinucleotide contextes for a mutation in the supplied alphabet
#'
#' @examples
#' context <- TrackSig:::generateContext(c("CG", "TA"))
#' dim(context)
#' head(context)
#'
#' @rdname helper_functions
#' @name generateContext
#' @export

generateContext <- function(alphabet){

  if (any(nchar(alphabet) != 2)){
    stop("Alphabet is malformed. Please provide alphabet as a list of complementary pairs")
  }

  allpha <- unlist(strsplit(alphabet, split=NULL))
  nTypes <- (length(allpha) - 1) * length(allpha)^3 * 1/2

  context <- data.frame()

  for (i in seq(1, length(allpha), by = 2)){

    midRef <- allpha[i]
    rest <- setdiff(allpha, midRef)
    repSize <- length(allpha)^2 - length(allpha)

    midSet <- cbind(rep(midRef, length.out = repSize), rep(rest, length.out=repSize),
                    paste0(sort(rep(allpha, repSize)), rep(midRef, length.out = repSize), rep(allpha, repSize)))
    context <- rbind(context, midSet)
  }

  stopifnot( dim(context)[1] == nTypes )

  return (context)
}

#' Helper functions for \code{TrackSig}
#'
#' Non-exported functions called by \code{TrackSig} functions. \cr
#' Not intended for end-user use.
#'
#' @rdname helper_functions
#' @name helper_functions
NULL


#' \code{list} Succinct object assignment
#' @rdname helper_functions
#'
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

#' \code{get_values_from_list} Retrieve values by name from list
#' @rdname helper_functions
get_values_from_list <- function(list, name_of_value, FUN=NULL, default = NULL, concat_function="cbind")
{
  warning("Called a depricated function.")
  concat_function <- match.fun(concat_function)
  res <- c()
  for (i in 1:length(list))
  {
    if (!is.null(FUN))
    {
      FUN <- match.fun(FUN)
      value <- sapply(list[[i]][names(list[[i]]) == name_of_value], FUN)

      if (is.null(value[[1]]) & !is.null(default))
      {
        value[[1]] <- default
      }
      res <- concat_function(res, value)
    } else
    {
      res <- concat_function(res, unlist(list[[i]][names(list[[i]]) == name_of_value]))
    }
  }

  return(res)
}

#' \code{toVerticalMatrix} Convert a list or vector-like object to a vertial matrix
#' @rdname helper_functions
toVerticalMatrix <- function(L)
{
  if (is.vector(L))
    return(matrix(L, ncol=1))
  else
    return(as.matrix(L))
}

IgnoreVectorOrMatrix <- function(x, FUN)
{
  warning("Called a depricated function.")
  if (is.vector(x))
  {
    return(x)
  } else if (is.matrix(x) | is.data.frame(x)) {
    FUN <- match.fun(FUN)
    return(FUN(x))
  } else
  {
    stop(paste("Unknown type of data:", head(x)))
  }
}




#' \code{toHorizontalMatrix} <man content>
#' @rdname helper_functions
toHorizontalMatrix <- function(L){
  warning("Called a depricated function.")
  if (is.vector(L))
    return(matrix(L, nrow=1))
  else
    return(as.matrix(L))
}


#' \code{truncate_to_range} <man content>
#' @rdname helper_functions
truncate_to_range <- function(mixtures, range_) {
  warning("Called a depricated function.")
  min = range_[1]
  max = range_[2]

  x <- mixtures
  col_names <- as.numeric(colnames(x))
  to_leave <- which(col_names <= max+0.01 & col_names >= min-0.01)

  x2 <- x[,to_leave, drop=F]
  colnames(x2) <- col_names[to_leave]
  return(list(x2,to_leave))
}






#[END]
