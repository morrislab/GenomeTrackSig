# AUTHOR: Yulia Rubanova
# Modified for package TrackSig by Cait Harrigan






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
