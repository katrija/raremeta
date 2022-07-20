#' Check if lengths of vectors are equal
#'
#' @param ... A collection of vectors, separated by a comma.
#'
#' @return A logical that indicates whether the lengths of the vectors are equal.
equalLength <- function(...) {
  z <- list(...)
  z <- z[sapply(z, function(x) length(x) > 0)]
  if (length(z) == 0L){
    return(TRUE)
  } else {
    lengthZ <- lengths(z)
    return(length(unique(lengthZ)) == 1L)
  }
}
