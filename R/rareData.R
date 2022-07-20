#' Constructor function for class rareData
#'
#' @param x a list containing the elements required for objects of the class rareData
#'
#' @return an object of class rareData
new_rareData <- function(x = list()) {
  stopifnot(is.list(x))

  structure(x,
            class = "rareData"
  )
}

#' validator function for class rareData
#'
#' @param x an object of class rareData
#'
#' @return an object of class rareData
validate_rareData <- function(x) {
  elements <- names(x)

  if (any(!(c(
    "k", "ksz", "k1sz", "k2sz", "kdz", "n1", "n2", "n", "nratio",
    "rf1", "rf2", "rf"
  ) %in% elements))) {
    stop("some list elements are missing.")
  }

  x
}

#' Helper function to generate objects of class rareData
#'
#' @param x a list containing the elements required for objects of the class rareData
#'
#' @return an object of class rareData
rareData <- function(x = list()) {
  validate_rareData(new_rareData(x))
}
