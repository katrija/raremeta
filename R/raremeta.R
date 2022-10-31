#' Constructor function for class raremeta
#'
#' @param x a list containing the elements required for objects of the class raremeta
#'
#' @return an object of class raremeta
new_raremeta <- function(x = list()) {
  stopifnot(is.list(x))

  structure(x,
            class = "raremeta"
  )
}

#' validator function for class raremeta
#'
#' @param x an object of class raremeta
#'
#' @return an object of class raremeta
validate_raremeta <- function(x) {
  elements <- names(x)

  if (any(!(c(
    "beta", "ci.lb", "ci.ub", "se", "vb", "p", "measure"
  ) %in% elements))) {
    stop("some list elements are missing.")
  }

  x
}

#' Helper function to generate objects of class raremeta
#'
#' @param x a list containing the elements required for objects of the class raremeta
#'
#' @return an object of class raremeta
raremeta <- function(x = list()) {
  validate_raremeta(new_raremeta(x))
}
