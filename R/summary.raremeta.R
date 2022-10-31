#' Summary Method for raremeta objects.
#'
#' Summary method for objects of class "raremeta".
#'
#' @param object an object of class raremeta
#' @param ... other arguments.
#'
#' @return An object of class raremeta
#' @export
summary.raremeta <- function(object, ...){

  if (!is.element("raremeta", class(object))){
    stop("object must be an object of class 'raremeta'.")
  }

  return(object)
}
