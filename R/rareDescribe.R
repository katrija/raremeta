rareDescribe <- function(ai, bi, ci, di, n1i, n2i,
                         data){

  # check if all arguments are given - function needs ai, n1i, ci, n2i or ai, bi, ci, di
  # if it is given, check if it is integer
  if (missing(ai)) {
    stop("ai must be specified.")
  }

  if (missing(bi) & missing(n1i)) {
    stop("bi or n1i must be specified.")
  }

  if (missing(ci)) {
    stop("ci must be specified.")
  }

  if (missing(di) & missing(n2i)) {
    stop("di or n2i must be specified.")
  }

  mf <- match.call()

  mf.ai <- mf[[match("ai", names(mf))]]
  mf.bi <- mf[[match("bi", names(mf))]]
  mf.ci <- mf[[match("ci", names(mf))]]
  mf.di <- mf[[match("di", names(mf))]]
  mf.n1i <- mf[[match("n1i", names(mf))]]
  mf.n2i <- mf[[match("n2i", names(mf))]]

  ai <- eval(mf.ai, data)
  bi <- eval(mf.bi, data)
  ci <- eval(mf.ci, data)
  di <- eval(mf.di, data)
  n1i <- eval(mf.n1i, data)
  n2i <- eval(mf.n2i, data)


  # check if data is a data frame. If not: stop
  if (!is.data.frame(data)) {
    stop("Argument 'data' must be a data frame.")
  }

  # check if all variables are not negative
  if (any(c(ai, bi, ci, di) < 0, na.rm = TRUE)) {
    stop("Arguments must not be negative.")
  }

  # check if length of all variables is the same
  if (!equalLength(ai, bi, ci, di, n1i, n2i)) {
    stop("Data vectors must be of the same length")
  }

  # in case bi and n1i (di and n2i) are given: check whether ai + bi = n1i (ci + di = n2i)
  if (!is.null(ai) & !is.null(bi) & !is.null(n1i)) {
    if (any(ai + bi != n1i)) {
      stop("ai and bi must add up to n1i.")
    }
  }

  if (!is.null(ci) & !is.null(di) & !is.null(n2i)) {
    if (any(ci + di != n2i)) {
      stop("ci and di must add up to n2i.")
    }
  }

  # in case n1i (n2i) is given, but bi (di) not, calculate bi (di)
  if (!is.null(n1i) & is.null(bi)) {
    bi <- n1i - ai
  }
  if (!is.null(n2i) & is.null(di)) {
    di <- n2i - ci
  }


  # give warning if there are missing values
  if (any(is.na(c(ai, bi, ci, di)))) {
    k1NA <- unique(c(which(is.na(ai)), which(is.na(bi))))
    k2NA <- unique(c(which(is.na(ci)), which(is.na(di))))

    if (length(k1NA) != 0) {
      warning(paste0("There are missing values in group 1 in studies: ", k1NA))
    }

    if (length(k2NA) != 0) {
      warning(paste0("There are missing values in group 2 in studies: ", k2NA))
    }
  }


  nEventsT <- sum(ai, na.rm = TRUE)
  nEventsC <- sum(ci, na.rm = TRUE)

  descEvents <- paste("There are", nEventsT, "events in the treatment group and", nEventsC, "in the control group.")

  return(descEvents)
}
