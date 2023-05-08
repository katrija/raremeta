#' Descriptives for a meta-analysis of a rare event
#'
#' Function to compute descriptives based on the 2x2 tables of the individual studies
#' for meta-analytic rare event data.
#'
#' @param ai data frame column to specify the number of events in group 1 (i.e., the treatment group).
#' @param bi data frame column to specify the number of non-events in group 1 (i.e., the treatment group).
#' @param ci data frame column to specify the number of events in group 2 (i.e., the control group).
#' @param di data frame column to specify number of non-events in group 2 (i.e., the control group).
#' @param n1i data frame column to specify the sample sizes in group 1 (i.e., the treatment group).
#' @param n2i data frame column to specify the sample sizes in group 2 (i.e., the control group).
#' @param data data frame.
#'
#' @return An object of the class `"rareData"`, which contains descriptives of the meta-analytic data.
#' @export
#'
#' @examples
#' data <- data.frame(
#'   ai = c(0, 3, 2, 0),
#'   bi = c(20, 18, 15, 19),
#'   ci = c(1, 4, 0, 0),
#'   di = c(19, 17, 16, 20)
#' )
#'
#' desc <- rareDescribe(ai = ai, bi = bi, ci = ci, di = di, data = data)
#' print(desc)
rareDescribe <- function(ai, bi, ci, di, n1i, n2i,
                         data) {

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

  # check if all arguments are given - function needs ai, n1i, ci, n2i or ai, bi, ci, di
  # if it is given, check if it is integer
  if (is.null(ai)) {
    stop("ai must be specified.")
  }

  if (is.null(bi) & is.null(n1i)) {
    stop("bi or n1i must be specified.")
  }

  if (is.null(ci)) {
    stop("ci must be specified.")
  }

  if (is.null(di) & is.null(n2i)) {
    stop("di or n2i must be specified.")
  }


  # check if data is a data frame. If not: stop
  if (!is.data.frame(data)) {
    stop("Argument 'data' must be a data frame.")
  }

  # check if all variables are not negative
  if (any(c(ai, bi, ci, di, n1i, n2i) < 0, na.rm = TRUE)) {
    stop("Arguments must not be negative.")
  }

  if (!is.null(ai) & !is.null(n1i) & any(ai > n1i, na.rm = TRUE)) {
    stop("ai cannot be larger than n1i.")
  }

  if (!is.null(ci) & !is.null(n2i) & any(ci > n2i, na.rm = TRUE)) {
    stop("ci cannot be larger than n2i.")
  }

  # in case bi and n1i (di and n2i) are given: check whether ai + bi = n1i (ci + di = n2i)
  if (!is.null(ai) & !is.null(bi) & !is.null(n1i)) {
    if (any(ai + bi != n1i, na.rm = TRUE)) {
      stop("ai and bi must add up to n1i.")
    }
  }



  if (!is.null(ci) & !is.null(di) & !is.null(n2i)) {
    if (any(ci + di != n2i, na.rm = TRUE)) {
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

  #if (any(c(ai, bi, ci, di) %% 1 != 0, na.rm = TRUE)) {
  #  warning("Some values in ai, bi, ci or di are not integers. Please check whether this is intended.")
  #}

  # count all studies
  k <- length(ai)

  # label and count 0 and 00 studies:
  sz1studies <- (ai == 0) & !(ai == 0 & ci == 0)
  sz2studies <- (ci == 0) & !(ai == 0 & ci == 0)
  dzstudies  <- (ai == 0) & (ci == 0)
  szstudies  <- sz1studies | sz2studies

  kdz  <- sum(dzstudies, na.rm = TRUE)
  ksz  <- sum(szstudies, na.rm = TRUE)
  k1sz <- sum(sz1studies, na.rm = TRUE)
  k2sz <- sum(sz2studies, na.rm = TRUE)

  # vector for sample size treatment group:
  n1i <- ai + bi

  # vector for sample size control group:
  n2i <- ci + di

  # vector for total sample size:
  ni <- n1i + n2i

  # vector of sample size ratios:
  nratioi <- n1i / n2i

  # summary statistics:
  n1 <- c(
    mean(n1i, na.rm = TRUE),
    stats::sd(n1i, na.rm = TRUE), stats::median(n1i, na.rm = TRUE),
    stats::quantile(n1i, c(0.25, 0.75), na.rm = TRUE),
    min(n1i, na.rm = TRUE), max(n1i, na.rm = TRUE),
    stats::IQR(n1i, na.rm = TRUE)
  )
  names(n1) <- c("mean", "sd", "median", "q25", "q75", "min", "max", "IQR")

  n2 <- c(
    mean(n2i, na.rm = TRUE),
    stats::sd(n2i, na.rm = TRUE), stats::median(n2i, na.rm = TRUE),
    stats::quantile(n2i, c(0.25, 0.75), na.rm = TRUE),
    min(n2i, na.rm = TRUE), max(n2i, na.rm = TRUE),
    stats::IQR(n2i, na.rm = TRUE)
  )
  names(n2) <- c("mean", "sd", "median", "q25", "q75", "min", "max", "IQR")

  n <- c(
    mean(ni, na.rm = TRUE),
    stats::sd(ni, na.rm = TRUE), stats::median(ni, na.rm = TRUE),
    stats::quantile(ni, c(0.25, 0.75), na.rm = TRUE),
    min(ni, na.rm = TRUE), max(ni, na.rm = TRUE),
    stats::IQR(ni, na.rm = TRUE)
  )
  names(n) <- c("mean", "sd", "median", "q25", "q75", "min", "max", "IQR")

  nratio <- c(
    mean(nratioi, na.rm = TRUE),
    stats::sd(nratioi, na.rm = TRUE), stats::median(nratioi, na.rm = TRUE),
    stats::quantile(nratioi, c(0.25, 0.75), na.rm = TRUE),
    min(nratioi, na.rm = TRUE), max(nratioi, na.rm = TRUE),
    stats::IQR(nratioi, na.rm = TRUE)
  )
  names(nratio) <- c("mean", "sd", "median", "q25", "q75", "min", "max", "IQR")

  # relative frequencies:
  rf1i <- ai / n1i
  rf2i <- ci / n2i

  rfi <- (ai + ci) / ni

  rf1 <- c(
    mean(rf1i, na.rm = TRUE), stats::sd(rf1i, na.rm = TRUE),
    stats::median(rf1i, na.rm = TRUE),
    stats::quantile(rf1i, c(0.25, 0.75), na.rm = TRUE),
    min(rf1i, na.rm = TRUE), max(rf1i, na.rm = TRUE),
    stats::IQR(rf1i, na.rm = TRUE)
  )
  names(rf1) <- c("mean", "sd", "median", "q25", "q75", "min", "max", "IQR")

  rf2 <- c(
    mean(rf2i, na.rm = TRUE), stats::sd(rf2i, na.rm = TRUE),
    stats::median(rf2i, na.rm = TRUE),
    stats::quantile(rf2i, c(0.25, 0.75), na.rm = TRUE),
    min(rf2i, na.rm = TRUE), max(rf2i, na.rm = TRUE),
    stats::IQR(rf2i, na.rm = TRUE)
  )
  names(rf2) <- c("mean", "sd", "median", "q25", "q75", "min", "max", "IQR")

  rf <- c(
    mean(rfi, na.rm = TRUE), stats::sd(rfi, na.rm = TRUE),
    stats::median(rfi, na.rm = TRUE),
    stats::quantile(rfi, c(0.25, 0.75), na.rm = TRUE),
    min(rfi, na.rm = TRUE), max(rfi, na.rm = TRUE),
    stats::IQR(rfi, na.rm = TRUE)
  )
  names(rf) <- c("mean", "sd", "median", "q25", "q75", "min", "max", "IQR")

  krare <- sum(rfi <= 0.05 & rfi > 0.01, na.rm = TRUE)
  kveryrare <- sum(rfi <= 0.01, na.rm = TRUE)

  k1rare <- sum(rf1i <= 0.05 & rf1i > 0.01, na.rm = TRUE)
  k1veryrare <- sum(rf1i <= 0.01, na.rm = TRUE)

  k2rare <- sum(rf2i <= 0.05 & rf2i > 0.01, na.rm = TRUE)
  k2veryrare <- sum(rf2i <= 0.01, na.rm = TRUE)

  # CREATE LIST
  res <- list(
    ai = ai, bi = bi, ci = ci, di = di,
    n1i = n1i, n2i = n2i, ni = ni, nratioi = nratioi,
    k = k, kdz = kdz, ksz = ksz, k1sz = k1sz, k2sz = k2sz,
    dzstudies = dzstudies, szstudies = szstudies,
    sz1studies = sz1studies, sz2studies = sz2studies,
    n1 = n1, n2 = n2, n = n, nratio = nratio,
    rf1 = rf1, rf2 = rf2, rf = rf,
    krare = krare, kveryrare = kveryrare,
    k1rare = k1rare, k1veryrare = k1veryrare,
    k2rare = k2rare, k2veryrare = k2veryrare
  )

  # make res class rareData:
  res <- rareData(res)

  return(res)
}
