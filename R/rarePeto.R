#' Conduct a meta-analysis using Peto's method
#'
#' Function to conduct a meta-analysis.
#' Effect size is the log odds ratio estimated from event counts in form of
#' 2x2 contingency tables using Peto's method (see Yusuf et al., 1985).
#'
#' @param x an object of class `"rareData"`.
#' @param ai data frame column to specify the number of events in group 1 (i.e., the treatment group).
#' @param bi data frame column to specify the number of non-events in group 1 (i.e., the treatment group).
#' @param ci data frame column to specify the number of events in group 2 (i.e., the control group).
#' @param di data frame column to specify number of non-events in group 2 (i.e., the control group).
#' @param n1i data frame column to specify the sample sizes in group 1 (i.e., the treatment group).
#' @param n2i data frame column to specify the sample sizes in group 2 (i.e., the control group).
#' @param data data frame.
#' @param level numeric nbetween 0 and 100 specifying the confidence interval
#' level (the default is 95)
#' @param digits integer specifying the number of decimal places to which the printed results
#' should be rounded (if unspecified, the default is 4).
#'
#' @details
#' # Details
#' ## Data input
#' The data input can be specified either through the arguments `ai`,`bi`,`ci`,`di`,`n1i`, and `n2i` (columns of the data frame `data`)
#' or through the argument `x`, which takes an object that results from applying the `rareDescribe()` function to the data
#' (i.e., the input for argument `x` must be an object of type `rareData`).
#' A `rareData` object can be produced from a data frame by applying the `rareDescribe()` function to it.
#' The `rareDescribe()` function pre-processes the data frame and stores the information required by the `rarePeto()` function in a list.
#' See `?rareDescribe` for more details.
#'
#' ## Effect size measure
#' The function is a meta-analytic method for estimating log odds ratios.
#' Data comes in form of event counts in form of 2x2 contingency tables.
#'
#' |                   | event | no event|
#' |:------------------|------:|--------:|
#' |group1 (treatment) | ai    | bi      |
#' |group2 (control)   | ci    | di      |
#'
#'
#'
#'
#' @return an object of class "raremeta".
#' The object is a list containing the following elements:
#' * `beta`, `b`: estimated effect size.
#' * `se`: standard error of the  estimator.
#' * `zval`: test statistics of the coefficients.
#' * `pval`: p-values corresponding to the test statistics.
#' * `ci.lb`: lower bound of the confidence intervals for the coefficients.
#' * `ci.ub`: upper bound of the confidence intervals for the coefficients.
#' * `k`: number of studies included in the analysis.
#' * `kdz`,`ksz`: number of double-zero and single-zero studies.
#' * `k1sz`, `k2sz`: number of single-zero studies where the zero is in group 1 or group 2.
#' * `ai`, `bi`, `ci`, `di`: original entries of the 2x2 tables for all studies.
#' * `ni`, `n1i`, `n2i`: original total and group sample sizes.
#' * ...
#'
#' @references
#' Yusuf, S., Peto, R., Lewis, J., Collins, R., & Sleight, P. (1985).
#' Beta blockade during and after myocardial infarction: an overview of
#' the randomized trials.
#' Progress in cardiovascular diseases, 27(5), 335-371.
#'
#' @export
#'
#' @examples
#' # load a data set
#' data(dat.nissen2007)
#' d <- dat.nissen2007
#'
#' # meta analysis using Peto's method
#' rarePeto(ai=miRosiglitazone, ci=miControl, n1i=nRosiglitazone, n2i=nControl, data=d)
#'
#' # same analysis with pre-processed data
#' x <- rareDescribe(ai=miRosiglitazone, ci=miControl, n1i=nRosiglitazone, n2i=nControl, data=d)
#' rarePeto(x)
#'
rarePeto <- function(x, ai, bi, ci, di, n1i, n2i, data, level = 95, digits = 4){


  ## argument checking ##

  # defining an object of class 'raredata' if raw data is put in
  if(missing(x)){
    if(missing(data))stop("Data must be specified via the argument 'x' or the argument 'data'.")

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

    x <- rareDescribe(ai=ai, bi=bi, ci=ci, di=di, n1i=n1i, n2i=n2i, data=data)
  }

  # check if x is an object of class rareData
  if(!inherits(x,"rareData")){
    stop("x must be an object of class 'rareData'. See ?rareDescribe for more details.")
  }

  # check if digits argument is valid
  if(length(digits) != 1 || digits%%1 != 0 || digits < 0){
    stop("'digits' must be an integer of length 1.")
  }

  # check if level argument is valid:
  if(!is.numeric(level) || length(level) > 1 || level < 0 || level > 100){
    stop("level must be a scalar between 0 and 100.")
  }

  mf <- match.call()

  # extract counts and sample sizes
  ai  <- x$ai
  bi  <- x$bi
  ci  <- x$ci
  di  <- x$di
  n1i <- x$n1i
  n2i <- x$n2i
  ni  <- n1i + n2i

  # calculating log(OR) estimate via O-E statistic
  measure <- "logOR"

  Ai    <- as.numeric(ai + ci)
  Bi    <- as.numeric(bi + di)
  Vi    <- (Ai * Bi * n1i * n2i) / (ni^2 *(ni - 1))
  sumVi <- sum(Vi)

  quant <- stats::qnorm((100-level)/200)

  if(sumVi == 0){
    stop("The data does not allow for the application of this method.\n
          There was no occurence of one of the outcomes in all of the studies.\n
          See e.g. ?rareCC() for possible continuity corrections.")
  }

  beta  <- sum((ai-Ai*(n1i/ni)))/sumVi
  se    <- sqrt(1/sumVi)
  zval  <- beta / se
  pval  <- 2*stats::pnorm(abs(zval), lower.tail=FALSE)
  ci.lb <- beta + quant * se
  ci.ub <- beta - quant * se


  res   <- append(list(beta = beta, b = beta, se = se, zval = zval,
                       pval = pval, ci.lb = ci.lb, ci.ub = ci.ub,
                       model = "rarePeto", measure = measure,
                       digits = digits),x)

  res <- raremeta(res)
  return(res)
}
