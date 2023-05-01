#' Conduct a meta-analysis using Peto's method
#'
#' Function to conduct a meta-analysis.
#' Effect size is the log odds ratio estimated from event counts in form of
#' 2x2 contingency tables using Peto's method (see Yusuf et al., 1985).
#'
#' @param x an object of class "raredata" (see `?rareDescribe` for more information)
#' @param ai Data frame column to specify the number of events in group 1 (i.e., the treatment group).
#' @param bi Data frame column to specify the number of non-events in group 1 (i.e., the treatment group).
#' @param ci Data frame column to specify the number of events in group 2 (i.e., the control group).
#' @param di Data frame column to specify number of non-events in group 2 (i.e., the control group).
#' @param n1i Data frame column to specify the sample sizes in group 1 (i.e., the treatment group).
#' @param n2i Data frame column to specify the sample sizes in group 2 (i.e., the control group).
#' @param data Data frame.
#' @param level numeric nbetween 0 and 100 specifying the confidence interval
#' level (the default is 95)
#' @param digits integer specifying the number of decimal places to which the printed results
#' should be rounded (if unspecified, the default is 4).
#'
#' @details
#' # Details
#' ## Data input
<<<<<<< HEAD
#' Data input can happen either through the parameter `x` (an object of type `rareData`)
#' or through the parameters `ai`,`bi`,`ci`, `di`, `n1i`, `n2i`, `data` (columns of a dataframe).
#' A `rareData` object can be produced from a data frame by applying the `rareDescribe()` function to it. The `rareDescribe()`
#' function pre-processes the data frame and stores the information required by the `rareMH()` function
||||||| 09c5356
#' The main input of the `rareMH()` function is a so-called `rareData` object. A `rareData` object
#' can be produced from a data frame by applying the `rareDescribe()` function to it. The `rareDescribe()`
#' function pre-processes the data frame and stores the information required by the `rareMH()` function
=======
#' The main input of the `rarePeto()` function is a so-called `rareData` object. A `rareData` object
#' can be produced from a data frame by applying the `rareDescribe()` function to it. The `rareDescribe()`
#' function pre-processes the data frame and stores the information required by the `rarePeto()` function
>>>>>>> 002bbf7b3861629eb0ed269627f9e31800d089a3
#' in a list. See `?rareDescribe` for more details.
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
#'
<<<<<<< HEAD
rarePeto <- function(x, ai, bi, ci, di, n1i, n2i, data, level = 95, digits = 4){
||||||| 09c5356
rarePeto <- function(x, level = 95, digits = 4){
=======
#' data <- data.frame(
#' ai = c(0, 3, 2, 0),
#' bi = c(20, 18, 15, 19),
#' ci = c(1, 4, 0, 0),
#' di = c(19, 17, 16, 20)
#' )
#'
#' x <- rareDescribe(ai = ai, bi = bi, ci = ci, di = di, data = data)
#'
#' rarePeto(x)
#'
rarePeto <- function(x, level = 95, digits = 4){
>>>>>>> 002bbf7b3861629eb0ed269627f9e31800d089a3


  ## argument checking ##

  # defining an object of class 'raredata' if raw data is put in
  if(missing(x)){
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

  Ai    <- ai + ci
  Bi    <- bi + di
  Vi    <- (Ai * Bi * n1i * n2i) / (ni^2 *(ni - 1))
  sumVi <- sum(Vi)

  quant <- stats::qnorm((100-level)/200)

  beta  <- sum((ci-Ai*(n1i/ni)))/sumVi
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
