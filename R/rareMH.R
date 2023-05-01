#' Conduct a meta-analysis using Mantel-Haenszel type estimators
#'
#' Function to conduct a meta-analysis.
#' Effect sizes are the log odds ratio, log relative risk and
#' risk difference estimated from event counts in form of 2x2 contingency
#' tables using Mantel-Haenszel type estimators.
#'
#'
<<<<<<< HEAD
#' @param x an object of class `"raredata"`
#' @param ai Data frame column to specify the number of events in group 1 (i.e., the treatment group).
#' @param bi Data frame column to specify the number of non-events in group 1 (i.e., the treatment group).
#' @param ci Data frame column to specify the number of events in group 2 (i.e., the control group).
#' @param di Data frame column to specify number of non-events in group 2 (i.e., the control group).
#' @param n1i Data frame column to specify the sample sizes in group 1 (i.e., the treatment group).
#' @param n2i Data frame column to specify the sample sizes in group 2 (i.e., the control group).
#' @param data Data frame.
#' @param measure measure character string specifying the effect size or outcome measure to be used
||||||| 09c5356
#' @param x an object of class `"raredata"`
#' @param measure measure character string specifying the effect size or outcome measure to be used
=======
#' @param x an object of class `"rareData"`
#' @param measure character string specifying the effect size or outcome measure to be used
>>>>>>> 002bbf7b3861629eb0ed269627f9e31800d089a3
#' (either `"logOR"` for the log odds ratio, `"logRR"` for the log relative risk,
#' or `"RD"` for the risk difference).
#' @param level numeric between 0 and 100 specifying the confidence interval level (the default is 95).
#' @param digits integer specifying the number of decimal places to which the printed results
#' should be rounded (if unspecified, the default is 4).
#'
#' @details
#' # Details
#' ## Data input
#' Data input can happen either through the parameter `x` (an object of type `rareData`)
#' or through the parameters `ai`,`bi`,`ci`, `di`, `n1i`, `n2i`, `data` (columns of a dataframe). A `rareData` object
#' can be produced from a data frame by applying the `rareDescribe()` function to it. The `rareDescribe()`
#' function pre-processes the data frame and stores the information required by the `rareMH()` function
#' in a list. See `?rareDescribe` for more details.
#'
#' ## Effect size measures
#' The function includes meta-analytic methods for log odds ratios, log risk ratios, and risk differences.
#' The effect size measure can be specified using the `measure` argument. The respective effect size,
#' along with an estimate of its sampling variance, is then calculated for each study based on the
#' entries of the study's 2x2 table:
#'
#' |                   | event | no event|
#' |:------------------|------:|--------:|
#' |group1 (treatment) | ai    | bi      |
#' |group2 (control)   | ci    | di      |
#'
#' ## Mantel-Haenszel type estimators
#' A prominent class of estimators for the common effect size under the the equal-efects model
#' was proposed by Mantel and Haenszel in 1959.
#' Mantel-Haenszel type estimators are particularly useful in meta-analysis of rare
#' events.
#' As a quotient of sums rather than a sum of quotients,
#' the Mantel-Haenszel estimator of the overall effect is oftentimes defined even though the standard
#' estimator for the corresponding effect size for some individual studies may not be defined.
#' The variance estimation is due to Greenland et al. (1985) and
#' Robins et al. (1986).
#' For a summary of the implemented estimators see e.g. sections 9.2-9.4 of Jewell
#' (2013).
#'
#' @return an object of class "raremeta".
#' The object is a list containing the following elements:
#' * `beta`, `b`: estimated effect size.
#' * `se`: estimated standard error of the  estimator.
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
#' Mantel, N., & Haenszel, W. (1959).
#' Statistical aspects of the analysis of data from retrospective studies of disease.
#' Journal of the national cancer institute, 22(4), 719-748.
#'
#' Greenland, S., & Robins, J. M. (1985).
#' Estimation of a common effect parameter from sparse follow-up data.
#' Biometrics, 55-68.
#'
#' Robins, J., Breslow, N., & Greenland, S. (1986).
#' Estimators of the Mantel-Haenszel variance consistent in both sparse data
#' and large-strata limiting models.
#' Biometrics, 311-323.
#'
#' Jewell, N. P. (2003). Statistics for epidemiology. Chapman and Hall/CRC.
#'
#' @export
#'

rareMH <- function(x, ai, bi, ci, di, n1i, n2i, data,
                   measure, level = 95,  digits = 4){

  ## argument checking ##

  # defining an object of class 'raredata' if raw data is put in
  if(missing(x)){
    x <- rareDescribe(ai=ai, bi=bi, ci=ci, di=di, n1i=n1i, n2i=n2i, data=data)
  }

  # check if x is an object of class rareData
  if(!inherits(x,"rareData")){
    stop("x must be an object of class 'rareData'. See ?rareDescribe for more details.")
  }

  # check if measure argument is specified
  if(missing(measure)){
    stop("'measure' argument must be specified.")
  }

  # check if measure argument is valid
  if(!is.element(measure, c("logOR", "logRR", "RD"))){
    stop("'measure' must be either 'logOR', 'logRR', or 'RD'.")
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

  quant <- stats::qnorm((100-level)/200)

  #calculating estimated effect sizes
  if(measure == "logOR"){
    Ai <- (ai+di)/ni
    Bi <- (ai*di)/ni
    Ci <- (bi+ci)/ni
    Di <- (bi*ci)/ni

    B <- sum(Bi)
    D <- sum(Di)

    if(B == 0 || D == 0){
      stop("The data does not allow for the application of this method.
           See e.g. ?rareCC() for possible continuity corrections.")
    }


    beta   <- log(B/D)
    se     <- sqrt(1/2 * (sum(Ai*Bi)/B^2 + sum(Ai*Di + Ci*Bi)/(B*D)
                    + sum(Ci*Di)/D^2))
    zval   <- beta / se
    pval   <- 2*stats::pnorm(abs(zval), lower.tail=FALSE)
    ci.lb  <- beta + quant * se
    ci.ub  <- beta - quant * se

  }

  if(measure == "logRR"){

    A <- sum(ai*n2i/ni)
    B <- sum(ci*n1i/ni)

    if(A == 0 || B == 0){
      stop("The data does not allow for application of this method.
           See e.g. ?rareCC() for possible continuity corrections.")
    }


      beta   <- log(A/B)
      se     <- sqrt(sum(n1i*n2i*(ai+ci)/ni^2 - ai*ci/ni)/(A*B))
      zval   <- beta / se
      pval   <- 2*stats::pnorm(abs(zval), lower.tail=FALSE)
      ci.lb  <- beta + quant * se
      ci.ub  <- beta - quant * se

  }

  if(measure == "RD"){

    beta   <- sum(ai*(n2i/ni) - ci*(n1i/ni))/sum(n1i*(n2i/ni))
    #se     <- sqrt(sum((ai*bi*n2i^3 + ci*di*n1i^3)/(n1i*n2i*ni^2))/
    #                   (sum(n1i*n2i/ni)^2))
    se     <- sqrt((beta * (sum(ci*(n1i/ni)^2 - ai*(n2i/ni)^2 + (n1i/ni)*(n2i/ni)*(n2i-n1i)/2))
     + sum(ai*(n2i-ci)/ni + ci*(n1i-ai)/ni)/2) / sum(n1i*(n2i/ni))^2)
    # alternative estimator from Sato (1989?)
    # equation from 'A new and improved confidence interval for the Mantel–Haenszel risk difference'
    # by B Klingenberg

    zval   <- beta / se
    pval   <- 2*stats::pnorm(abs(zval), lower.tail=FALSE)
    ci.lb  <- beta + quant * se
    ci.ub  <- beta - quant * se
  }

  res <- list(
    # model information
    model = "rareMH",
    b = beta,
    beta = beta,
    se = se,
    zval = zval,
    pval = pval,
    ci.lb = ci.lb,
    ci.ub = ci.ub,
    digits = digits,
    #b.exp = exp(beta),
    #beta.exp = exp(beta),
    #ci.lb.exp = exp(ci.lb),
    #ci.ub.exp = exp(ci.ub)
    k = x$k,
    kdz = x$kdz,
    ksz = x$ksz,
    k1sz = x$k1sz,
    k2sz = x$k2sz,
    ai = ai,
    bi = bi,
    ci = ci,
    di = di,
    ni = n1i+n2i,
    n1i = n1i,
    n2i = n2i,
    # arguments:
    measure = measure,
    level = level,
    digits = digits,
    # package version and call:
    version = utils::packageVersion("raremeta"),
    call = mf
  )

  res <- raremeta(res)
  return(res)
}

#Variance Estimators
#references to the original papers
#OR: Robins et. al 1986 (p. 312)
#RR and RD: Greenland et. al 1985 (p. 63)







