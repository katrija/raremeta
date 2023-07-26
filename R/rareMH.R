#' Conduct a meta-analysis using Mantel-Haenszel type estimators
#'
#' Function to conduct a meta-analysis.
#' Effect sizes are the log odds ratio, log relative risk and
#' risk difference estimated from event counts in form of 2x2 contingency
#' tables using Mantel-Haenszel type estimators.
#'
#'
#' @param x an object of class `"rareData"`.
#' @param ai data frame column to specify the number of events in group 1 (i.e., the treatment group).
#' @param bi data frame column to specify the number of non-events in group 1 (i.e., the treatment group).
#' @param ci data frame column to specify the number of events in group 2 (i.e., the control group).
#' @param di data frame column to specify number of non-events in group 2 (i.e., the control group).
#' @param n1i data frame column to specify the sample sizes in group 1 (i.e., the treatment group).
#' @param n2i data frame column to specify the sample sizes in group 2 (i.e., the control group).
#' @param data data frame
#' @param measure character string specifying the effect size or outcome measure to be used
#' (either `"logOR"` for the log odds ratio, `"logRR"` for the log relative risk,
#' or `"RD"` for the risk difference).
#' @param level numeric between 0 and 100 specifying the confidence interval level (the default is 95).
#' @param digits integer specifying the number of decimal places to which the printed results
#' should be rounded (if unspecified, the default is 4).
#' @param correct logical specifying whether the data shall be continuity corrected before application of `rareMH()`. Default is `FALSE`
#' @param cc character string specifying the type of continuity correction to be used
#' (either `"none"`, `"constant"`, `"tacc"` or `"empirical"`). Default is `"constant"`.
#' @param ccval scalar or numerical vector specifying the value of the continuity correction if
#' `cc = "constant"`. Must be a scalar or a vector of length equal to the number of studies.
#' Default is `ccval = 0.5`. If a scalar is specified, the value is added to all studies for
#' which the number of events is zero in at least one of the groups. This behavior can be changed
#' using the argument `ccto`. `ccval` is overwritten by tccval and cccval if both arguments are specified.
#' @param tccval scalar or numerical vector specifying the value of the continuity correction
#' applied to the observations from the treatment group (group 1) if `cc = "constant"`. Must be a scalar or a vector
#' of length equal to the number of studies. If `cc = "constant"` and `tccval` is not specified, `tccval` is
#' set to the value of `ccval` internally.
#' @param cccval scalar or numerical vector specifying the value of the continuity correction
#' applied to the observations from the control group (group 2) if `cc = "constant"`. Must be a scalar or a vector
#' of length equal to the number of studies. If `cc = "constant"` and `cccval` is not specified, `cccval` is
#' set to the value of `ccval` internally.
#' @param ccsum numeric value specifying the value of the sum of the continuity correction applied to the
#' observations from the treatment group and the continuity correction applied to the observations from
#' the control group. Default is `ccsum = 1`. Currently, setting this argument to a different number only has
#' an effect when `cc = "tacc"` or `cc = "empirical"`.
#' @param ccto character string indicating to which studies the continuity correction should
#' be applied. Either `"only0"`, for which the continuity correction is applied to all studies
#' for which the number of events is zero in at least one of the groups, `"all"`, for which the
#' continuity correction is applied to all studies, or `"if0all"`, for which the continuity
#' correction is applied to all studies if any of the individual studies has zero events in at
#' least one of the groups.
#' @param method character string specifying whether a fixed- or a random-effects model should be fitted when
#' calculating the empirical and treatment-arm continuity correction. See ´?rareCC()´ for more details.
#' A fixed-effects model is fitted when using `method = "FE"` . A random-effects model is fitted
#' by setting `method` equal to one of the following: `"DL"`, `"HE"`, `"SJ"`, `"ML"`, `"REML"`, `"EB"`, `"HS"`,
#' `"PM"`, `"IPM"`, `"GENQ"`, `"PMM"` or `"GENQM"`.
#' @param drop00 logical indicating whether double-zero studies (i.e., studies with no events or only events in both groups)
#' should be excluded prior to calculating the studies' effect sizes and sampling variances.

#'
#' @details
#' # Details
#' ## Data input
#' The data input can be specified either through the arguments `ai`,`bi`,`ci`,`di`,`n1i`, and `n2i` (columns of the data frame `data`)
#' or through the argument `x`, which takes an object that results from applying the `rareDescribe()` function to the data
#' (i.e., the input for argument `x` must be an object of type `rareData`).
#' A `rareData` object can be produced from a data frame by applying the `rareDescribe()` function to it.
#' The `rareDescribe()` function pre-processes the data frame and stores the information required by the `rareMH()` function in a list.
#' See `?rareDescribe` for more details.
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
#' A prominent class of estimators for the common effect size under the the equal-effects model
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
#' @examples
#' # load a data set
#' data(dat.nissen2007)
#' d <- dat.nissen2007
#'
#' # meta analysis of the log odds ratio using the Mantel-Haenszel method
#' rareMH(ai=miRosiglitazone, ci=miControl, n1i=nRosiglitazone, n2i=nControl, data=d, measure="logOR")
#'
#' # same analysis with pre-processed data
#' x <- rareDescribe(ai=miRosiglitazone, ci=miControl, n1i=nRosiglitazone, n2i=nControl, data=d)
#' rareMH(x, measure="logOR")
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
                   measure, level = 95,  digits = 4,
                   correct = FALSE, cc = "constant", ccval = 0.5, tccval, cccval,
                   ccsum = 1, ccto = "only0", method = "FE"){

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

  # continuity correction if wanted
  # extract counts and sample sizes
  if(!is.logical(correct)) stop("argument 'correct' must be a logical")

  if(correct == TRUE){

    x <- rareCC(x, cc = cc, ccval = ccval, tccval = tccval, cccval = cccval,
                ccsum = ccsum, ccto = ccto, measure = measure, method = method)
  }

  if(!is.null(x$cc)){

    ai  <- x$ai.cc
    bi  <- x$bi.cc
    ci  <- x$ci.cc
    di  <- x$di.cc
    n1i <- x$n1i.cc
    n2i <- x$n2i.cc
    ni  <- n1i + n2i

  }else{
  ai  <- x$ai
  bi  <- x$bi
  ci  <- x$ci
  di  <- x$di
  n1i <- x$n1i
  n2i <- x$n2i
  ni  <- n1i + n2i
  }



  quant <- stats::qnorm((100-level)/200)

  #calculating estimated effect sizes
  if(measure == "logOR"){
    Ai <- (ai+di)/ni
    Bi <- (ai*di)/ni
    Ci <- (bi+ci)/ni
    Di <- (bi*ci)/ni

    B <- sum(Bi, na.rm = TRUE)
    D <- sum(Di, na.rm = TRUE)


    if(B == 0 || D == 0){
      stop("The data does not allow for the application of this method.\n
           This might be the result of no occurence of one of the outcomes in all of the studies.\n
           See e.g. ?rareCC() for possible continuity corrections.")
    }


    beta   <- log(B/D)
    se     <- sqrt(1/2 * (sum(Ai*Bi, na.rm = TRUE)/B^2 + sum(Ai*Di + Ci*Bi, na.rm = TRUE)/(B*D)
                    + sum(Ci*Di, na.rm = TRUE)/D^2))
    zval   <- beta / se
    pval   <- 2*stats::pnorm(abs(zval), lower.tail=FALSE)
    ci.lb  <- beta + quant * se
    ci.ub  <- beta - quant * se

  }

  if(measure == "logRR"){

    A <- sum(ai*n2i/ni, na.rm = TRUE)
    B <- sum(ci*n1i/ni, na.rm = TRUE)

    if(A == 0 || B == 0){
      stop("The data does not allow for application of this method. \n
           This might be the result of no events in one of the groups in all of the studies.\n
           See e.g. ?rareCC() for possible continuity corrections.")
    }


      beta   <- log(A/B)
      se     <- sqrt(sum(n1i*n2i*(ai+ci)/ni^2 - ai*ci/ni, na.rm = TRUE)/(A*B))
      zval   <- beta / se
      pval   <- 2*stats::pnorm(abs(zval), lower.tail=FALSE)
      ci.lb  <- beta + quant * se
      ci.ub  <- beta - quant * se

  }

  if(measure == "RD"){

    beta   <- sum(ai*(n2i/ni) - ci*(n1i/ni), na.rm = TRUE)/sum(n1i*(n2i/ni), na.rm = TRUE)
    se     <- sqrt((beta * (sum(ci*(n1i/ni)^2 - ai*(n2i/ni)^2 + (n1i/ni)*(n2i/ni)*(n2i-n1i)/2, na.rm = TRUE))
     + sum(ai*(n2i-ci)/ni + ci*(n1i-ai)/ni, na.rm = TRUE)/2) / sum(n1i*(n2i/ni), na.rm = TRUE)^2)

    #se     <- sqrt(sum((ai*bi*n2i^3 + ci*di*n1i^3)/(n1i*n2i*ni^2), na.rm = TRUE)/
    #                   (sum(n1i*n2i/ni, na.rm = TRUE)^2))
    # alternative estimator from Sato (1989)
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







