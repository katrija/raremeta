#' Conduct a meta-analysis of a rare event using the inverse variance model
#'
#' Function to conduct a meta-analysis of a rare event using the fixed- or
#' random-effects model of the inverse variance approach.
#'
#' @param x an object of class `"rareData"`.
#' @param measure character string specifying the effect size or outcome measure to be used
#' (either `"logOR"` for the log odds ratio, `"logRR"` for the log relative risk,
#' or `"RD"` for the risk difference).
#' @param method character string specifying whether a fixed- or a random-effects model should be fitted.
#' A fixed-effects model is fitted when using `method = "FE"`. A random-effects model is fitted
#' by setting `method` equal to one of the following: "DL", "HE", "SJ", "ML", "REML", "EB", "HS", "PM", or "GENQ".
#' Default is `"SJ"`.
#' @param cc character string specifying the type of continuity corrections to be used
#' (either `"constant"`, `"reciprocal"` or `"empirical"`). Default is "constant". See 'Details'.
#' @param ccval scalar or numerical vector specifying the value of the continuity correction if
#' `cc = "constant"`. Must be a scalar or a vector of length equal to the number of studies.
#' Default is `ccval = 0.5`. If a scalar is specified, the value is added to all studies for
#' which the number of events is zero in at least one of the groups. This behavior can be changed
#' using the argument `ccto`.
#' @param ccto character string indicating to which studies the continuity correction should
#' be applied. Either `"only0"`, for which the continuity correction is applied to all studies
#' for which the number of events is zero in at least one of the groups, `"all"`, for which the
#' continuity correction is applied to all studies, or `"if0all"`, for which the continuity
#' correction is applied to all studies if any of the individual studies has zero events in at
#' least one of the groups.
#' @param drop00 logical indicating whether double-zero studies (i.e., studies with no events or
#' only events in both groups) should be excluded when calculating the outcome measures for the
#' individual studies.
#' @param test character string specifying how test statistics and confidence intervals for the
#' fixed effects should be computed (either `"z"`, for Wald-type tests and CIs, or `"hksj"`, for
#' tests and CIs based on the method by Knapp and Hartung (2003) and Sidik and Jonkman (2002)).
#' @param digits integer specifying the number of decimal places to which the printed results
#' should be rounded (if unspecified, the default is 4).
#' @param verbose logical indicating whether output should be generated on the progress of model
#' fitting (the default is `FALSE`). Can also be an integer. Values > 1 generate more verbose output.
#' @param control optional list of control values for the iterative algorithms. If unspecified, default
#' values are defined inside the functions.
#' @param ... additional arguments.
#'
#' @return an object of class raremeta. The object is a list containing the following elements:
#'
#' ...
#'
#' @details
#'
#' ...
#'
#' @references
#' Hartung, J., and Knapp, G. (2001). A refined method for the meta-analysis of controlled clinical trials
#' with binary outcome. Statistics in Medicine 20, 3875-3889. doi: 10.1002/sim.1009
#'
#' Sidik, K., and Jonkman, J. N. (2002). A simple confidence interval for meta-analysis.
#' Statistics in Medicine 21, 3153-3159. doi: 10.1002/sim.1262
#'
#' @export
#'
rareIV <- function(x, measure, method, cc, ccval = 0.5, ccto = "only0",
                   drop00 = TRUE, test="z", digits = 3, verbose=FALSE, control,
                   ...){

  # check if x is an object of class rareData
  if(!class(x) == "rareData"){
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

  # check if method argument is specified
  if(missing(method)){
    stop("'method' argument must be specified.")
  }

  # check if method argument is valid
  if(!is.element(method, c("FE","EE","CE","HS","HSk","HE","DL","DLIT","GENQ",
                           "GENQM","SJ","SJIT","PM","MP","PMM","ML","REML","EB",
                           "IPM"))){
    stop("unknown 'method' specified See ?rareIV for more details.")
  }

  # extract counts and sample sizes
  ai <- x$ai
  bi <- x$bi
  ci <- x$ci
  di <- x$di
  n1i <- x$n1i
  n2i <- x$n2i

  # check if cc is specified (in case there are zero-studies)
  if(any(c(ai,bi,ci,di) == 0) & missing(cc)){
    stop("Some studies have zero events. \n
         You must specify the 'cc' argument to determine how they are handled.\n
         In case you want to exclude all zero-studies, set 'cc' equal to 'none'.")
  }

  if(!is.logical(drop00)){
    stop("'drop00' must be a logical")
  }

  # check if cc argument is valid
  if(!is.element(cc, c("none", "constant", "reciprocal", "empirical"))){
    stop("'cc' must be either 'none', 'constant', 'reciprocal', or 'empirical'")
  }


  # if cc argument exists, check that ccto argument is valid
  if(!missing(cc)){
    if(!is.element(ccto, c("only0", "all", "if0all"))){
      stop("'ccto' must be either 'only0', 'all', or 'if0all'")
    }
  }

  # check if ccval argument is valid
  if(cc == "constant"){
    if(drop00 == FALSE){
      if(!is.element(length(ccval), c(1,x$k))){
        stop("'ccval' must have length 1 or length equal to the number of studies.")
      }
    }
    if(drop00 == TRUE){
      if(!is.element(length(ccval), c(1,x$k,x$k-x$kdz))){
        stop("'ccval' must have length 1 or length equal to the number of studies (in- or excluding double-zero studies).")
      }
    }
  }

  # check if test argument is valid
  if(!is.element(test, c("z", "knha", "hksj"))){
    stop("'test' must be either 'z', 'knha', or 'hksj'")
  }

  # check if digits argument is valid
  if(length(digits) != 1 | digits%%1 != 0 | digits < 0){
    stop("'digits' must be an integer of length 1.")
  }

  # check if ccval is non-negative
  if(cc == "constant" & any(ccval < 0)){
    stop("'ccval' must be non-negative.")
  }

  # empirical cc currently not supported for RD
  if(cc == "empirical" & measure == "RD"){
    stop("continuity correction of type 'empirical' is currently not supported for measure 'RD'.")
  }

  # check that there are non-zero studies when the empirical cc shall be applied
  if(cc == "empirical" & all(ai == 0 | bi == 0 | ci == 0 | di == 0)){
    stop("continuity correction of type 'empirical' can only be applied if there is at least one non-zero study.")
  }

  if(cc == "empirical" & method == "IPM"){
    stop("method = 'IPM' is currently not defined for cc = 'empirical'.")
  }

  # convert measure to the metafor notation
  metafor_measure <- sub("log", "", measure)

  # convert test to metafor notation if needed
  if(test == "hksj"){
    test = "knha"
  }

  # remove double-zero studies if desired:
  if(drop00 == TRUE){
    remove <- (ai == 0 & ci == 0) | (bi == 0 & di == 0)
    ai <- ai[!remove]
    bi <- bi[!remove]
    ci <- ci[!remove]
    di <- di[!remove]
    n1i <- n1i[!remove]
    n2i <- n2i[!remove]

    if(length(ccval) > length(ai)){
      ccval <- ccval[!remove]
    }
  }

  # specify studies to be continuity corrected:
  if(ccto == "only0"){
    ccstudies <- (ai == 0 | ci == 0 | bi == 0 | di == 0)
  }

  if(ccto == "all" | (ccto == "if0all" & any(ai == 0 | ci == 0 | bi == 0 | di == 0) )){
    ccstudies <- rep(TRUE, length(ai))
  }

  if(cc == "constant"){
    tcc <- ccval
    ccc <- ccval

    if(length(ccval) == 1){
      tcc <- rep(tcc, length(ai))
      ccc <- rep(ccc, length(ai))
    }
  }

  if(cc == "reciprocal"){
    rinv <- n2i/n1i

    tcc <- 1/(rinv+1)
    ccc <- rinv/(rinv+1)
  }

  if(cc == "empirical" & is.element(measure, c("logOR", "logRR"))){
    rinv <- n2i/n1i
    nozero <- which((ai != 0) & (bi != 0) & (ci != 0) & (di != 0))
    fit_nozero <- metafor::rma(ai = ai[nozero], bi = bi[nozero],
                               ci = ci[nozero], di = di[nozero],
                               measure = metafor_measure, method = method)
    prior <- exp(as.numeric(fit_nozero$beta))

    tcc <- prior/(rinv+prior)
    ccc <- rinv/(rinv+prior)

  }

  ai_cc <- ai; bi_cc <- bi; ci_cc <- ci; di_cc <- di

  # apply continuity correction:
  ai_cc[ccstudies] <- ai[ccstudies]+tcc[ccstudies]
  bi_cc[ccstudies] <- bi[ccstudies]+tcc[ccstudies]
  ci_cc[ccstudies] <- ci[ccstudies]+ccc[ccstudies]
  di_cc[ccstudies] <- di[ccstudies]+ccc[ccstudies]

  n1i_cc <- ai_cc+bi_cc
  n2i_cc <- ci_cc+di_cc

  # treat the case of method = "IPM" (improved Paule-Mandel)
  # calculate the IPM estimate and assign it to tau2
  # then, rma will handle tau2 as known and use the IPM estimate in all its calculations
  if(method == "IPM"){
    metafor_method = "REML" # workaround - method needs to be specified properly when rma is used later on

    if(measure != "logOR"){
      stop("method = 'IPM' is only defined for measure = 'logOR'.")
    }

    if(!(cc == "constant" & ccto == "all")){
      warning("method = 'IPM' was developed for cc = 'constant' and ccto = 'all'.
              Using method = 'IPM' with other values for cc and ccto has never been evaluated empirically.")
    }

    yi <- log(ai_cc*di_cc/(bi_cc*ci_cc))
    vi <- 1/ai_cc+1/bi_cc+1/ci_cc+1/di_cc

    mu_hat <- mean(log(ci_cc/(n2i_cc-ci_cc)))
    theta_hat <- mean(yi)

    ipm_optim <- function(tau2){
      k <- length(yi)

      sigma2i <- 1/n1i_cc*(exp(-mu_hat-theta_hat+tau2/2)+2+exp(mu_hat+theta_hat+tau2/2))+1/n2i_cc*(exp(-mu_hat)+2+exp(mu_hat))
      # note that here, it is n1i+1 and no n1i because the continuity correction is added here
      # so it is crucial to make sure that the cc has not yet been applied to n1i!

      weights_ipm <- 1/(sigma2i+tau2)

      theta_w <- sum(weights_ipm*yi)/sum(weights_ipm)

      f_tau2 <- sum(weights_ipm*(yi-theta_w)^2)-(k-1)

      return(f_tau2)
    }

    if(ipm_optim(0) < 0){
      tau2 <- 0
    }else{
      tau2 <- stats::uniroot(ipm_optim, c(0,100))$root
    }
  }else{
    # what to do if method != "IPM"
    metafor_method = method
    tau2 = NULL
  }

  # run metafor
  res <- metafor::rma(ai = ai_cc, bi = bi_cc,
                      ci = ci_cc, di = di_cc,
                      measure = metafor_measure, method = metafor_method,
                      test = test,
                      tau2 = tau2,
                      to = "none", # prevent application of further continuity corrections
                      digits = digits, verbose = verbose)

  # make results list
  # UNDER CONSTRUCTION

  # apply class "raremeta"
  # UNDER CONSTRUCTION

  return(res)
}
