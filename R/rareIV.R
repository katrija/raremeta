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
#' by setting `method` equal to one of the following: "DL", "HE", "SJ", "ML", "REML", "EB", "HS", "PM", "IPM", or "GENQ".
#' Default is `"SJ"`.
#' @param cc character string specifying the type of continuity corrections to be used
#' (either `"constant"`, `"tacc"` or `"empirical"`). Default is "constant". See 'Details'.
#' @param ccval scalar or numerical vector specifying the value of the continuity correction if
#' `cc = "constant"`. Must be a scalar or a vector of length equal to the number of studies.
#' Default is `ccval = 0.5`. If a scalar is specified, the value is added to all studies for
#' which the number of events is zero in at least one of the groups. This behavior can be changed
#' using the argument `ccto`. `ccval` is overwritten by tccval and cccval if both arguments are specified.
#' @param tccval scalar or numerical vector specifying the value of the continuity correction
#' applied to the observations from the treatment group if `cc = "constant"`. Must be a scalar or a vector
#' of length equal to the number of studies. If `cc = "constant"` and `tccval` is not specified, `tccval` is
#' set to the value of `ccval` internally.
#' @param cccval scalar or numerical vector specifying the value of the continuity correction
#' applied to the observations from the control group if `cc = "constant"`. Must be a scalar or a vector
#' of length equal to the number of studies. If `cc = "constant"` and `cccval` is not specified, `cccval` is
#' set to the value of `ccval` internally.
#' @param ccto character string indicating to which studies the continuity correction should
#' be applied. Either `"only0"`, for which the continuity correction is applied to all studies
#' for which the number of events is zero in at least one of the groups, `"all"`, for which the
#' continuity correction is applied to all studies, or `"if0all"`, for which the continuity
#' correction is applied to all studies if any of the individual studies has zero events in at
#' least one of the groups.
#' @param drop00 logical indicating whether double-zero studies (i.e., studies with no events or
#' only events in both groups) should be excluded when calculating the outcome measufit for the
#' individual studies.
#' @param test character string specifying how test statistics and confidence intervals for the
#' fixed effects should be computed (either `"z"`, for Wald-type tests and CIs, or `"hksj"`, for
#' tests and CIs based on the method by Knapp and Hartung (2003) and Sidik and Jonkman (2002).
#' Specifying `test = "knha"` instead of `test = "hksj"` will produce the same results).
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
rareIV <- function(x, measure, method, cc, ccval = 0.5, tccval, cccval, ccto = "only0",
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
  if(!is.element(cc, c("none", "constant", "tacc", "empirical"))){
    stop("'cc' must be either 'none', 'constant', 'tacc', or 'empirical'")
  }


  # if cc argument exists, check that ccto argument is valid
  if(!missing(cc)){
    if(!is.element(ccto, c("only0", "all", "if0all"))){
      stop("'ccto' must be either 'only0', 'all', or 'if0all'")
    }
  }

  # check if ccval argument is valid
  if(cc == "constant"){

    if(drop00 == FALSE & !is.element(length(ccval), c(1,x$k))){
      stop("'ccval' must have length 1 or length equal to the number of studies.")
    }

    if(drop00 == TRUE & !is.element(length(ccval), c(1,x$k,x$k-x$kdz))){
      stop("'ccval' must have length 1 or length equal to the number of studies (in- or excluding double-zero studies).")
    }

    if(!missing(tccval)){

      if((!missing(tccval) & missing(cccval))|(!missing(cccval) & missing(tccval))){
        stop("Please specify both 'tccval' and 'cccval'.")
      }

      if(!equalLength(tccval, cccval)){
        stop("'tccval' and 'cccval' must have equal length.")
      }

      if(drop00 == FALSE & !is.element(length(tccval), c(1,x$k))){
        stop("'ccval' must have length 1 or length equal to the number of studies.")
      }

      if(drop00 == TRUE & !is.element(length(tccval), c(1,x$k,x$k-x$kdz))){
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

  # check whether tccval and cccval are specified; if not, set them to ccval:
  if(missing(tccval)){
    tccval <- ccval
    cccval <- ccval
  }

  mf <- match.call()

  remove <- rep(FALSE, length(ai))

  # remove double-zero studies if desired:
  if(drop00 == TRUE){
    remove <- (ai == 0 & ci == 0) | (bi == 0 & di == 0)
    ai <- ai[!remove]
    bi <- bi[!remove]
    ci <- ci[!remove]
    di <- di[!remove]
    n1i <- n1i[!remove]
    n2i <- n2i[!remove]

    if(length(tccval) > length(ai)){
      tccval <- tccval[!remove]
      cccval <- cccval[!remove]
    }
  }

  ccstudies <- rep(FALSE, length(ai))

  # specify studies to be continuity corrected:
  if(ccto == "only0"){
    ccstudies <- (ai == 0 | ci == 0 | bi == 0 | di == 0)
  }

  if(ccto == "all" | (ccto == "if0all" & any(ai == 0 | ci == 0 | bi == 0 | di == 0) )){
    ccstudies <- rep(TRUE, length(ai))
  }

  tcc <- rep(0,length(ai))
  ccc <- rep(0,length(ci))

  if(cc == "constant"){

    if(length(tccval) == 1){
      tcc[ccstudies] <- tccval
      ccc[ccstudies] <- cccval
    }else{
      tcc <- tccval
      ccc <- cccval
    }

  }

  # continuity corrections for logOR and logRR:

  if(cc == "tacc" & is.element(measure, c("logOR", "logRR"))){
    rinv <- n2i[ccstudies]/n1i[ccstudies]

    tcc[ccstudies] <- 1/(rinv+1)
    ccc[ccstudies] <- rinv/(rinv+1)
  }

  if(cc == "empirical" & is.element(measure, c("logOR", "logRR"))){
    rinv <- n2i[ccstudies]/n1i[ccstudies]
    nozero <- which((ai != 0) & (bi != 0) & (ci != 0) & (di != 0))
    fit_nozero <- metafor::rma(ai = ai[nozero], bi = bi[nozero],
                               ci = ci[nozero], di = di[nozero],
                               measure = metafor_measure, method = method)
    prior <- exp(as.numeric(fit_nozero$beta))

    tcc[ccstudies] <- prior/(rinv+prior)
    ccc[ccstudies] <- rinv/(rinv+prior)

  }

  ai_cc <- ai; bi_cc <- bi; ci_cc <- ci; di_cc <- di

  # apply continuity correction:
  ai_cc <- ai+tcc
  bi_cc <- bi+tcc
  ci_cc <- ci+ccc
  di_cc <- di+ccc

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
  fit <- metafor::rma(ai = ai_cc, bi = bi_cc,
                      ci = ci_cc, di = di_cc,
                      measure = metafor_measure, method = metafor_method,
                      test = test,
                      tau2 = tau2,
                      to = "none", # prevent application of further continuity corrections
                      digits = digits, verbose = verbose,
                      control = control)


  # make results list
  # UNDER CONSTRUCTION
  res <- list(
    # model information:
    b = fit$b,
    beta = fit$beta,
    se = fit$se,
    zval = fit$zval,
    pval = fit$pval,
    ci.lb = fit$ci.lb,
    ci.ub = fit$ci.ub,
    vb = fit$vb,
    tau2 = fit$tau2,
    se.tau2 = fit$se.tau2,
    I2 = fit$I2,
    H2 = fit$H2,
    R2 = fit$R2,
    vt = fit$vt,
    QE = fit$QE,
    QEp = fit$QEp,
    QM = fit$QM,
    QMdf = fit$QMdf,
    QMp = fit$QMp,
    fit.stats = fit$fit.stats,
    p = fit$p,
    # study numbers:
    k = fit$k,
    k.all = x$k,
    kdz = x$kdz,
    ksz = x$ksz,
    k1sz = x$k1sz,
    k2sz = x$k2sz,
    ids = 1:length(ai),
    incl.studies = !remove,
    # effect sizes and sampling variances:
    yi = fit$yi,
    vi = fit$vi,
    # model matrix:
    X = fit$X,
    # counts and sample sizes::
    ai = ai,
    bi = bi,
    ci = ci,
    di = di,
    ai.cc = ai_cc,
    bi.cc = bi_cc,
    ci.cc = ci_cc,
    di.cc = di_cc,
    ni = n1i+n2i,
    n1i = n1i,
    n2i = n2i,
    n1i.cc = n1i_cc,
    n2i.cc = n2i_cc,
    # continuity corrections:
    tcc = tcc,
    ccc = ccc,
    cc.studies = ccstudies,
    # arguments:
    measure = measure,
    method = method,
    cc = cc,
    ccval = ccval,
    tccval = tccval,
    cccval = cccval,
    ccto = ccto,
    drop00 = drop00,
    test = test,
    digits = digits,
    # control = control,
    # package version and call:
    version = utils::packageVersion("raremeta"),
    call = mf
  )

  # apply class "raremeta"
  # UNDER CONSTRUCTION

  return(res)
}
