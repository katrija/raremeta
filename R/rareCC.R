#' Apply continuity correction to your data
#'
#' Function to apply different kinds of continuity corrections to meta-analytic data.
#'
#'
#' @param x an object of class `"rareData"`.
#' @param ai data frame column to specify the number of events in group 1 (i.e., the treatment group).
#' @param bi data frame column to specify the number of non-events in group 1 (i.e., the treatment group).
#' @param ci data frame column to specify the number of events in group 2 (i.e., the control group).
#' @param di data frame column to specify number of non-events in group 2 (i.e., the control group).
#' @param n1i data frame column to specify the sample sizes in group 1 (i.e., the treatment group).
#' @param n2i data frame column to specify the sample sizes in group 2 (i.e., the control group).
#' @param data data frame.
#' @param cc character string specifying the type of continuity correction to be used
#' (either `"constant"`, `"tacc"` or `"empirical"`). Default is "constant". See 'Details'.
#' @param ccval scalar or numerical vector specifying the value of the continuity correction if
#' `cc = "constant"`. Must be a scalar or a vector of length equal to the number of studies.
#' Default is `ccval = 0.5`. If a scalar is specified, the value is added to all studies for
#' which the number of events is zero in at least one of the groups. This behavior can be changed
#' using the argument `ccto`. The argument `ccval` is overwritten by `tccval` and `cccval` if both arguments are specified.
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
#' @param ccto character string indicating which studies the continuity correction should
#' be applied to. Either `"only0"`, for which the continuity correction is applied to all studies
#' for which the number of events is zero in at least one of the groups, `"all"`, for which the
#' continuity correction is applied to all studies, or `"if0all"`, for which the continuity
#' correction is applied to all studies if any of the individual studies has zero events in at
#' least one of the groups.
#' @param drop00 logical indicating whether double-zero studies (i.e., studies with no events or only events in both groups)
#' should be excluded prior to calculating the studies' effect sizes and sampling variances.
#' @param measure character string specifying the effect size or outcome measure to be used
#' (either `"logOR"` for the log odds ratio, `"logRR"` for the log relative risk,
#' or `"RD"` for the risk difference). Important when selecting continuity corrections dependent
#' on the estimated summary effect size of some subset of the studies involved in the analyisis,
#' i.e., `cc = "empirical"`.
#' @param method character string specifying whether a fixed- or a random-effects model should be fitted.
#' A fixed-effects model is fitted when using `method = "FE"` . A random-effects model is fitted
#' by setting `method` equal to one of the following: `"DL"`, `"HE"`, `"SJ"`, `"ML"`, `"REML"`, `"EB"`, `"HS"`,
#' `"PM"`, `"IPM"`, `"GENQ"`, `"PMM"` or `"GENQM"`. Important when selecting continuity corrections dependent
#' on the estimated summary effect size of some subset of the studies involved in the analyisis,
#' i.e., `cc = "empirical"`.
#'
#' @details
#' # Details
#' ## Data input
#' The data input can be specified either through the arguments `ai`,`bi`,`ci`,`di`,`n1i`, and `n2i` (columns of the data frame `data`)
#' or through the argument `x`, which takes an object that results from applying the `rareDescribe()` function to the data
#' (i.e., the input for argument `x` must be an object of type `rareData`).
#' A `rareData` object can be produced from a data frame by applying the `rareDescribe()` function to it.
#' The `rareDescribe()` function pre-processes the data frame and stores the information required by the `rareCC()` function in a list.
#' See `?rareDescribe` for more details.
#'
#' ## Types of continuity correction
#' This function offers three kinds of continuity correction.
#'
#' ### Constant continuity correction
#' When setting `cc = "constant"`, a constant value will be added to all
#' cells of all studies specified via the `ccto` argument. The default is `ccval = 0.5`.
#' Through specification of `tccval` and `cccval` different values can be specified for the
#' cells of the treatment group and the control group, respectively.
#' By inputting vectors in the afforementioned arguments, different values can be added
#' to different studies (adding the `k`-th value to the `k`-th study),
#' making it possible for the user to introduce their own continuity correction.
#'
#' ### Treatment-arm continuity correction
#' When setting `cc = "tacc"` (treatment-arm continuity correction), the size of
#' both of study arms is used to calculate a value for continuity correction.
#' For a precise definition see Sweeting et al. (2004).
#'
#' ### Empirical correction
#' When setting `cc = "empirical"`, the estimated summary effect size
#' (logOR, logRR or RD) for all studies which enable the estimation of the corresponding individual effect
#' size is used to calculate a value for continuity correction.
#' This means that there must be at least one study enabling the estimation of
#' the corresponding individual effect size and the `measure`- and `method` argument must be
#' specified.
#' For a precise definition see Sweeting et al. (2004).
#' When it comes to model fitting, there is the possibility to fit fixed-effects models (also known as equal-effects models)
#' and random-effects models using the inverse variance approach.
#' A fixed-effects model is fitted when `method` is set to `"FE"` (or `"EE"`). A random-effects model
#' is fitted when `method` is set to either `"DL"`, `"HE"`, `"SJ"`, `"ML"`, `"REML"`, `"EB"`, `"HS"`,
#' `"PM"`, `"IPM"`, `"GENQ"`, `"PMM"` or `"GENQM"`.
#' Currently, the model is fitted by applying the `rma` function from the `metafor` package, see Viechtbauer (2010).
#'
#' ## Estimation of the between-study variance in random-effects meta-analysis
#' Different estimators can be used to estimate the between-study variance, tau^2, in random-effects
#' meta-analysis. The estimator to be used is specified via the `methods` argument.
#' * `"DL"`: DerSimonian-Laird estimator
#' * `"HE"`: Hedges estimator
#' * `"SJ"`: Sidik-Jonkman estimator
#' * `"ML"`: maximum likelihood estimator
#' * `"REML"`: restricted maximum likelihood estimator
#' * `"EB"`: empirical Bayes estimator
#' * `"HS"`: Hunter-Schmidt estimator
#' * `"PM"`: Paule-Mandel estimator
#' * `"IPM"`: improved Paule-Mandel estimator
#' * `"GENQ"`: generalized Q-statistic estimator
#' * `"PMM"`: median-unbiased Paule-Mandel estimator
#' * `"GENQM"`: median-unbiased generalized Q-statistic estimator
#'
#' Most of these estimators are described in Zhang et al. (2021). For details on the improved Paule-Mandel estimator,
#' see also Bhaumik et al. (2012). The median-unbiased Paule-Mandel estimator and the median-unbiased generalized
#' Q-statistic estimator are described in Viechtbauer (2021).
#'
#' @return an object of class "rareData". The object is a list containing the following elements:
#' * `ai`, `bi`, `ci`, `di`: original entries of the 2x2 tables for all studies.
#' * `ai.cc`, `bi.cc`, `ci.cc`, `di.cc`: entries of the 2x2 tables for all studies after
#' application of the specified continuity correction.
#' * `ni`, `n1i`, `n2i`: original total and group sample sizes.
#' * `ni.cc`, `n1i.cc`, `n2i.cc`: total and group sample sizes after application of the
#' specified continuity correction.
#' * `tcc`, `ccc`: value of the specified continuity correction for the treatment group (group 1) and control
#' group (group 2).
#' * `cc.studies`: vector which indicates whether the continuity correction was applied
#' to a study.
#' * `remove`: logical vector indicating which studies were removed before
#' aplication of the continuity correction
#' * `kdz`,`ksz`: number of double-zero and single-zero studies.
#' * `k1sz`, `k2sz`: number of single-zero studies where the zero is in group 1 or group 2.
#' * `cc` : the type of continuity correction applied
#' * `drop00`: logical indicating whether double-zero studies were omitted
#' * ...
#'
#'
#'
#'
#' @export
#'
#' @examples
#' # introducing a data set
#' data("dat.nissen2007")
#' d <- dat.nissen2007
#'
#' # applying the default continuity correction (i.e. cc = "constant" with ccval = 0.5)
#' d.constant <- rareCC(ai=miRosiglitazone, n1i=nRosiglitazone, ci=miControl, n2i=nControl, data= d)
#' d.constant
#'
#' # applying the treatment-arm continuity correction with the correction in the treatment group and the control group summing up to 0.1
#' d.tacc <- rareCC(ai=miRosiglitazone, n1i=nRosiglitazone, ci=miControl, n2i=nControl, data= d, cc = "tacc", ccsum = 0.1, measure = "logOR")
#' d.tacc
#'
#' #applying the empirical continuity correction for the log odds ratio for the fixed effects model
#' d.emp <- rareCC(ai=miRosiglitazone, n1i=nRosiglitazone, ci=miControl, n2i=nControl, data= d, cc = "empirical", measure = "logOR", method = "FE")
#' d.emp
#'
#' # same analysis with pre-processed data
#' x <- rareDescribe(ai=miRosiglitazone, n1i=nRosiglitazone, ci=miControl, n2i=nControl, data= d)
#' x.emp <- rareCC(x, cc="empirical", measure="logOR", method="FE")
#'
#' @references
#' Borenstein, M., Hedges, L. V., Higgins, J. P., & Rothstein, H. R. (2021). Introduction to meta-analysis.
#' John Wiley & Sons.
#'
#' Bhaumik, D. K., Amatya, A., Normand, S.-L. T., Greenhouse, J., Kaizar, E., Neelon, B., & Gibbons, R. D. (2012).
#' Meta-analysis of rare binary adverse event data. Journal of the American Statistical
#' Association 107, 555–567. doi: 10.1080/01621459.2012.664484
#'
#' Sweeting, M. J., Sutton, A. J., & Lambert, P. C. (2004). What to add to nothing? Use and avoidance of
#' continuity corrections in meta-analysis of sparse data. Statistics in Medicine, 23, 1351–1375. doi: 10.1002/sim.1761
#'
#' Viechtbauer, W. (2010). Conducting Meta-Analyses in R with the metafor Package.
#' Journal of Statistical Software, 36(3), 1–48. https://doi.org/10.18637/jss.v036.i03
#'
#' Viechtbauer, W. (2021). Median-unbiased estimators for the amount of heterogeneity in meta-analysis. European Congress of Methodology,
#' Valencia, Spain. https://www.wvbauer.com/lib/exe/fetch.php/talks:2021_viechtbauer_eam_median_tau2.pdf
#'
#' Zhang, C., Chen, M., & Wang, X. (2020). Statistical methods for quantifying between-study heterogeneity in meta-analysis
#' with focus on rare binary events. Statistics and its interface, 13(4), 449. doi: 10.4310/sii.2020.v13.n4.a3


rareCC <- function(x, ai, bi, ci, di, n1i, n2i, data,
                   cc = "constant", ccval = 0.5, tccval, cccval, ccsum = 1,
                   ccto = "only0", drop00 = TRUE, measure, method = "FE"){

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

  # check if x was already continuity corrected
  if(!is.null(x$cc)){
    stop("This dataset was already continuity corrected before.")
  }

  # extract counts and sample sizes
  ai <- x$ai
  bi <- x$bi
  ci <- x$ci
  di <- x$di
  n1i <- x$n1i
  n2i <- x$n2i


  # check if cc argument is valid
  if(!is.element(cc, c("none", "constant", "tacc", "empirical"))){
    stop("'cc' must be either 'none', 'constant', 'tacc', or 'empirical'.")
  }

  # check if drop00 argument is valid
  if(!is.logical(drop00)){
    stop("'drop00' must be a logical.")
  }

  # check if ccto argument is valid (if cc shall be applied)
  if(cc != "none" &&  !is.element(ccto, c("only0", "all", "if0all", "none"))){
    stop("'ccto' must be either 'only0', 'all', 'none', or 'if0all'.")
  }

  # check ccval, tccval, cccval arguments (if cc == "constant")
  if(cc == "constant"){

    if(drop00 == FALSE && !is.element(length(ccval), c(1,x$k))){
      stop("'ccval' must have length 1 or length equal to the number of studies.")
    }

    if(drop00 == TRUE && !is.element(length(ccval), c(1,x$k,x$k-x$kdz))){
      stop("'ccval' must have length 1 or length equal to the number of studies (in- or excluding double-zero studies).")
    }

    if(any(ccval < 0)){
      stop("'ccval' must be non-negative.")
    }

    if(!missing(tccval) || !missing(cccval)){

      if((!missing(tccval) && missing(cccval)) || (missing(tccval) && !missing(cccval))){
        stop("Please specify both 'tccval' and 'cccval'.")
      }

      if(!equalLength(tccval, cccval)){
        stop("'tccval' and 'cccval' must have equal length.")
      }

      if(drop00 == FALSE && !is.element(length(tccval), c(1,x$k))){
        stop("'tccval' must have length 1 or length equal to the number of studies.")
      }

      if(drop00 == TRUE && !is.element(length(tccval), c(1,x$k,x$k-x$kdz))){
        stop("'tccval' must have length 1 or length equal to the number of studies (in- or excluding double-zero studies).")
      }

      if(any(c(tccval, cccval) < 0)){
        stop("All values in 'tccval' and 'cccval' must be non-negative.")
      }
    }
  }

  # check if measure argument is specified (if needed)
  if(cc == "empirical" && missing(measure)){
    stop("To apply the the empirical continuity correction, the 'measure' argument must be specified.")
  }

  # check if measure argument is valid (if needed)
  if(cc == "empirical" && !is.element(measure, c("logOR", "logRR", "RD"))){
    stop("To apply the empirical continuity correction, 'measure' must be either 'logOR', 'logRR' or 'RD'.")
  }

  # check if method argument is valid (if needed)
  if(cc == "empirical"
     && !is.element(method, c("FE", "EE", "DL", "HE", "HS", "HSk", "SJ", "ML",
                               "REML", "EB", "PM", "GENQ", "PMM", "GENQM"))){
    stop("'method' argument must be one of the following: 'FE','EE','DL','HE','HS',
         'HSK','SJ','ML','REML','EB','PM','GENQ','PMM','GENQM'.")
  }

  # tacc currently not supported for RD
  if(cc == "tacc" && measure == "RD"){
    stop("continuity correction of type 'tacc' is currently not supported for measure 'RD'.")
  }

  # empirical cc currently not supported for RD
  if(cc == "empirical" && measure == "RD"){
    stop("continuity correction of type 'empirical' is currently not supported for measure 'RD'.")
  }

  # check if there are non-zero studies when the empirical cc shall be applied
  if(cc == "empirical" && all(ai == 0 | bi == 0 | ci == 0 | di == 0)){
    stop("Continuity correction of type 'empirical' can only be applied if there is at least one non-zero study.")
  }

  # check ccsum argument:
  if(!is.numeric(ccsum) || length(ccsum) > 1 || ccsum < 0){
    stop("ccsum must be a scalar larger than 0.")
  }

  # check whether tccval and cccval are specified; if not, set them to ccval:
  if(missing(tccval)){
    tccval <- ccval
    cccval <- ccval
  }

  ## remove studies if desired ##

  remove <- rep(FALSE, length(ai))


  if(drop00 == TRUE) remove <- x$dzstudies

  if(cc == "none") remove <- (x$szstudies | x$dzstudies)

  if(cc == "none" && (sum(remove, na.rm = TRUE) > 0)){
    warning("The data contains studies with no event in at least one of the cells.
    Since cc = 'none' (no continuity correction) was selected, these studies were removed from the analysis.")
  }

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


  # specify studies to be continuity corrected

  ccstudies <- rep(FALSE, length(ai))

  if(cc == "none"){
     ccto  <- "none"
     ccval <- 0
  }

  if(ccto == "only0"){
    ccstudies <- (ai == 0 | ci == 0 | bi == 0 | di == 0)
  }

  if(ccto == "all" || (ccto == "if0all" && any(ai == 0 | ci == 0 | bi == 0 | di == 0) )){
    ccstudies <- rep(TRUE, length(ai))
  }

  # calculate the continuity correction

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
  if(cc == "tacc" && is.element(measure, c("logOR", "logRR"))){
    rinv <- n2i[ccstudies]/n1i[ccstudies]

    tcc[ccstudies] <- ccsum*1/(rinv+1)
    ccc[ccstudies] <- ccsum*rinv/(rinv+1)
  }

  # convert measure to the metafor notation
  if(cc == "empirical")
    metafor_measure <- sub("log", "", measure)

  if(cc == "empirical" && is.element(measure, c("logOR", "logRR"))){
    rinv <- n2i[ccstudies]/n1i[ccstudies]
    nozero <- which((ai != 0) & (bi != 0) & (ci != 0) & (di != 0))
    fit_nozero <- metafor::rma(ai = ai[nozero], bi = bi[nozero],
                               ci = ci[nozero], di = di[nozero],
                               measure = metafor_measure, method = method)
    prior <- exp(as.numeric(fit_nozero$beta))

    tcc[ccstudies] <- ccsum*prior/(rinv+prior)
    ccc[ccstudies] <- ccsum*rinv/(rinv+prior)

  }



  ai.cc <- ai; bi.cc <- bi; ci.cc <- ci; di.cc <- di

  # apply continuity correction
  ai.cc <- ai+tcc
  bi.cc <- bi+tcc
  ci.cc <- ci+ccc
  di.cc <- di+ccc

  if(any(ai.cc == 0 | bi.cc == 0 | ci.cc == 0 |di.cc == 0)){
    stop("There can not be studies with no event in either of the cells after the application of the continuity correction. \n
      See '?rareCC()' for more options.")
  }

  ## adding description of continuity corrected data ##

  data.cc <- data.frame(ai.cc = ai.cc, bi.cc = bi.cc, ci.cc = ci.cc, di.cc = di.cc)
  x.cc <- rareDescribe(ai.cc, bi.cc, ci.cc, di.cc, data = data.cc)

  out <- append(list(ai.cc = x.cc$ai, bi.cc = x.cc$bi, ci.cc = x.cc$ci,
                     di.cc = x.cc$di, n1i.cc = x.cc$n1i, n2i.cc = x.cc$n2i,
                     ni.cc = x.cc$ni, nratioi.cc = x.cc$nratioi,
                     nratio.cc = x.cc$nratio, k.cc = x.cc$k, kdz.cc = x.cc$kdz,
                     ksz.cc = x.cc$ksz, k1sz.cc = x.cc$k1sz, k2sz.cc = x.cc$k2sz,
                     n1.cc = x.cc$n1, n2.cc = x.cc$n2, n.cc = x.cc$n,
                     rf1.cc = x.cc$rf1, rf2.cc = x.cc$rf2, rf.cc = x.cc$rf,
                     krare.cc = x.cc$krare, kveryrare.cc = x.cc$kveryrare,
                     k1rare.cc = x.cc$k1rare, k1veryrare.cc = x.cc$k1veryrare,
                     k2rare.cc = x.cc$k2rare, k2veryrare.cc = x.cc$k2veryrare,
                     cc = cc, ccto = ccto, drop00 = drop00, ccstudies = ccstudies,
                     ccc = ccc, tcc = tcc, ccval = ccval, remove = remove), x)

  #report measure and method argument used in 'empirical' continuity correction
  if(cc == "empirical"){
    out <- append(out, list(method.cc = method, measure.cc = measure))
  }

  out <- rareData(out)
  return(out)

}

