#' Conduct a meta-analysis of a rare event using the inverse variance model
#'
#' Function to conduct a meta-analysis of a rare event using the fixed- or
#' random-effects model of the inverse variance approach. See below for more details
#' on these models and their application in meta-analyses of rare events.
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
#' @param method character string specifying whether a fixed- or a random-effects model should be fitted.
#' A fixed-effects model is fitted when using `method = "FE"` . A random-effects model is fitted
#' by setting `method` equal to one of the following: `"DL"`, `"HE"`, `"SJ"`, `"ML"`, `"REML"`, `"EB"`, `"HS"`,
#' `"PM"`, `"IPM"`, `"GENQ"`, `"PMM"` or `"GENQM"`.
#' @param cc character string specifying the type of continuity corrections to be used
#' (either `"constant"`, `"tacc"` or `"empirical"`). Default is "constant". See 'Details'.
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
#' @param drop00 logical indicating whether double-zero studies (i.e., studies with no events or only events in both groups)
#' should be excluded prior to calculating the studies' effect sizes and sampling variances.
#' @param weighted logical specifying whether weighted (default) or unweighted estimation of the pooled
#' effect size should be used.
#' @param weights numeric specifying user-defined weights to be used when fitting the model.
#' @param level numeric between 0 and 100 specifying the confidence interval level (the default is 95)
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
#' @details
#' # Details
#' ## Data input
#' The data input can be specified either through the arguments `ai`,`bi`,`ci`,`di`,`n1i`, and `n2i` (columns of the data frame `data`)
#' or through the argument `x`, which takes an object that results from applying the `rareDescribe()` function to the data
#' (i.e., the input for argument `x` must be an object of type `rareData`).
#' A `rareData` object can be produced from a data frame by applying the `rareDescribe()` function to it.
#' The `rareDescribe()` function pre-processes the data frame and stores the information required by the `rareIV()` function in a list.
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
#' ## Handling single-zero and double-zero studies
#' Single-zero studies are studies for which one entry of the 2x2 table is zero. Double-zero studies are
#' studies for which two entries of the same column of the 2x2 table are zero.
#' The function includes a variety of arguments to handle single-zero and double-zero studies. Per default,
#' double-zero studies are currently excluded from the analysis (this behavior might be changed in the future).
#' The inclusion of double-zero studies can be enforced by setting the argument `drop00` to `FALSE`.
#' If the data includes at least one single-zero study, the function throws an error if the user
#' did not specify whether or not a continuity correction shall be applied.
#' By setting `cc = "none"`, any zero-studies (studies with at least one zero-cell)
#' which remain in the data are excluded from the analysis. If it is desired that zero-studies are
#' included, the user needs to specify which continuity correction shall be used, and to which
#' studies it shall be applied. Per default, the continuity correction is only applied to zero-studies, while
#' studies for which all cells are larger than zero are left uncorrected. This behavior can be changed
#' using the argument `ccto`. Per default, the constant value 0.5 (`cc = "constant"`, `ccval = 0.5`) is added to
#' all cells of the studies specified by `ccto`. This continuity correcton was desribed by Gart and Zweifel (1967).
#' Alternative continuity corrections which were described by Sweeting et al. (2004) can be applied by setting `cc` to `"tacc"`
#' (for the treatment-arm continuity correction), and to `"empirical"` for the empirical continuity correction.
#' Per default sum of the corrections for treatment and control groups is set to `1`, but this can be changed by setting the
#' the argument `ccsum` to a different value.
#' It is possible to set the continuity correction to a user-defined value (or a vector of user-defined values) using the
#' argument `ccval` (if the value). If the user wants to specify different values for the treatment and the control group,
#' this is possible via the arguments `tccval` and `cccval`.
#'
#' ### Differences between effect size measures in the application of continuity corrections
#' When either the log odds ratio or the log risk ratio is used as an effect size measure, both the effect sizes and
#' their sampling variances are calculated based on the continuity-corrected 2x2 table.
#' When the effect size measure is the risk difference, the continuity-corrected 2x2 table is only used in the
#' calculation of the sampling variances.
#'
#' ## Model specification
#' The function can be used to fit fixed-effects models (also known as equal-effects models)
#' and random-effects models using the inverse variance approach. Currently, it is not possible to include moderators in any of these models.
#' A fixed-effects model is fitted when `method` is set to `"FE"` (or `"EE"`). A random-effects model
#' is fitted when `method` is set to either `"DL"`, `"HE"`, `"SJ"`, `"ML"`, `"REML"`, `"EB"`, `"HS"`,
#' `"PM"`, `"IPM"`, `"GENQ"`, `"PMM"` or `"GENQM"`. See below for details on heterogeneity estimation in random-effects meta-analysis.
#' For a basic introduction to fixed-effects and random-effects meta-analysis, please refer to Borenstein et al. (2021)
#' In usual applications of fixed- and random-effects meta-analyses, weighted estimation is used, i.e.,
#' the pooled effect size is calculated as a weighted average of the study effect sizes, where the weights are
#' defined as the inverse variances of the effect sizes. It is possible to switch to unweighted estimation
#' via `weighted = FALSE`.
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
#' ## Tests and confidence intervals for model coefficients
#' Currently, tests and confidence intervals for model coefficients are obtained from the output of the
#' function `rma()` from the `metafor()` package. Note that setting the argument `test`
#' to `"hksj"` produces the same result as setting it to `"knha"`. See `?metafor::rma` for further details.
#'
#' ## Tests for residual heterogeneity
#' Currently, the results of the test for residual heterogeneity are obtained from the output of the function
#' `rma()` from the `metafor()` package. See `?metafor::rma` for further details.
#'
#'
#' @return an object of class "raremeta". The object is a list containing the following elements:
#' * `model`: name of the model used for conducting the meta-analysis.
#' * `beta`: estimated coefficients of the model.
#' * `se`: standard errors of the  coefficients.
#' * `zval`: test statistics of the coefficients.
#' * `pval`: p-values corresponding to the test statistics.
#' * `ci.lb`: lower bound of the confidence intervals for the coefficients.
#' * `ci.ub`: upper bound of the confidence intervals for the coefficients.
#' * `vb`: variance-covariance matrix of the estimated coefficients.
#' * `tau2`: estimated amount of (residual) heterogeneity. Always `0` when `method = "FE"`.
#' * `se.tau2`: standard error of the estimated amount of (residual) heterogeneity.
#' * `I2`: value of I^2 (total heterogeneity/total variability).
#' * `H2`: value of H^2 (total variability/sampling variability).
#' * `R2`: value of R^2.
#' * `QE`: test statistic of the test for (residual) heterogeneity.
#' * `QEp`: p-value corresponding to the test statistic.
#' * `fit.stats`: a list with log-likelihood, deviance, AIC, BIC, and AICc values under
#' the unrestricted and restricted likelihood.
#' * `p`: number of coefficients in the model (including the intercept).
#' * `k`: number of studies included in the analysis.
#' * `k.all`: total number of studies (before exclusion).
#' * `kdz`,`ksz`: number of double-zero and single-zero studies.
#' * `k1sz`, `k2sz`: number of single-zero studies where the zero is in group 1 or group 2.
#' * `yi`, `vi`: vectors containing the estimated effect sizes and their estimated sampling
#' variances for all study.
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
#' * ...
#'
#' @examples
#' # introducing a data set
#' data(dat.nissen2007)
#' d <- dat.nissen2007
#'
#' # estimating the log relative risk in the fixed-effects model
#' rareIV(ai=miRosiglitazone, ci=miControl, n1i=nRosiglitazone, n2i=nControl, data=d, measure="logOR", method="FE", cc="constant")
#'
#' # estimating the log relative risk in a random effects-model
#' rareIV(ai=miRosiglitazone, ci=miControl, n1i=nRosiglitazone, n2i=nControl, data=d, measure="logRR", method="DL", cc="constant")
#'
#' # same analysis with pre-processed data
#' x <- rareDescribe(ai=miRosiglitazone, ci=miControl, n1i=nRosiglitazone, n2i=nControl, data=d)
#' rareIV(x, measure="logRR", method="DL", cc="constant")
#'
#'
#' @references
#' Borenstein, M., Hedges, L. V., Higgins, J. P., & Rothstein, H. R. (2021). Introduction to meta-analysis.
#' John Wiley & Sons.
#'
#' Bhaumik, D. K., Amatya, A., Normand, S.-L. T., Greenhouse, J., Kaizar, E., Neelon, B., & Gibbons, R. D. (2012).
#' Meta-analysis of rare binary adverse event data. Journal of the American Statistical
#' Association 107, 555–567. doi: 10.1080/01621459.2012.664484
#'
#' Gart, John J, and James R Zweifel. 1967. On the bias of various estimators of the logit and
#' its variance with application to quantal bioassay. Biometrika, 54, 181–187. doi:10.1093/BIOMET/54.1-2.181
#'
#' Hartung, J., and Knapp, G. (2001). A refined method for the meta-analysis of controlled clinical trials
#' with binary outcome. Statistics in Medicine, 20, 3875-3889. doi: 10.1002/sim.1009
#'
#' Sidik, K., and Jonkman, J. N. (2002). A simple confidence interval for meta-analysis.
#' Statistics in Medicine, 21, 3153-3159. doi: 10.1002/sim.1262
#'
#' Sweeting, M. J., Sutton, A. J., & Lambert, P. C. (2004). What to add to nothing? Use and avoidance of
#' continuity corrections in meta-analysis of sparse data. Statistics in Medicine, 23, 1351–1375. doi: 10.1002/sim.1761
#'
#' Viechtbauer, W. (2021). Median-unbiased estimators for the amount of heterogeneity in meta-analysis. European Congress of Methodology,
#' Valencia, Spain. https://www.wvbauer.com/lib/exe/fetch.php/talks:2021_viechtbauer_eam_median_tau2.pdf
#'
#' Zhang, C., Chen, M., & Wang, X. (2020). Statistical methods for quantifying between-study heterogeneity in meta-analysis
#' with focus on rare binary events. Statistics and its interface, 13(4), 449. doi: 10.4310/sii.2020.v13.n4.a3
#'
#' @export
#'
rareIV <- function(x, ai, bi, ci, di, n1i, n2i, data,
                   measure, method, cc, ccval = 0.5, tccval, cccval, ccsum = 1,
                   ccto = "only0",
                   drop00 = TRUE, weighted = TRUE, weights,
                   level = 95,
                   test="z", digits = 4, verbose=FALSE, control,
                   ...){


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

  # check if cc is specified
  # (in case there are zero-studies and continuity correction was not applied beforehand)
  if(any(c(ai,bi,ci,di) == 0) & is.null(x$cc) & missing(cc)){
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

    # check if ccval is non-negative
    if(any(ccval < 0)){
      stop("'ccval' must be non-negative.")
    }

    if(!missing(tccval) | !missing(cccval)){

      if((!missing(tccval) & missing(cccval))|(!missing(cccval) & missing(tccval))){
        stop("Please specify both 'tccval' and 'cccval'.")
      }

      if(!equalLength(tccval, cccval)){
        stop("'tccval' and 'cccval' must have equal length.")
      }

      if(drop00 == FALSE & !is.element(length(tccval), c(1,x$k))){
        stop("'tccval' must have length 1 or length equal to the number of studies.")
      }

      if(drop00 == TRUE & !is.element(length(tccval), c(1,x$k,x$k-x$kdz))){
        stop("'tccval' must have length 1 or length equal to the number of studies (in- or excluding double-zero studies).")
      }

      if(any(c(tccval, cccval) < 0)){
        stop("All values in 'tccval' and 'cccval' must be non-negative.")
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

  # tacc currently not supported for RD
  if(cc == "tacc" & measure == "RD"){
    stop("continuity correction of type 'tacc' is currently not supported for measure 'RD'.")
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

  # check ccsum argument:
  if(!is.numeric(ccsum) | length(ccsum) > 1 | ccsum < 0){
    stop("ccsum must be a scalar larger than 0.")
  }

  # check if level argument is valid:
  if(!is.numeric(level) | length(level) > 1 | level < 0 | level > 100){
    stop("level must be a scalar between 0 and 100.")
  }

  # convert measure to the metafor notation
  metafor_measure <- sub("log", "", measure)

  # convert test to metafor notation if needed
  if(test == "hksj"){
    test = "knha"
  }


  if(missing(cc)){
    cc = "none"
  }

  # check whether tccval and cccval are specified; if not, set them to ccval:
  if(missing(tccval)){
    tccval <- ccval
    cccval <- ccval
  }

  udweights <- FALSE

  # check whether weights are defined
  if(!missing(weights)){
    udweights <- TRUE

    if(drop00 == FALSE & !is.element(length(weights), c(1, x$k))|drop00 == TRUE & !is.element(length(weights), c(1, x$k, x$k-x$kdz))){
      stop("'weights' must be of length 1 or of length equal to the number of studies.")
    }
  }

  # calculate effect sizes with the specified continuity correction:
  es <- rareES(x, measure = measure, cc = cc, ccval = ccval, tccval = tccval, cccval = cccval,
               ccsum = ccsum, ccto = ccto, drop00 = drop00)

  yi <- es$yi
  vi <- es$vi
  ai.cc <- es$ai.cc
  bi.cc <- es$bi.cc
  ci.cc <- es$ci.cc
  di.cc <- es$di.cc
  n1i.cc <- es$n1i.cc
  n2i.cc <- es$n2i.cc

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

    mu_hat <- mean(log(ci.cc/(n2i.cc-ci.cc)))
    theta_hat <- mean(yi)

    ipm_optim <- function(tau2){
      k <- length(yi)

      sigma2i <- 1/n1i.cc*(exp(-mu_hat-theta_hat+tau2/2)+2+exp(mu_hat+theta_hat+tau2/2))+1/n2i.cc*(exp(-mu_hat)+2+exp(mu_hat))

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
    metafor_method <- method
    tau2 <- NULL
  }

  # run metafor
  if(udweights == FALSE){
    fit <- metafor::rma(yi = yi, vi = vi,
                        method = metafor_method,
                        weighted = weighted,
                        test = test,
                        tau2 = tau2,
                        to = "none", # prevent application of further continuity corrections
                        level = level,
                        digits = digits, verbose = verbose,
                        control = control)
  }

  if(udweights == TRUE){
    fit <- metafor::rma(yi = yi, vi = vi,
                        method = metafor_method,
                        weighted = weighted,
                        weights = weights,
                        test = test,
                        tau2 = tau2,
                        to = "none", # prevent application of further continuity corrections
                        level = level,
                        digits = digits, verbose = verbose,
                        control = control)
  }


  mf <- match.call()

  # make results list
  # UNDER CONSTRUCTION
  res <- list(
    # model information:
    model = "rareIV",
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
    #incl.studies = !attr(es, "remove"),
    # effect sizes and sampling variances:
    yi = yi,
    vi = vi,
    # model matrix:
    X = fit$X,
    # counts and sample sizes::
    ai = ai,
    bi = bi,
    ci = ci,
    di = di,
    ai.cc = ai.cc,
    bi.cc = bi.cc,
    ci.cc = ci.cc,
    di.cc = di.cc,
    ni = n1i+n2i,
    n1i = n1i,
    n2i = n2i,
    n1i.cc = n1i.cc,
    n2i.cc = n2i.cc,
    # continuity corrections:
    tcc = es$tcc,
    ccc = es$ccc,
    cc.studies = es$ccstudies,
    # arguments:
    measure = measure,
    method = method,
    cc = cc,
    ccval = ccval,
    tccval = tccval,
    cccval = cccval,
    ccsum = ccsum,
    ccto = ccto,
    drop00 = drop00,
    weighted = weighted,
    level = level,
    test = test,
    digits = digits,
    # control = control,
    # package version and call:
    version = utils::packageVersion("raremeta"),
    call = mf
  )

  # apply class "raremeta"
  res <- raremeta(res)

  return(res)
}
