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
                   drop00 = TRUE, test="z", digits, verbose=FALSE, control,
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
  if(!(measure %in% c("logOR", "logRR", "RD"))){
    stop("'measure' must be either 'logOR', 'logRR', or 'RD'.")
  }

  # check if method argument is specified
  if(missing(method)){
    stop("'method' argument must be specified.")
  }

  # check if method argument is valid

  # check if cc is specified (in case there are zero-studies)
  if(any(c(x$ai, x$bi, x$ci, x$di) == 0) & missing(cc)){
    stop("Some studies have zero events. \n
         You must specify the 'cc' argument to determine how they are handled.\n
         If you want to exclude all zero-studies, set 'cc' equal to 'none'.")
  }

  # check if cc argument is valid
  if(!(cc %in% c("constant", "reciprocal", "empirical"))){
    stop("'cc' must be either 'constant', 'reciprocal', or 'empirical'")
  }


  # check if ccto argument is valid
  if(!(ccto %in% c("only0", "all", "if0all"))){
    stop("'ccto' must be either 'only0', 'all', or 'if0all'")
  }

  # check if ccval argument is valid
  if(drop00 == FALSE){
    if(!(length(ccval) %in% c(1,x$k))){
      stop("'ccval' must have length 1 or length equal to the number of studies.")
    }
  }

  if(drop00 == TRUE){
    if(!(length(ccval) %in% c(1,x$k,x$k-x$kdz))){
      stop("'ccval' must have length 1 or length equal to the number of studies (in- or excluding double-zero studies).")
    }
  }

  # check if test = "knha". if so, set test to "hksj"
  if(test == "knha"){
    warning("`test = 'knha'` is no valid option. Setting `test = 'hksj'` instead.")
    test = "hksj"
  }

  # check if test argument is valid
  if(!(test %in% c("z", "hksj"))){
    stop("'test' must be either 'z' or 'hksj'")
  }

  # check if digits argument is valid
  if(length(digits) != 1 | digits%%1 != 0){
    stop("'digits' must be an integer of length 1.")
  }

  # apply continuity correction (and drop double-zero studies if drop00 = TRUE)
  ai_cc <- x$ai
  bi_cc <- x$bi
  ci_cc <- x$ci
  di_cc <- x$di

  # convert measure to the metafor notation
  metafor_measure <- sub("log", "", measure)

  # run metafor
  res <- metafor::rma(ai = ai_cc, bi = bi_cc,
                      ci = ci_cc, di = di_cc,
                      measure = metafor_measure, method = method,
                      test = test,
                      to = "none", # prevent application of further continuity corrections
                      digits = digits, verbose = verbose)


  return(res)
}
