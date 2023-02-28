#' Calculate effect sizes for meta-analyses of rare events.
#'
#' @param x an object of class `"rareData"`.
#' @param measure character string specifying the effect size or outcome measure to be used
#' (either `"logOR"` for the log odds ratio, `"logRR"` for the log relative risk,
#' or `"RD"` for the risk difference).
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
#' @param drop00 logical indicating whether double-zero studies (i.e., studies with no events or
#' only events in both groups) should be excluded when calculating the outcome measufit for the
#' individual studies.
#' @param method Needs to be specified when `cc = "empirical"`. Meta-analytic method to be used when
#' determining the prior used for the empirical continuity correction. Current default is `method = "FE"`.
#' See Sweeting et al. (2004) for details.
#' @param ... additional arguments
#'
#' @details
#' # Details
#' ## Data input
#' The main input of the `rareES()` function is a so-called `rareData` object. A `rareData` object
#' can be produced from a data frame by applying the `rareDescribe()` function to it. The `rareDescribe()`
#' function pre-processes the data frame and stores the information required by the `rareES()` function
#' in a list. See `?rareDescribe` for more details.
#' ## Effect size measures
#' The function includes methods for calculating log odds ratios, log risk ratios, and risk differences.
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
#' If the data include at least one single-zero study, the function throws an error if the user
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
#' @references
#' Gart, John J, and James R Zweifel. 1967. On the bias of various estimators of the logit and
#' its variance with application to quantal bioassay. Biometrika, 54, 181–187. doi:10.1093/BIOMET/54.1-2.181
#'
#' Sweeting, M. J., Sutton, A. J., & Lambert, P. C. (2004). What to add to nothing? Use and avoidance of
#' continuity corrections in meta-analysis of sparse data. Statistics in Medicine, 23, 1351–1375. doi: 10.1002/sim.1761
#'
#' @return a matrix with eight columns. The first six columns (`"ai.cc"`, `"bi.cc"`, `"ci.cc"`, `"di.cc"`, `"n1i.cc"`, `"n2i.cc"`)
#' contain the continuity corrected counts and sample sizes. The last two columns (`"yi"` and `"vi"`) contain
#' the effect sizes and sampling variances, respectively. Information on how the effect sizes were obtained (in particular
#' with regard to the use of continuity corrections) is stored in attributes.
#'
#' @export
#'
#' @examples
#'
#' data <- data.frame(
#' ai = c(0, 3, 2, 0),
#' bi = c(20, 18, 15, 19),
#' ci = c(1, 4, 0, 0),
#' di = c(19, 17, 16, 20)
#' )
#'
#' x <- rareDescribe(ai = ai, bi = bi, ci = ci, di = di, data = data)
#'
#' rareES(x, measure = "logOR", cc = "constant", ccval = 0.5)
rareES <- function(x, measure, cc, ccval = 0.5, tccval, cccval, ccsum = 1,
                   ccto = "only0", drop00 = TRUE,
                   method = "FE", ...){

  # check if x is an object of class rareData
  if(!inherits(x,"rareData")){
    stop("x must be an object of class 'rareData'. See ?rareDescribe for more details.")
  }

  # apply continuity correction if needed (argument checks are done in rareCC)
  x <- rareCC(x, cc = cc, ccval = ccval, tccval = tccval, cccval = cccval,
              ccsum = ccsum, ccto = ccto, drop00 = drop00, measure = measure,
              method = method)

  # extract counts and sample sizes
  ai.cc   <- x[,"ai.cc"]
  bi.cc   <- x[,"bi.cc"]
  ci.cc   <- x[,"ci.cc"]
  di.cc   <- x[,"di.cc"]
  n1i.cc  <- x[,"n1i.cc"]
  n2i.cc  <- x[,"n2i.cc"]

  # calculate effect sizes and sampling variances:
  if(measure == "logOR"){
    yi <- log(ai.cc*di.cc/(bi.cc*ci.cc))
    vi <- 1/ai.cc+1/bi.cc+1/ci.cc+1/di.cc
  }

  if(measure == "logRR"){
    yi <- log((ai.cc/n1i.cc)/(ci.cc/n2i.cc))
    vi <- 1/ai.cc-1/n1i.cc+1/ci.cc-1/n2i.cc
  }

  if(measure == "RD"){
    yi <- ai/n1i - ci/n2i
    vi <- (ai.cc*(n1i.cc-ai.cc))/(n1i.cc^3)+(ci.cc*(n2i.cc-ci.cc))/(n2i.cc^3)
  }

  out <- cbind(x, yi, vi)

  #adding attributes
  #attr(out, "measure") <- measure
  #attr(out, "cc") <- cc
  #attr(out, "ccto") <- ccto
  #attr(out, "drop00") <- drop00

  ##The following attributes are local to the rareCC function
  ##attr(out, "ccstudies") <- ccstudies
  ##attr(out, "ccc") <- ccc
  ##attr(out, "tcc") <- tcc
  ##attr(out, "remove") <- remove

  #colnames(out) <- c("ai.cc", "bi.cc", "ci.cc", "di.cc", "n1i.cc", "n2i.cc", "yi", "vi")


  #Output
  return(out)
}
