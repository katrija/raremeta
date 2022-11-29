#' Testing for between-study homogeneity in meta-analyses of rare events
#'
#' @param x an object of class `"rareData"`, or `"raremeta"`.
#' @param measure character string specifying the effect size or outcome measure to be used
#' (either `"logOR"` for the log odds ratio, `"logRR"` for the log relative risk,
#' or `"RD"` for the risk difference).
#' @param method  character string specifying whether the test shall be conducted in  a
#' fixed- or a random-effects framework. The test is based on a fixed-effects model when specifying `method = "FE"` .
#' The test is carried out in a random-effects framework if `method` equal to one of the following: `"DL"`, `"HE"`,
#' `"SJ"`, `"ML"`, `"REML"`, `"EB"`, `"HS"`, `"PM"`, `"IPM"`, `"GENQ"`, `"PMM"` or `"GENQM"`.
#' @param test type of test for between-study homogeneity, either `"q"` (for a Q test),
#' `"lrt"` (for a likelihood-ratio test), `"score"` (for a score test), `"wald"` (for a Wald test).
#' Default is `q`.
#' @param type specific type of test to be used. For `test = "q"`,
#' `type` can be set to `"standard"`, `"modified"`, `"jackson"`, `"gart"`, `"bliss"`, or
#' `"kd"` (short for Kulinskaya-Dollinger). For `test = "lrt"`, `test = "score"` or `test = "wald"`,
#' `method` can be set to `"ML"` or `"REML"`.
#' @param approx logical to specify whether an approximate distribution shall be used if possible
#' @param cc character string specifying the type of continuity corrections to be used
#' (either `"constant"`, `"tacc"` or `"empirical"`). Default is "constant".
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
#' (this option is available for `type = "score"` and `method = "cond"` or `method = "BD"`).
#' @return a list including the value of the test statistic and the associated **p** value.
#'
#' @references
#' Zhang, C., Chen, M., & Wang, X. (2020). Statistical methods for quantifying between-study heterogeneity in meta-analysis
#' with focus on rare binary events. Statistics and its interface, 13(4), 449. doi: 10.4310/sii.2020.v13.n4.a3
#'
#' @export
#'
#' @examples
#'
rareHet <- function(x, measure = "logOR", method = "FE",
                    test = "q", type = "standard",
                    approx = FALSE,
                    cc, ccval = 0.5, tccval, cccval, ccsum = 1,
                    ccto = "only0", drop00 = TRUE){

  if(missing(x)){
    stop("x must be specified.")
  }

  if(!is.element(class(x), c("rareData", "raremeta"))){
    stop("x must be either of class 'rareData', or of class 'raremeta'.")
  }

  if(test == "q" & !is.element(type, c("standard", "modified", "jackson", "bliss", "kd"))){
    stop("For `type = 'q'`, type must be chosen from one of the following:
         'standard', 'modified', 'jackson', 'bhaumik', 'gart', 'bliss', 'kd'." )
  }

  if(test == "lrt" & is.element(type, c("ML", "REML", "uncondFE", "uncondRE", "cond"))){
    stop("For `type = 'lrt'`, type must be chosen from one of the following:
         'ML', 'REML', 'uncondFE', 'uncondRE', 'cond'.")
  }

  if(test == "score" & is.element(type, c("ML", "REML", "uncondFE", "uncondRE", "cond", "BD"))){
    stop("For `type = 'score'`, type must be chosen from one of the following:
         'ML', 'REML', 'uncondFE', 'uncondRE', 'cond', 'BD'.")
  }

  if(test == "wald" & is.element(type, c("ML", "REML", "FE"))){
    stop("For `type = 'wald'`, method must be chosen from one of the following:
         'ML', 'REML', 'FE'.")
  }

  # check whether tccval and cccval are specified; if not, set them to ccval:
  if(missing(tccval)){
    tccval <- ccval
    cccval <- ccval
  }

  # fit inverse variance model (such that results can be used in the calculation of the Q-Test)
  ivfit <- rareIV(x, measure = measure, method = method,
                  cc = cc, ccval = ccval, tccval = tccval, cccval = cccval, ccsum = 1,
                  ccto = ccto,
                  drop00 = drop00)

  yi <- ivfit$yi
  vi <- ivfit$vi
  wi <- 1/vi
  beta <- ivfit$beta
  k <- ivfit$k

  ai <- ivfit$ai
  ci <- ivfit$ci
  n1i <- ivfit$n1i
  n2i <- ivfit$n2i

  if(test == "q"){

    if(type == "standard"){
      Qval <- sum(wi*(yi-beta)^2)
      pval <- stats::pchisq(Qval, k-1, lower.tail = FALSE)
    }

    if(type == "modified"){
      vi_mod <- (n1i-16*(1-2*ai/n1i)^2)/n1i*(1/(ai+0.5)+1/(n1i-ai+0.5))+(n2i-16*(1-2*ci/n2i)^2)/n2i*(1/(ci+0.5)+1/(n2i-ci+0.5))
      wi_mod <- 1/vi_mod
      }


  }



}
