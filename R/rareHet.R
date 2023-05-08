#' Testing for between-study homogeneity in meta-analyses of rare events
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
#' @param method  character string specifying whether the test shall be conducted in  a
#' fixed- or a random-effects framework. The test is based on a fixed-effects model when specifying `method = "FE"` .
#' The test is carried out in a random-effects framework if `method` equal to one of the following: `"DL"`, `"HE"`,
#' `"SJ"`, `"ML"`, `"REML"`, `"EB"`, `"HS"`, `"PM"`, `"IPM"`, `"GENQ"`, `"PMM"` or `"GENQM"`.
#' @param test type of test for between-study homogeneity, either `"q"` (for a Q test),
#' `"lrt"` (for a likelihood-ratio test), `"score"` (for a score test), `"wald"` (for a Wald test).
#' Default is `q`.
#' @param type specific type of test to be used. For `test = "q"`,
#' `type` can be set to `"standard"`, `"modified"`, `"jackson1", "jackson2"`, `"gart"`, `"bliss"`, or
#' `"kd"` (short for Kulinskaya-Dollinger). For `test = "lrt"`, `test = "score"` or `test = "wald"`,
#' `method` can be set to `"ML"` or `"REML"`. See Details.
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
#' ## Data input
#' Data input can happen either through the parameter `x` (an object of type `rareData`)
#' or through the parameters `ai`,`bi`,`ci`, `di`, `n1i`, `n2i`, `data` (columns of a dataframe).
#' A `rareData` object can be produced from a data frame by applying the `rareDescribe()` function to it.
#' The `rareDescribe()` function pre-processes the data frame and stores the information required by the `rareES()` function
#' in a list. See `?rareDescribe` for more details.
#'
#' @references
#' Zhang, C., Chen, M., & Wang, X. (2020). Statistical methods for quantifying between-study heterogeneity in meta-analysis
#' with focus on rare binary events. Statistics and its interface, 13(4), 449. doi: 10.4310/sii.2020.v13.n4.a3
#'
#' @examples
rareHet <- function(x, ai, bi, ci, di, n1i, n2i, data,
                    measure = "logOR", method = "FE",
                    test = "q", type = "standard",
                    approx = FALSE,
                    cc, ccval = 0.5, tccval, cccval, ccsum = 1,
                    ccto = "only0", drop00 = TRUE) {

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


  if (!is.element(class(x), c("rareData", "raremeta"))) {
    stop("x must be either of class 'rareData', or of class 'raremeta'.")
  }

  if (test == "q" & !is.element(type, c("standard", "gart", "modified", "bhaumik", "jackson1", "jackson2", "bliss", "kd"))) {
    stop("For `type = 'q'`, type must be chosen from one of the following:
         'standard', 'modified', 'jackson1', 'jackson2', 'bhaumik', 'gart', 'bliss', 'kd'.")
  }

  if (test == "lrt" & is.element(type, c("ML", "REML", "uncondFE", "uncondRE", "cond"))) {
    stop("For `type = 'lrt'`, type must be chosen from one of the following:
         'ML', 'REML', 'uncondFE', 'uncondRE', 'cond'.")
  }

  if (test == "score" & is.element(type, c("ML", "REML", "uncondFE", "uncondRE", "cond", "BD"))) {
    stop("For `type = 'score'`, type must be chosen from one of the following:
         'ML', 'REML', 'uncondFE', 'uncondRE', 'cond', 'BD'.")
  }

  if (test == "wald" & is.element(type, c("ML", "REML", "FE"))) {
    stop("For `type = 'wald'`, method must be chosen from one of the following:
         'ML', 'REML', 'FE'.")
  }

  if(test == "q" & type == "gart"){
    warning("You can obtain the Q-test according to Gart (1966) with the following settings:
            `test = 'q'`, `type = 'standard'`, `method = 'FE'`, `ccto = 'constant', `ccval = 0.5`, and `to = 'all'.
             Results were obtained by resetting the arguments to these values.")

    test <- "q"; type = "standard"; method = "FE"; cc = "constant"; ccto = "all"; ccval = 0.5
  }

  # check whether tccval and cccval are specified; if not, set them to ccval:
  if (missing(tccval)) {
    tccval <- ccval
    cccval <- ccval
  }

  # fit fixed-effects inverse variance model
  # (such that results can be used in the calculation of the Q-Test)
  fit <- rareIV(x,
    measure = measure, method = method,
    cc = cc, ccval = ccval, tccval = tccval, cccval = cccval, ccsum = ccsum,
    ccto = ccto,
    drop00 = drop00
  )

  yi <- fit$yi
  vi <- fit$vi
  beta <- as.numeric(fit$beta)
  tau2 <- as.numeric(fit$tau2)
  k <- fit$k

  ai <- fit$ai
  ci <- fit$ci
  n1i <- fit$n1i
  n2i <- fit$n2i
  ni <- n1i+n2i

  incl.studies <- fit$incl.studies

  if (test == "q") {
    if (type == "standard") {
      wi <- 1 / vi
      betaFE <- sum(wi*yi)/sum(wi)
      stat <- sum(wi * (yi - betaFE)^2)
      pval <- stats::pchisq(stat, k - 1, lower.tail = FALSE)
    }

    if (type == "modified") {
      vi_mod <- (n1i - 16 * (1 - 2 * ai / n1i)^2) / n1i * (1 / (ai + 0.5) + 1 / (n1i - ai + 0.5)) + (n2i - 16 * (1 - 2 * ci / n2i)^2) / n2i * (1 / (ci + 0.5) + 1 / (n2i - ci + 0.5))
      vi_mod <- vi_mod[incl.studies]

      if(any(vi_mod < 0)){
        stop("Some of the modified sampling variances are negative. Please choose a different test and/or type.")
      }

      wi_mod <- 1 / (vi_mod+tau2)

      stat <- sum(wi_mod*(yi-beta)^2)
      pval <-stats::pchisq(stat, k - 1, lower.tail = FALSE)
    }

    if(type == "jackson1"){

      if(method != "FE"){
        warning("You are using a method which is not recommended for this type of test.
                Consider choosing `method = 'FE'` for reliable results.")
      }

      wi <- 1 / (vi+tau2)
      stat <- sum(wi*(yi-beta)^2)
      W <- diag(wi)
      w <- matrix(wi, ncol = 1)
      wsum <- sum(wi)
      A <- W-(1/wsum)*w%*%t(w)

      Y <- matrix(yi, ncol = 1)
      Q <- t(Y)%*%A%*%Y

      Sigma2 <- diag(vi)
      S <- sqrt(Sigma2)%*% A %*% sqrt(Sigma2)
      lambda <- eigen(S)$values

      pval <- CompQuadForm::farebrother(Q, lambda[-k])$Qq
    }

    if(type == "jackson2"){

      if(method != "FE"){
        warning("You are using a method which is not recommended for this type of test.
                Consider choosing `method = 'FE'` for reliable results.")
      }

      wi <- 1 / sqrt(vi+tau2)
      stat <- sum(wi*(yi-beta)^2)
      W <- diag(wi)
      w <- matrix(wi, ncol = 1)
      wsum <- sum(wi)
      A <- W-(1/wsum)*w%*%t(w)

      Y <- matrix(yi, ncol = 1)
      Q <- t(Y)%*%A%*%Y

      Sigma2 <- diag(vi)
      S <- sqrt(Sigma2)%*% A %*% sqrt(Sigma2)
      lambda <- eigen(S)$values

      pval <- CompQuadForm::farebrother(Q, lambda[-k])$Qq
    }

    if(type == "bliss"){
      wi <- 1 / vi
      QFE <- sum(wi * (yi - beta)^2)

      nn <- sum(ni-2)/k
      stat <- k-1 + sqrt((nn-4)/(nn-1))*((nn-2)/nn*QFE-(k-1))
      pval <-stats::pchisq(stat, k - 1, lower.tail = FALSE)
    }

    if(type == "bhaumik"){
      wi = 1/vi
      betaFE = sum(wi*yi)/sum(wi)
      QFE = sum(wi*(yi-betaFE)^2)
      Var.QFE <- 2*sum(wi*(vi+1/sum(wi)-2*wi*vi/sum(wi))^2)
      Var.lnQ <- Var.QFE*(1/(k-1))^2

      wi <- 1/(vi+tau2)
      QRE <- sum(wi*(yi-beta)^2)
      stat <- (k-1)*(log(QRE)-log(k-1))/sqrt(Var.lnQ)
      pval <- stats::pnorm(stat, lower.tail = F) # double-check
    }




  }

  res <- list(test = test,
              type = type,
              measure = measure,
              method = method,
              stat = stat,
              pval = pval)

  return(res)

}
