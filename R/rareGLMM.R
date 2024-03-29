#' Conduct a meta-analysis using a generalized linear (mixed) model (GLMM)
#'
#'@description
#' Function to conduct a meta-analysis of a rare event using a generalized linear
#' or generalized linear mixed model. See below for more details on these
#' models and their application in meta-analyses of rare events.
#'
#' @param x an object of class `"rareData"`.
#' @param ai data frame column to specify the number of events in group 1 (i.e., the treatment group).
#' @param bi data frame column to specify the number of non-events in group 1 (i.e., the treatment group).
#' @param ci data frame column to specify the number of events in group 2 (i.e., the control group).
#' @param di data frame column to specify number of non-events in group 2 (i.e., the control group).
#' @param n1i data frame column to specify the sample sizes in group 1 (i.e., the treatment group).
#' @param n2i data frame column to specify the sample sizes in group 2 (i.e., the control group).
#' @param data data frame.
#' @param measure character string specifying the effect size or outcome measure to be used
#' (either `"logOR"` for the log odds ratio or `"logRR"` for the log relative risk). See below for more details.
#' @param intercept character string specifying whether to fit a model with fixed, study-specific
#' intercepts (`"fixed"`) or a random intercept (`"random"`). See below for more details.
#' @param slope character string specifying whether to fit a model with a fixed slope (`"fixed"`) or a
#' random slope (`"random"`). A model with a fixed slope is a model with homogeneous effects, that is, a
#' fixed-effects model, while a model with a random slope is a model with heterogeneous effects, that is, a
#' random effects model. See below for more details.
#' @param conditional logical specifying whether to estimate a conditional generalized linear mixed model.
#' Default is `FALSE`.
#' @param approx logical specifying whether to use the approximate version of the conditional generalized linear
#' mixed model (i.e., the conditional binomial model). Only relevant when `conditional = TRUE` and `measure = "logOR`.
#' Default is `FALSE`. See below for more details.
#' @param cor logical specifying whether random effects should be modeled as correlated or uncorrelated.
#' Default is `cor = FALSE`. This argument is only relevant if `intercept = "random"` and `slope = "random"`.
#' See below for more details.
#' @param coding numeric specifying the coding scheme used for the random effects structure. Values between 0
#' and 1 can be specified. Default is `coding = 1/2`. Given that `cor = FALSE` is specified, the default
#' option implies equal variances in the two groups. Values closer to 0 imply a larger variance in the
#' control group (group 2), while values closer to 1 imply a larger variance in the treatment group (group 1).
#' See below for more details.
#' @param drop00 logical indicating whether double-zero studies (i.e., studies with no events or only events in both groups)
#' should be excluded prior to calculating the studies' effect sizes and sampling variances.
#' @param level numeric between 0 and 100 specifying the confidence interval level (the default is 95).
#' @param test character string specifying how test statistics and confidence intervals for the
#' fixed effects should be computed (currently, only `"z"`, for Wald-type tests is available).
#' @param digits integer specifying the number of decimal places to which the printed results
#' should be rounded (if unspecified, the default is 4).
#' @param verbose logical indicating whether output should be generated on the progress of model
#' fitting (the default is `FALSE`). Can also be an integer. Values > 1 generate more verbose output.
#' @param control optional list of control values for the iterative algorithms. If unspecified, default
#' values are defined inside the functions.
#' @param ... additional arguments.
#'
#' @details
#' \loadmathjax{}
#'
#' ## Data input
#' The data input can be specified either through the arguments `ai`,`bi`,`ci`,`di`,`n1i`, and `n2i` (columns of the data frame `data`)
#' or through the argument `x`, which takes an object that results from applying the `rareDescribe()` function to the data
#' (i.e., the input for argument `x` must be an object of type `rareData`).
#' A `rareData` object can be produced from a data frame by applying the `rareDescribe()` function to it.
#' The `rareDescribe()` function pre-processes the data frame and stores the information required by the `rareGLMM()` function in a list.
#' See `?rareDescribe` for more details.
#'
#' ## Effect size measures
#' The function includes different versions of the generalized linear (mixed) model. Per default, the function fits an
#' unconditional model. For information on fitting conditional generalized linear (mixed) models, see further below.
#'
#' The regression equation for the unconditional fixed-effects model (GLM) can be expressed as
#'
#' \mjseqn{g(\pi_{ij}) = \alpha_i + \theta \cdot x_{ij}},
#'
#' and the model equation for the random-effects model (GLMM) can be expressed as
#'
#' \mjseqn{g(\pi_{ij}) = \alpha_i+ \theta \cdot x_{ij} + \epsilon_{i} \cdot z_{ij}},
#'
#' where \mjseqn{i = 1, ..., k } is the study index and \mjseqn{j = 1, 2} is the group index,
#' \mjseqn{x_{i1} = 1}, \mjseqn{x_{i2} = 0}, and \mjseqn{z_{ij}} is defined by the `coding` argument.
#' Specifically, for `coding = z`, \mjseqn{z_{i1} = z} and \mjseqn{z_{i2} = z-1}. Default is \mjseqn{z = 1/2}.
#' `coding = 1` corresponds to Models 2 and 4 in Jackson et al. (2017), and `coding = 1/2` corresponds to Models 3 and 5 in
#' the same study.
#'
#' For `measure = "logOR"`, a GL(M)M with a binomial within-study distribution is used for the numbers of
#' events in each group, with sample size \mjseqn{n_{ij}} and probability \mjseqn{\pi_{ij}}. For this model,
#' \mjseqn{g(\cdot)} corresponds to the logit function, i.e., \mjseqn{g(\pi_{ij}) = \log \left( \frac{\pi_{ij}}{1-\pi_{ij}} \right)}.
#' See Stijnen et al. (2011), for further details.
#'
#' For `measure = "logRR"`, a GL(M)M with a Poisson within-study distribution is used for the numbers of
#' events in each group, with rate \mjseqn{n_{ij}\pi_{ij}}. For this model,
#' \mjseqn{g(\cdot)} corresponds to the natural logarithm. See Böhning et al. (2015), for
#' further details.
#'
#' It is currently not possible to fit GL(M)Ms for `measure = "RD"`, but this functionality may be
#' enabled in future versions of the package.
#'
#' ## Baseline and effect heterogeneity
#' Baseline heterogeneity refers to the way the intercept is modeled. By specifying
#' `intercept = "fixed"`, the intercepts \mjseqn{\alpha_i} are modeled as fixed, study-specific
#' intercepts. By specifying `intercept = "random"`, it is assumed that \mjseqn{\alpha_i \sim N(\alpha, \sigma_{\alpha}^2)}.
#'
#' Effect heterogeneity refers to the way the slope is modeled. By specifying
#' `slope = "fixed"`, a fixed slope is modeled, corresponding to a fixed-effects meta-analysis.
#' By specifying `slope = "random`, it is assumed that \mjseqn{\epsilon_i \sim N(0, \tau^2)},
#' corresponding to a random-effects meta-analysis with between-study variance \mjseqn{\tau^2}.
#'
#' ## Correlated random effects
#' When both `intercept = "random"` and `slope = "random"` is specified, the argument `cor` can be
#' used to specify whether random-effects are assumed to be correlated. Per default, this is not
#' the case, i.e. `cor = FALSE`.
#'
#' ## Conditional generalized linear (mixed) models
#'
#' When specifying `conditional = TRUE`, a conditional generalized linear (mixed) model is fitted.
#' For `measure = "logOR"` and `intercept = "fixed"`, this is the fixed-effects hypergeometric model described by van Houwelingen et al. (1993).
#' The hypergeometric-normal model described in the same article, which assumes a normal distribution for the true effect sizes (i.e., the intercepts
#' of the hypergeometric model) can be fitted by setting `intercept = "random"`.
#' For `measure = "logRR"`, the conditional binomial(-normal) model described by Cai et al. (2010) is fitted. The same model is fitted
#' when specifying `measure = "logOR` and `approx = TRUE`, as the conditional binomial(-normal) model can be viewed as an approximation
#' to the hypergeometric(-normal) model (see Stijnen et al., 2010, for further details).
#'
#'
#' @return an object of class "raremeta". The object is a list containing the following elements:
#' * `model`: name of the model used for conducting the meta-analysis.
#' * `beta`: estimated coefficients of the model.
#' * `se`: standard errors of the  coefficients.
#' * `zval`: test statistics of the coefficients.
#' * `zval`: p-values corresponding to the test statistics.
#' * `ci.lb`: lower bound of the confidence intervals for the coefficients.
#' * `ci.ub`: upper bound of the confidence intervals for the coefficients.
#' * `vb`: variance-covariance matrix of the estimated coefficients.
#' * `tau2`: estimated amount of (residual) heterogeneity. Always `0` when `method = "FE"`.
#' * `LRT.Chisq`: Test statistic of the likelihood ratio test testing for homogeneity.
#' * `LRT.df`: Degrees of freedom of the likelihood ratio test testing for homogeneity.
#' * `LRT.pval`: p-value of the likelihood ratio test testing for homogeneity.
#' * `fit.stats`: a list with log-likelihood, deviance, AIC, BIC, and AICc values under
#' the unrestricted likelihood.
#' * `p`: number of coefficients in the model (including the intercept).
#' * `k`: number of studies included in the analysis.
#' * `k.all`: total number of studies (before exclusion).
#' * `kdz`,`ksz`: number of double-zero and single-zero studies.
#' * `k1sz`, `k2sz`: number of single-zero studies where the zero is in group 1 or group 2.
#' * `ai`, `bi`, `ci`, `di`: original entries of the 2x2 tables for all studies.
#' * `ni`, `n1i`, `n2i`: original total and group sample sizes.
#' * ...
#' @references
#' Böhning, D., Mylona, K., & Kimber, A. (2015). Meta-analysis of clinical trials with rare
#' events. Biometrical Journal, 57 (4), 633–648. doi: 10.1002/bimj.201400184
#'
#' Cai, T., Parast, L., & Ryan, L. (2010). Meta-analysis for rare events. Statistics in
#' Medicine, 29 (20), 2078–2089. doi: 10.1002/sim.3964
#'
#' Jackson, D., Law, M., Stijnen, T., Viechtbauer, W., & White, I. R. (2018). A comparison
#' of seven random-effects models for meta-analyses that estimate the summary odds
#' ratio. Statistics in Medicine, 37 (7), 1059–1085. doi: 10.1002/sim.7588
#'
#' Stijnen, T., Hamza, T. H., & Özdemir, P. (2010). Random effects meta-analysis of event
#' outcome in the framework of the generalized linear mixed model with applications
#' in sparse data. Statistics in Medicine, 29 (29), 3046–3067. doi: 10.1002/sim.4040
#'
#' van Houwelingen, H. C., Zwinderman, K. H., & Stijnen, T. (1993). A bivariate approach
#' to meta-analysis. Statistics in Medicine, 12 (24), 2273–2284. doi: 10.1002/sim.4780122405
#'
#' @examples
#'
#' # load a dataset
#' data(dat.nissen2007)
#' d <- dat.nissen2007
#'
#' # GLMM for the log odds ratio with fixed intercept
#' rareGLMM(ai=miRosiglitazone, ci=miControl, n1i=nRosiglitazone, n2i=nControl, data=d, measure="logOR", intercept="fixed")
#'
#' # same analysis with pre-processed data
#' x <- rareDescribe(ai = ai, bi = bi, ci = ci, di = di, data = data)
#' rareGLMM(x, measure="logOR", intercept="fixed")
#'
#' @export
#' @import mathjaxr
rareGLMM <- function(x, ai, bi, ci, di, n1i, n2i, data, measure,
                     intercept = "fixed", slope = "random", conditional = FALSE, approx = FALSE,
                     cor = FALSE, coding = 1 / 2,
                     drop00 = FALSE,
                     level = 95,
                     test = "z", digits = 4, verbose = FALSE, control,
                     ...) {

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
  if (!inherits(x, "rareData")) {
    stop("x must be an object of class 'rareData'. See ?rareDescribe for more details.")
  }

  # check if measure argument is specified
  if (missing(measure)) {
    stop("'measure' argument must be specified.")
  }

  # check if measure argument is valid
  if (!is.element(measure, c("logOR", "logRR"))) {
    stop("'measure' must be either 'logOR' or 'logRR'.")
  }

  if (!is.element(slope, c("fixed", "random"))) {
    stop("'slope' must be either 'fixed' or 'random'.")
  }

  if (!is.element(intercept, c("fixed", "random"))) {
    stop("'intercept' must be either 'fixed' or 'random'.")
  }

  if (!is.element(conditional, c(TRUE, FALSE))) {
    stop("'conditional' must be either TRUE or FALSE.")
  }

  if (intercept == "random" & slope == "random" & !is.element(cor, c(TRUE, FALSE))) {
    stop("'cor' must be either TRUE or FALSE.")
  }

  # extract counts and sample sizes
  ai <- x$ai
  bi <- x$bi
  ci <- x$ci
  di <- x$di
  n1i <- x$n1i
  n2i <- x$n2i

  if (!is.logical(drop00)) {
    stop("'drop00' must be a logical")
  }

  # check if level argument is valid:
  if (!is.numeric(level) || length(level) > 1 || level < 0 || level > 100) {
    stop("level must be a scalar between 0 and 100.")
  }

  # check if test argument is valid
  if (!is.element(test, c("z"))) {
    stop("'test' must be 'z'.")
  }

  # check if digits argument is valid
  if (length(digits) != 1 || digits %% 1 != 0 || digits < 0) {
    stop("'digits' must be an integer of length 1.")
  }

  mf <- match.call()

  level <- (100 - level) / 100

  # remove NAs:
  missings <- (is.na(x$ai)|is.na(x$bi)|is.na(x$ci)|is.na(x$di))

  if(sum(missings) > 0){
    warning(paste("Your data contain", sum(missings), "studies with missing values. Those studies were removed from the analysis."))
    ai <- ai[!missings]
    bi <- bi[!missings]
    ci <- ci[!missings]
    di <- di[!missings]
    n1i <- n1i[!missings]
    n2i <- n2i[!missings]
  }

  # remove double-zero studies if desired:
  if (drop00 == TRUE) {
    remove <- (ai == 0 & ci == 0) | (bi == 0 & di == 0)
    ai <- ai[!remove]
    bi <- bi[!remove]
    ci <- ci[!remove]
    di <- di[!remove]
    n1i <- n1i[!remove]
    n2i <- n2i[!remove]
  }


  dataLong <- data.frame(
    y = c(ai, ci),
    n = c(n1i, n2i),
    id = factor(c(1:length(ai), 1:length(ci))),
    group = c(rep(1, length(ai)), rep(0, length(ci)))
  )

  dataLong$groupRE <- ifelse(dataLong$group == 1, coding, coding-1)

  if(conditional == FALSE){
    # Define family --------------------------------------------------------------
    if (measure == "logOR") {
      fam <- stats::binomial(link = "logit")
    }

    if (measure == "logRR") {
      fam <- stats::poisson(link = "log")
    }

    # Define model formula -------------------------------------------------------
    if (measure == "logOR") {
      if (intercept == "random" & slope == "random" & conditional == FALSE) {
        if (cor == FALSE) {
          m <- cbind(y, n - y) ~ 1 + group + (1 + groupRE || id)
        } else {
          m <- cbind(y, n - y) ~ 1 + group + (1 + groupRE | id)
        }
        mFE <- cbind(y, n - y) ~ 1 + group + (1 | id)
        mSAT <- cbind(y, n - y) ~ 1 + group + id:group + (1 | id)
      }

      if (intercept == "fixed" & slope == "random" & conditional == FALSE) {
        m <- cbind(y, n - y) ~ -1 + id + group + (0 + groupRE | id)
        mFE <- cbind(y, n - y) ~ -1 + id + group
        mSAT <- cbind(y, n - y) ~ -1 + id + group + id:group
      }

      if (intercept == "random" & slope == "fixed" & conditional == FALSE) {
        m <- cbind(y, n - y) ~ 1 + group + (1 | id)
      }

      if (intercept == "fixed" & slope == "fixed" & conditional == FALSE) {
        m <- cbind(y, n - y) ~ -1 + id + group
      }
    }

    if (measure == "logRR") {
      if (intercept == "random" & slope == "random" & conditional == FALSE) {
        if (cor == FALSE) {
          m <- y ~ 1 + group + offset(log(n)) + (1 + groupRE || id)
        } else {
          m <- y ~ 1 + group + offset(log(n)) + (1 + groupRE | id)
        }
        mFE <- y ~ 1 + group + offset(log(n)) + (1 | id)
        mSAT <- y ~ 1 + group + id:group + offset(log(n)) +(1 | id)
      }

      if (intercept == "fixed" & slope == "random" & conditional == FALSE) {
        m <- y~ -1 + id + group + offset(log(n)) + (0 + groupRE | id)
        mFE <- y~ -1 + id + group + offset(log(n))
        mSAT <- y ~ -1 + id + group + id:group + offset(log(n))
      }

      if (intercept == "random" & slope == "fixed" & conditional == FALSE) {
        m <- y~1 + group + offset(log(n)) + (1 | id)
      }

      if (intercept == "fixed" & slope == "fixed" & conditional == FALSE) {
        m <- y~ -1 + id + group + offset(log(n))
      }
    }

    # Fit ML model ---------------------------------------------------------------
    if (intercept == "fixed" & slope == "fixed") {
      fitML <- try(
        stats::glm(m,
                   data = dataLong,
                   family = fam
        ), silent = TRUE
      )
    } else {
      fitML <- try(
        lme4::glmer(m,
                    data = dataLong,
                    family = fam
        ), silent = TRUE
      )
    }

    if(inherits(fitML, "try-error")){
      stop("Unable to fit model.")
    }

    llML <- stats::logLik(fitML)

    # Fit FE and SAT model -----------------------------------------------------
    llFE <- LRT.Chisq <- LRT.df <- LRT.pval <- NA

    if (slope == "random") {
      if (intercept == "fixed") {
        fitFE <- try(
          stats::glm(mFE,
                     data = dataLong,
                     family = fam
          ), silent = TRUE
        )

        fitSAT <- try(
          stats::glm(mSAT,
                     data = dataLong,
                     family = fam
          ), silent = TRUE
        )
      } else {
        fitFE <- try(
          lme4::glmer(mFE,
                      data = dataLong,
                      family = fam
          ), silent = TRUE
        )

        fitSAT <- try(
          lme4::glmer(mSAT,
                      data = dataLong,
                      family = fam
          ), silent = TRUE
        )
      }

      # LRT --------------------------------------------------------------------
      if(inherits(fitFE, "try-error")||inherits(fitSAT, "try-error")){
        warning("Unable to fit fixed-effects model or saturated model.\n Results of the Likelihood-ratio test for homogeneity cannot be obtained.")
      }else{
        llFE <- stats::logLik(fitFE)
        llSAT <- stats::logLik(fitSAT)

        LRT.Chisq <- as.numeric(2 * (llSAT - llFE))
        LRT.df <- attributes(llSAT)$df - attributes(llFE)$df
        LRT.pval <- stats::pchisq(LRT.Chisq, df = LRT.df, lower.tail = FALSE)
      }
    }

    # Output generation ----------------------------------------------------------

    if (inherits(fitML, "glmerMod") & slope == "random") {
      beta <- lme4::fixef(fitML)
      vb <- as.matrix(stats::vcov(fitML))
      sigma2 <- lme4::VarCorr(fitML)

      if(cor == FALSE){
        tau2 <- sigma2[[length(sigma2)]][1]
      }else{
        tau2 <- sigma2[[1]][2,2]
      }

      singular <- lme4::isSingular(fitML)
      conv <- ifelse(fitML@optinfo$conv$opt == 0, TRUE, FALSE)
    }

    if (inherits(fitML, "glmerMod") & slope == "fixed") {
      beta <- lme4::fixef(fitML)
      vb <- as.matrix(stats::vcov(fitML))
      sigma2 <- lme4::VarCorr(fitML)
      tau2 <- 0

      singular <- lme4::isSingular(fitML)
      conv <- ifelse(fitML@optinfo$conv$opt == 0, TRUE, FALSE)
    }

    if (inherits(fitML, "glm")) {
      beta <- fitML$coefficients
      vb <- as.matrix(stats::vcov(fitML))
      sigma2 <- NA
      tau2 <- 0

      conv <- fitML$converged
      singular <- NA
    }


    se <- sqrt(diag(vb))

    zval <- beta / se
    pval <- 2 * stats::pnorm(abs(zval), lower.tail = FALSE)

    zcrit <- stats::qnorm(1 - level / 2)
    ci.lb <- beta - zcrit * se
    ci.ub <- beta + zcrit * se

    X <- stats::model.matrix(fitML)
    p <- ifelse(all(X[1, ] == 1), ncol(X) - 1, ncol(X))

    AICc <- -2 * llML + 2 * (p + 1) * max(2 * nrow(dataLong), p + 3) / (max(nrow(dataLong), 3) - p)
    fit.stats <- rbind(llML, stats::deviance(fitML), stats::AIC(fitML), stats::BIC(fitML), AICc)
    rownames(fit.stats) <- c("llML", "dev", "AIC", "BIC", "AICc")
    colnames(fit.stats) <- c("ML")
  }

  if(conditional == TRUE){

    llFE <- LRT.Chisq <- LRT.df <- LRT.pval <- NA
    sigma2 <- singular <- NA

    if(measure == "logOR" && approx == FALSE){

      if(intercept == "fixed"){

        fitMH <- rareMH(x, measure = "logOR")
        parms <- c(fitMH$beta)

        fit <- try(stats::nlminb(parms, objective = .negllnchg, gradient = NULL, hessian = NULL,
                                 ai = ai, bi = bi, ci = ci, di = di, intercept = "fixed"), silent = TRUE)

        if(inherits(fit, "try-error")){
          stop("Unable to fit model.")
        }

        conv <- ifelse(fit$convergence == 0, TRUE, FALSE)

        beta <- fit$par
        names(beta) <- measure

        tau2 <- NA

        hessian <- numDeriv::hessian(.negllnchg, fit$par, ai = ai, bi = bi, ci = ci, di = di, intercept = "fixed")
        vb <- try(solve(hessian), silent = TRUE)

        if(inherits(vb, "try-error")){
          warning("Standard errors could not be obtained.")
          se <- zval <- pval <- zcrit <- ci.lb <- ci.ub <- NA
        }else{
          se <- sqrt(diag(vb))
          zval <- beta/se
          pval <- 2 * stats::pnorm(abs(zval), lower.tail = FALSE)

          zcrit <- stats::qnorm(1 - level / 2)
          ci.lb <- beta - zcrit * se
          ci.ub <- beta + zcrit * se
        }

        llML <- -fit$objective
        p <- 1
        X <- as.matrix(rep(1, length(ai)))

        #AICc <- -2 * llML + 2 * (p + 1) * max(2 * length(ai), p + 3) / (max(length(ai), 3) - p)
        fit.stats <- rbind(llML, NA, NA, NA, NA )
        rownames(fit.stats) <- c("llML", "dev", "AIC", "BIC", "AICc")
        colnames(fit.stats) <- c("ML")


      }


      if(intercept == "random"){

        # Versuch eigener Implementierung (work in progress):

        # fitMH <- rareMH(x, measure = "logOR")
        # parms <- c(fitMH$beta)

        # fitFE <- try(stats::nlminb(parms, objective = .negllnchg, gradient = NULL, hessian = NULL,
        #                            ai = ai, bi = bi, ci = ci, di = di, intercept = "fixed"), silent = TRUE)

        # fitIV <- rareIV(x, measure = "logOR", method = "IPM", cc = "constant", ccval = 0.5, ccto = "all", drop00 = FALSE)
        # parms <- c(fitMH$beta, fitIV$tau2)

        # fitML <- try(stats::nlminb(par = parms, fn = .negllnchg, method = "BFGS",
        #                           ai = ai, bi = bi, ci = ci, di = di, intercept = "random"),
        #              silent = TRUE)

        fitML <- try(metafor::rma.glmm(ai = ai, bi = bi, ci = ci, di = di,
                                       model = "CM.EL", method = "ML", measure = "OR",
                                       control = list(optimizer = "BFGS")), silent = TRUE)


        if(inherits(fitML, "try-error")){
          stop("Unable to fit model.")
        }

        conv <- NA
        # conv <- ifelse(fitML$convergence == 0, TRUE, FALSE)

        beta <- fitML$beta
        # beta <- fitML$par[1]
        names(beta) <- measure

        tau2 <- fitML$tau2
        # tau2 <- fitML$par[2]

        # hessian <- numDeriv::hessian(.negllnchg, fitML$par, ai = ai, bi = bi, ci = ci, di = di, intercept = "random")
        vb <- fitML$vb

        # if(inherits(vb, "try-error")){
        #   warning("Standard errors could not be obtained.")
        #   se <- zval <- pval <- zcrit <- ci.lb <- ci.ub <- NA
        # }else{
          se <- fitML$se
          # se <- sqrt(diag(vb))[1]
          zval <- beta/se
          pval <- 2 * stats::pnorm(abs(zval), lower.tail = FALSE)

          zcrit <- stats::qnorm(1 - level / 2)
          ci.lb <- beta - zcrit * se
          ci.ub <- beta + zcrit * se
       #  }

        # llML <- -fitML$objective
        p <- 1
        X <- as.matrix(rep(1, length(ai)))

        fit.stats <- matrix(fitML$fit.stats[,"ML"], ncol = 1)
        rownames(fit.stats) <- c("llML", "dev", "AIC", "BIC", "AICc")
        colnames(fit.stats) <- c("ML")

        llML <- fit.stats[1,1]

        # if(inherits(fitFE, "try-error")){
        #   warning("Unable to fit fixed-effects model")
        # }else{

          llFE <- NA
          # llFE <- -fitFE$objective

          LRT.Chisq <- fitML$QE.LRT
          LRT.df <- fitML$QE.df

          # LRT.df <- length(fitML$par)-length(fitFE$par)
          LRT.pval <- stats::pchisq(LRT.Chisq, df = LRT.df, lower.tail = FALSE)
       # }
      }


    }

    if(measure == "logRR" || approx == TRUE){

      off <- log(n1i/n2i)

      if(intercept == "fixed"){

        fit <- try(stats::glm(cbind(ai, ci)~1, offset = off,
                          family = binomial(link = "logit")), silent = TRUE)

        if(inherits(fit, "try-error")){
          stop("Unable to fit model.")
        }

        conv <- fit$converged

        beta <- fit$coefficients
        names(beta) <- measure

        tau2 <- NA

        vb <- as.matrix(stats::vcov(fit))

        if(inherits(vb, "try-error")){
          warning("Standard errors could not be obtained.")
          se <- zval <- pval <- zcrit <- ci.lb <- ci.ub <- NA
        }else{
          se <- sqrt(diag(vb))[1]
          zval <- beta/se
          pval <- 2 * stats::pnorm(abs(zval), lower.tail = FALSE)

          zcrit <- stats::qnorm(1 - level / 2)
          ci.lb <- beta - zcrit * se
          ci.ub <- beta + zcrit * se
        }

        llML <- stats::logLik(fit)

        X <- stats::model.matrix(fit)
        p <- ifelse(all(X[1, ] == 1), ncol(X) - 1, ncol(X))

        AICc <- -2 * llML + 2 * (p + 1) * max(2 * length(ai), p + 3) / (max(length(ai), 3) - p)
        fit.stats <- rbind(llML, stats::deviance(fit), stats::AIC(fit), stats::BIC(fit), AICc)
        rownames(fit.stats) <- c("llML", "dev", "AIC", "BIC", "AICc")
        colnames(fit.stats) <- c("ML")
      }


      if(intercept == "random"){
        study <- 1:length(ai)

        fitFE <- try(stats::glm(cbind(ai, ci)~1, offset = off,
                            family = binomial(link = "logit")), silent = TRUE)

        fitML <- try(lme4::glmer(cbind(ai, ci)~1+(1|study), offset = off,
                             family = binomial(link = "logit")), silent = TRUE)

        fitSAT <- try(stats::glm(cbind(ai, ci)~1+factor(study), offset = off,
                                 family = binomial(link = "logit")), silent = TRUE)

        if(inherits(fitML, "try-error")){
          stop("Unable to fit model.")
        }

        conv <- ifelse(fitML@optinfo$conv$opt == 0, TRUE, FALSE)

        beta <- lme4::fixef(fitML)
        names(beta) <- measure

        sigma2 <- lme4::VarCorr(fitML)
        tau2 <- sigma2[[1]][1]

        vb <- as.matrix(stats::vcov(fitML))

        if(inherits(vb, "try-error")){
          warning("Standard errors could not be obtained.")
          se <- zval <- pval <- zcrit <- ci.lb <- ci.ub <- NA
        }else{
          se <- sqrt(diag(vb))[1]
          zval <- beta/se
          pval <- 2 * stats::pnorm(abs(zval), lower.tail = FALSE)

          zcrit <- stats::qnorm(1 - level / 2)
          ci.lb <- beta - zcrit * se
          ci.ub <- beta + zcrit * se
        }

        llML <- stats::logLik(fitML)

        X <- stats::model.matrix(fitML)
        p <- ifelse(all(X[1, ] == 1), ncol(X) - 1, ncol(X))

        AICc <- -2 * llML + 2 * (p + 1) * max(2 * length(ai), p + 3) / (max(length(ai), 3) - p)
        fit.stats <- rbind(llML, stats::deviance(fitML), stats::AIC(fitML), stats::BIC(fitML), AICc)
        rownames(fit.stats) <- c("llML", "dev", "AIC", "BIC", "AICc")
        colnames(fit.stats) <- c("ML")

        if(inherits(fitFE, "try-error")||inherits(fitSAT, "try-error")){
          warning("Unable to fit fixed-effects model or saturated model.\n Results of the Likelihood-ratio test for homogeneity cannot be obtained.")
        }else{
          llFE <- stats::logLik(fitFE)
          llSAT <- stats::logLik(fitSAT)
          LRT.Chisq <- as.numeric(2 * (llSAT - llFE))
          LRT.df <- attributes(llSAT)$df - attributes(llFE)$df
          LRT.pval <- stats::pchisq(LRT.Chisq, df = LRT.df, lower.tail = FALSE)
        }

      }
    }
  }

  # make results list
  # UNDER CONSTRUCTION
  res <- list(
    # model information:
    model = "rareGLMM",
    b = beta,
    beta = beta,
    se = se,
    zval = zval,
    pval = pval,
    ci.lb = ci.lb,
    ci.ub = ci.ub,
    vb = vb,
    sigma2 = sigma2,
    tau2 = tau2,
    # LRT:
    LRT.Chisq = LRT.Chisq,
    LRT.df = LRT.df,
    LRT.pval = LRT.pval,
    # se.tau2 = fit$se.tau2,
    # I2 = fit$I2,
    # H2 = fit$H2,
    # R2 = fit$R2,
    # vt = fit$vt,
    # QE = fit$QE,
    # QEp = fit$QEp,
    # QM = fit$QM,
    # QMdf = fit$QMdf,
    # QMp = fit$QMp,
    fit.stats = fit.stats,
    p = p,
    # convergence information:
    conv = conv,
    singular = singular,
    # study numbers:
    k = nrow(dataLong) / 2,
    k.all = x$k,
    kdz = x$kdz,
    ksz = x$ksz,
    k1sz = x$k1sz,
    k2sz = x$k2sz,
    ids = 1:length(ai),
    # effect sizes and sampling variances:
    # yi = yi,
    # vi = vi,
    # model matrix:
    X = X,
    # counts and sample sizes::
    ai = ai,
    bi = bi,
    ci = ci,
    di = di,
    ni = n1i + n2i,
    n1i = n1i,
    n2i = n2i,
    # arguments:
    measure = measure,
    intercept = intercept,
    slope = slope,
    conditional = conditional,
    approx = approx,
    coding = coding,
    cor = cor,
    drop00 = drop00,
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
