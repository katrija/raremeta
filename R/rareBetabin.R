#' Conduct a meta-analysis using the beta-binomial model
#'
#' Function to conduct a meta-analysis of a rare events using a beta-binomial model.
#' See below for more details on this model and its application in meta-analyses of
#' rare events.
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
#' (either `"logOR"` for the log odds ratio or `"logRR"` for the log relative risk).
#' @param common_rho logical specifying whether a common intraclass correlations shall be assumed
#' for the two groups (`TRUE`), or whether the intraclass correlation shall be allowed
#' to differ between the two groups (`FALSE`). See below for more detail.s
#' @param drop00 logical indicating whether double-zero studies (i.e., studies with no events or
#' only events in both groups) should be excluded when calculating the outcome measufit for the
#' individual studies.
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
#' Data input can happen either through the parameters `ai`,`bi`,`ci`,`di`,`n1i`,`n2i` (columns of the data frame `data`)
#' or pre-processed throuth the parameter `x` (an object of type `rareData`).
#' A `rareData` object can be produced from a data frame by applying the `rareDescribe()` function to it.
#' The `rareDescribe()` function pre-processes the data frame and stores the information required by the `rareBetabin()` function in a list.
#' See `?rareDescribe` for more details.
#'
#' ## Effect size measures
#' The function includes different versions of the beta-binomial model (see Kuss (2014), for a description).
#' The regression equation used in all these models can be expressed as
#'
#' \mjseqn{g(\mu_j) = \alpha + \theta \cdot x_{ij}},
#'
#' where \mjseqn{i = 1, ..., k } is the study index and \mjseqn{j = 1, 2} is the group index,
#' \mjseqn{x_{i1} = 1} and \mjseqn{x_{i2} = 0}.
#'
#' All models assume that the number of events in each group follows a beta-binomial distribution
#' with mean \mjseqn{n_{ij}\mu_j} and intraclass correlation \mjseqn{\rho_j}.
#'
#' For `measure = "logOR"`, \mjseqn{g(\cdot)} corresponds to the logit function, i.e., \mjseqn{g(\mu_j) = \log \left( \frac{\mu_j}{1-\mu_j} \right)}.
#'
#' For `measure = "logRR"`, \mjseqn{g(\cdot)} correspond to the natural logarithm.
#'
#' It is currently not possible to use `rareBetabin` with `measure = "RD"`,
#' but this functionality may be implemented in future versions of this package.
#'
#' ## Group-specific intraclass correlation
#' Per default, it is assumed that \mjseqn{\rho_1 = \rho_2}, i.e. the intraclass correlations
#' are assumed to be equal for both groups. It is possible to fit the beta-binomial
#' with different intraclass correlations by setting `common_rho` to `FALSE`.
#' Internally, different intraclass correlations are modeled via the regression equation
#'
#' \mjseqn{\log\left(\frac{\rho_j}{1-\rho_j}\right) = \zeta + \gamma \cdot x_{ij}},
#'
#' where \mjseqn{x_{i1} = 1} and \mjseqn{x_{i2} = 0}. Estimates
#' for \mjseqn{\rho_1} and \mjseqn{\rho_2} are then obtained by back-transforming
#' the results to the original scale.
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
#' * `rho`: estimated intraclass correlation. If `common_rho = FALSE`, rho is a vector
#' which contains both intraclass correlations.
#' * `LRT.Chisq`: Test statistic of the likelihood ratio test testing for homogeneity.
#' * `LRT.df`: Degrees of freedom of the likelihood ratio test testing for homogeneity.
#' * `LRT.pval`: p-value of the likelihood ratio test testing for homogeneity.
#' * `fit.stats`: a list with log-likelihood, deviance, AIC, BIC, and AICc values under
#' the unrestricted and restricted likelihood.
#' * `p`: number of coefficients in the model (including the intercept).
#' * `k`: number of studies included in the analysis.
#' * `k.all`: total number of studies (before exclusion).
#' * `kdz`,`ksz`: number of double-zero and single-zero studies.
#' * `k1sz`, `k2sz`: number of single-zero studies where the zero is in group 1 or group 2.
#' * `ai`, `bi`, `ci`, `di`: original entries of the 2x2 tables for all studies.
#' * `ni`, `n1i`, `n2i`: original total and group sample sizes.
#' * ...
#'
#' @references
#' Kuss, O. (2014). Statistical methods for meta-analyses including information from studies
#' without any events-add nothing to nothing and succeed nevertheless. Statistics in
#' Medicine, 34 (7), 1097–1116. doi: 10.1002/sim.6383
#'
#' @examples
#'
#' # indtroduce the data
#' data <- data.frame(
#' ai = c(0, 3, 2, 0),
#' bi = c(20, 18, 15, 19),
#' ci = c(1, 4, 0, 0),
#' di = c(19, 17, 16, 20)
#' )
#'
#' # estimating the log relative risk assuming common intraclass correlations
#' mRR <- rareBetabin(ai=ai, bi=bi, ci=ci, di=di, data=data, measure="logRR", common_rho=TRUE)
#' mRR
#'
#'
#' # estimating the log odds ratio assuming differing intraclass correlations
#' # data is pre-processed by use of the `rareDescribe()` function
#'
#' x   <- rareDescribe(ai=ai, bi=bi, ci=ci, di=di, data=data)
#' mOR <- rareBetabin(x, measure="logOR", common_rho=FALSE)
#' mOR
#' @export
#'
rareBetabin <- function(x, ai, bi, ci, di, n1i, n2i, data, measure,
                        common_rho = TRUE,
                        drop00 = FALSE,
                        level = 95,
                        test = "z", digits = 4, verbose = FALSE, control,
                        ...) {

  # defining an object of class 'raredata' if raw data is put in
  if(missing(x)){
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

  nobs <- nrow(dataLong)
  nobs1 <- length(ai)
  nobs2 <- length(ci)

  X <- cbind(rep(1, nobs), c(rep(1, nobs1), rep(0, nobs2)))
  p <- 2

  if (common_rho == TRUE) {
    Xr <- cbind(rep(1, nobs))
    pr <- 1
    br.init <- -1
  } else {
    Xr <- cbind(c(rep(0, nobs1), rep(1, nobs2)), c(rep(1, nobs1), rep(0, nobs2)))
    pr <- 2
    br.init <- c(-1, -1)
  }

  parms <- p + pr
  b.init <- c(-2, 0)

  # ML fit ---------------------------------------------------------------------
  params.init <- c(b.init, br.init)

  if(measure == "logOR"){
    fam <- stats::binomial(link = "logit")
    mFE <- cbind(y, n-y) ~ 1 + group
  }

  if(measure == "logRR"){
    fam <- stats::binomial(link = "log")
    mFE <- cbind(y, n-y) ~ 1 + group
  }

  if(measure == "RD"){
    fam <- stats::binomial(link = "identity")
    mFE <- cbind(y, n-y) ~ 1 + group
  }

  invlink <- fam$linkinv
  link <- fam$link

  logit <- stats::binomial(link = "logit")$linkfun
  expit <- stats::binomial(link = "logit")$linkinv

  # log-likelihood ------------------------
  negll <- function(params, yi, mi, X, Xr) {
    b.init <- params[1:ncol(X)]
    br.init <- params[(ncol(X) + 1):(ncol(X) + ncol(Xr))]

    mu <- invlink(X %*% b.init)
    rho <- expit(Xr %*% br.init)

    rho <- ifelse(rho == 0, 0.00001,
      ifelse(rho == 1, 0.99999, rho)
    )
    mu <- ifelse(mu == 0, 0.00001,
      ifelse(mu == 1, 0.99999, mu)
    )

    llcomp <- sapply(1:nrow(X), function(obs) {
      yi_obs <- yi[obs]
      mi_obs <- mi[obs]
      ni_obs <- yi_obs + mi_obs
      mu_obs <- mu[obs]
      rho_obs <- rho[obs]

      alpha_obs <- mu_obs * (1 - rho_obs) / rho_obs
      beta_obs <- (1 - mu_obs) * (1 - rho_obs) / rho_obs

      if (yi_obs > 0) {
        comp1 <- sum(log(alpha_obs + (0:(yi_obs - 1))))
      } else {
        comp1 <- 0
      }

      if (mi_obs > 0) {
        comp2 <- sum(log(beta_obs + (0:(mi_obs - 1))))
      } else {
        comp2 <- 0
      }

      if (ni_obs > 0) {
        comp3 <- sum(log(alpha_obs + beta_obs + (0:(ni_obs - 1))))
      } else {
        comp3 <- 0
      }

      comp <- lchoose(ni_obs, yi_obs) + comp1 + comp2 - comp3

      return(comp)
    })

    negll <- -sum(llcomp)

    return(negll)
  }

  fitML <- try(
    stats::nlminb(
      start = params.init, objective = negll,
      yi = dataLong$y, mi = dataLong$n - dataLong$y,
      X = X, Xr = Xr
    ),
    silent = TRUE
  )

  if(inherits(fitML, "try-error")){
    stop("Unable to fit model.")
  }

  hessian <- numDeriv::hessian(
    func = negll, x = fitML$par,
    yi = dataLong$y, mi = dataLong$n - dataLong$y,
    X = X, Xr = Xr
  )

  llML <- -fitML$objective

  conv <- ifelse(fitML$convergence == 0, TRUE, FALSE)

  b <- fitML$par[1:p]
  br <- fitML$par[(p+1):(p+pr)]

  hessian.singular <- (det(hessian) == 0)

  if(!hessian.singular){
    vb <- solve(hessian)
    ses <- sqrt(diag(vb))
  }else{
    ses <- rep(NA, nrow(hessian))
  }

  se <- ses[1:p]
  ser <- ses[(p+1):(p+pr)]

  if(test == "z"){
    crit = stats::qnorm(1-level/2)
    zval = b/se
    pval = 2 * stats::pnorm(abs(zval), lower.tail = FALSE)
  }

  if(test == "t"){
    crit = stats::qt(1-level/2, nobs-2)
    zval = b/se
    pval = 2 * stats::pt(abs(zval), df = nobs-2, lower.tail = FALSE)
  }


  ci.lb <- c(b-crit*se)
  ci.ub <- c(b+crit*se)

  rho <- as.numeric(expit(br))

  fitFE <- try(
    stats::glm(mFE,
      data = dataLong,
      family = fam
    ),
    silent = TRUE
  )

  LRT.Chisq <- LRT.df <- LRT.pval <- NA

  if(inherits(fitFE, "try-error")){
    warning("Unable to fit fixed-effects model.")
  }else{
    llFE <- stats::logLik(fitFE)

    LRT.Chisq <- as.numeric(2 * (llML - llFE))
    LRT.df <- parms - attributes(llFE)$df
    LRT.pval <- stats::pchisq(LRT.Chisq, df = LRT.df, lower.tail = FALSE)
  }


  # make results list
  # UNDER CONSTRUCTION
  res <- list(
    # model information:
    model = "rareBetabin",
    b = b,
    beta = b,
    se = se,
    zval = zval,
    pval = pval,
    ci.lb = ci.lb,
    ci.ub = ci.ub,
    vb = vb,
    rho = rho,
    # LRT:
    LRT.Chisq = LRT.Chisq,
    LRT.df = LRT.df,
    LRT.pval = LRT.pval,
    # fit.stats = fit.stats,
    p = parms,
    # convergence information:
    conv = conv,
    # singular = singular,
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
    common_rho = common_rho,
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
