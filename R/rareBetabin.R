rareBetabin <- function(x, measure,
                        common_rho = TRUE,
                        drop00 = FALSE,
                        level = 95,
                        test = "z", digits = 4, verbose = FALSE, control,
                        ...) {
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
    mFE <- y ~ 1 + group
  }

  if(measure == "RD"){
    fam <- stats::binomial(link = "identity")
    mFE <- y ~ 1 + group
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
