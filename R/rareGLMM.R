#' Conduct a meta-analysis using a generalized linear mixed model (GLMM)
#'
#' @param x
#' @param measure
#' @param intercept
#' @param slope
#' @param conditional
#' @param drop00
#' @param level
#' @param test
#' @param digits
#' @param verbose
#' @param control
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
rareGLMM <- function(x, measure, intercept = "fixed", slope = "random", conditional = FALSE,
                     drop00 = FALSE,
                     level = 95,
                     test="z", digits = 4, verbose=FALSE, control,
                     ...){

  # check if x is an object of class rareData
  if(!inherits(x,"rareData")){
    stop("x must be an object of class 'rareData'. See ?rareDescribe for more details.")
  }

  # check if measure argument is specified
  if(missing(measure)){
    stop("'measure' argument must be specified.")
  }

  # check if measure argument is valid
  if(!is.element(measure, c("logOR", "logRR"))){
    stop("'measure' must be either 'logOR' or 'logRR'.")
  }

  # extract counts and sample sizes
  ai <- x$ai
  bi <- x$bi
  ci <- x$ci
  di <- x$di
  n1i <- x$n1i
  n2i <- x$n2i

  if(!is.logical(drop00)){
    stop("'drop00' must be a logical")
  }

  # check if level argument is valid:
  if(!is.numeric(level) || length(level) > 1 || level < 0 || level > 100){
    stop("level must be a scalar between 0 and 100.")
  }

  # check if test argument is valid
  if(!is.element(test, c("z"))){
    stop("'test' must be 'z'.")
  }

  # check if digits argument is valid
  if(length(digits) != 1 || digits%%1 != 0 || digits < 0){
    stop("'digits' must be an integer of length 1.")
  }

  mf <- match.call()

  level <- (100-level)/100

  # remove double-zero studies if desired:
  if(drop00 == TRUE){
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
    group = c(rep(1,length(ai)), rep(0, length(ci)))
  )

  dataLong$group2 = ifelse(dataLong$group == 1, 0.5, -0.5)


  if(measure == "logOR"){

    if(intercept == "random" & slope == "random" & conditional == FALSE){
      fit <- lme4::glmer(cbind(y,n-y)~1+group+(1+group||id), data = dataLong,
                         family = binomial(link = "logit"))
    }

    if(intercept == "random" & slope == "random" & conditional == FALSE){
      fit <- lme4::glmer(cbind(y,n-y)~1+id+group+(0+group|id), data = dataLong,
                         family = binomial(link = "logit"))
    }


  }

  if(measure == "logRR"){

    if(intercept == "random" & slope == "random" & conditional == FALSE){
      fit <- lme4::glmer(y~1+group+offset(log(n))+(1+group||id), data = dataLong,
                         family = poisson(link = "log"))
    }

    if(intercept == "random" & slope == "random" & conditional == FALSE){
      fit <- lme4::glmer(y~1+id+group+offset(log(n))+(0+group|id), data = dataLong,
                         family = poisson(link = "log"))
    }

  }

  beta <- lme4::fixef(fit)
  vb <- as.matrix(vcov(fit))
  sigma2 <- lme4::VarCorr(fit)
  tau2 <- sigma2[[length(sigma2)]][1]

  se <- sqrt(diag(vb))

  names(se) <- NULL
  names(beta) <- NULL

  zval <- beta/se
  pval <- 2*stats::pnorm(abs(zval), lower.tail=FALSE)

  zcrit <- stats::qnorm(1-level/2)
  ci.lb <- beta - zcrit*se
  ci.ub <- beta + zcrit*se

  X <- stats::model.matrix(fit)
  p <- ifelse(all(X[1,]==1), ncol(X)-1, ncol(X))

  AICc <- -2*logLik(fit)+2*(p+1)*max(2*nrow(dataLong), p+3)/(max(nrow(dataLong), 3)-p)
  fit.stats <- rbind(stats::logLik(fit), stats::deviance(fit), stats::AIC(fit), stats::BIC(fit), AICc)
  rownames(fit.stats) <- c("ll", "dev", "AIC", "BIC", "AICc")
  colnames(fit.stats) <- c("ML")

  singular <- lme4::isSingular(fit)
  conv <- fit@optinfo$conv$opt

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
    # se.tau2 = fit$se.tau2,
    #I2 = fit$I2,
    #H2 = fit$H2,
    #R2 = fit$R2,
    #vt = fit$vt,
    #QE = fit$QE,
    #QEp = fit$QEp,
    #QM = fit$QM,
    #QMdf = fit$QMdf,
    #QMp = fit$QMp,
    fit.stats = fit.stats,
    p = p,
    # convergence information:
    conv = conv,
    singular = singular,
    # study numbers:
    k = nrow(dataLong),
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
    ni = n1i+n2i,
    n1i = n1i,
    n2i = n2i,
    # arguments:
    measure = measure,
    intercept = intercept,
    slope = slope,
    conditional = conditional,
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
