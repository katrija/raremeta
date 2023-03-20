#implementing Peto's method for estimating log(OR) in FE context

rarePeto <- function(x, level = 95, digits = 4){


  # argument checking

  # check if x is an object of class rareData
  if(!inherits(x,"rareData")){
    stop("x must be an object of class 'rareData'. See ?rareDescribe for more details.")
  }

  # check if digits argument is valid
  if(length(digits) != 1 || digits%%1 != 0 || digits < 0){
    stop("'digits' must be an integer of length 1.")
  }

  # check if level argument is valid:
  if(!is.numeric(level) || length(level) > 1 || level < 0 || level > 100){
    stop("level must be a scalar between 0 and 100.")
  }

  mf <- match.call()

  # extract counts and sample sizes
  ai  <- x$ai
  bi  <- x$bi
  ci  <- x$ci
  di  <- x$di
  n1i <- x$n1i
  n2i <- x$n2i
  ni  <- n1i + n2i

  # calculating log(OR) estimate via O-E statistic
  Ai    <- ai + ci
  Bi    <- bi + di
  Vi    <- (Ai * Bi * n1i * n2i) / (ni^2 *(ni - 1))
  sumVi <- sum(Vi)

  quant <- stats::qnorm((100-level)/200)

  beta  <- sum((ci-Ai*(n1i/ni)))/sumVi
  se    <- sqrt(1/sumVi) #works only under the premise that the ORs are 1?
  zval  <- beta / se
  pval  <- 2*stats::pnorm(abs(zval), lower.tail=FALSE)
  ci.lb <- beta + quant * se
  ci.ub <- beta - quant * se


  res   <- append(x, list(beta = beta, b = beta, se = se, zval = zval,
                          pval = pval, ci.lb = ci.lb, ci.ub = ci.ub))


  return(res)
}
