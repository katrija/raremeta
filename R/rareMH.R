##Mantel-Haenszel estimator for OR in the FE model
#
#Input is data of type rareDate obtained by use of the rareDescribe() function
#Confidence intervals and significance tests are conducted using
#Wald-type tests and CIs [true?]

rareMH <- function(x, measure, level = 95,  digits = 4){

  # check if x is an object of class rareData
  if(!class(x) == "rareData"){
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

  # check if digits argument is valid
  if(length(digits) != 1 | digits%%1 != 0 | digits < 0){
    stop("'digits' must be an integer of length 1.")
  }

  # check if level argument is valid:
  if(!is.numeric(level) | length(level) > 1 | level < 0 | level > 100){
    stop("level must be a scalar between 0 and 100.")
  }

  # extract counts and sample sizes
  ai  <- x$ai
  bi  <- x$bi
  ci  <- x$ci
  di  <- x$di
  n1i <- x$n1i
  n2i <- x$n2i
  ni  <- n1i + n2i

  #calculating effect sizes
  if(measure == "logOR"){
    Ai <- (ai+di)/ni
    Bi <- (ai*di)/ni
    Ci <- (bi+ci)/ni
    Di <- (bi*ci)/ni

    B <- sum(Bi)
    D <- sum(Di)

    if(B == 0 || D == 0){
      stop("The data does not allow for the application of this method.
           See e.g. ?rareDescribe() for possible continuity corrections.")
    }

    else{
    beta   <- log(B/D)
    se     <- sqrt(1/2 * (sum(Ai*Bi)/B^2 + sum(Ai*Di + Ci*Bi)/(B*D)
                    + sum(Ci*Di)/D^2))
    zval   <- beta / se
    pval   <- 2*pnorm(abs(zval), lower.tail=FALSE)
    ci.lb  <- beta + qnorm((100-level)/200) * se
    ci.ub  <- beta - qnorm((100-level)/200) * se
    }
  }

  if(measure == "logRR"){

    A <- sum(ai*n2i/ni)
    B <- sum(ci*n1i/ni)

    if(A == 0 || B == 0){
      stop("The data does not allow for application of this method.
           See e.g. ?rareDescribe() for possible continuity corrections.")
    }

    else{
      beta   <- log(A/B)
      se     <- sqrt(sum(n1i*n2i*(ai+ci)/ni^2 - ai*ci/ni)/(A*B))
      zval   <- beta / se
      pval   <- 2*pnorm(abs(zval), lower.tail=FALSE)
      ci.lb  <- beta + qnorm((100-level)/200) * se
      ci.ub  <- beta - qnorm((100-level)/200) * se
    }
  }

  if(measure == "RD"){

    beta   <- sum(ai*(n2i/ni) - ci*(n1i/ni))/sum(n1i*(n2i/ni))
    se     <- sqrt(sum((ai*bi*n2i^3 + ci*di*n1i^3)/(n1i*n2i*ni^2))/
                       (sum(n1i*n2i/ni)^2))
    zval   <- beta / se
    pval   <- 2*pnorm(abs(zval), lower.tail=FALSE)
    ci.lb  <- beta + qnorm((100-level)/200) * se
    ci.ub  <- beta - qnorm((100-level)/200) * se
  }

  res <- list(
    # model information
    model = "rareMH",
    b = beta,
    beta = beta,
    se = se,
    vb = NA, #change print method to allow for dropping vb
    zval = zval,
    pval = pval,
    ci.lb = ci.lb,
    ci.ub = ci.ub
    #b.exp = exp(beta),
    #beta.exp = exp(beta),
    #ci.lb.exp = exp(ci.lb),
    #ci.ub.exp = exp(ci.ub)
  )

  #res <- raremeta(res)
  return(res)
}



##QUESTIONS

#(1)
#Do we incorporate different tests than Wald-type tests?
#Do we even call it Wald-type testing if we do not consider ML-estimates
#(The estimator is normal anyway but out of different(?) reasons)

#code we might need if we decide to incorporate different types of tests:
##check if test argument is valid
#if(!is.element(test, c("z", "knha", "hksj"))){
#  stop("'test' must be either 'z', 'knha', or 'hksj'")
#}

## convert test to metafor notation if needed
#if(test == "hksj"){
#  test = "knha"
#}

#(2)
#Parts of the rareIV code maybe needed here as well

# convert measure to the metafor notation
#metafor_measure <- sub("log", "", measure)
#mf <- match.call() #What does this do?

#time.start <- proc.time()
#time.end <- proc.time()
#res$time <- unname(time.end - time.start)[3]

#(3)
#references to the original papers

#(4)
#How to structure the output? Incorporate our 'print()'-method 'raremeta'?
#Add rounding by 'digits' argument




