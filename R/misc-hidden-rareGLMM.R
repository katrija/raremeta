#' Non-central hypergeometric likelihood for a single study
#'
#' @noRd
.llnchgi <- function(thetai, ai, bi, ci, di, theta, tau2,
                     intercept = "fixed"){
  mi <- ai+ci
  n1i <- ai+bi
  n2i <- ci+di

  # truncate values which are too large/small:
  pow <- 12

  thetai[thetai < log(10^-pow)] <- log(10^-pow)
  thetai[thetai > log(10^pow)]  <- log(10^pow)

  len <- length(thetai)

  lli <- sapply(1:len, function(i) try(MCMCpack::dnoncenhypergeom(x = ai, n1 = n1i, n2 = n2i, m1 = mi, psi = exp(thetai[i]))))

  if(intercept == "random"){
    lli <- lli*dnorm(thetai, mean = theta, sd = sqrt(tau2))
  }

  return(lli)
}

#' Negative log-likelihood for the hypergeometric-normal model
#'
#' @noRd
.negllnchg <- function(parms, ai, bi, ci, di, intercept = "fixed"){

  theta <- parms[1]

  k <- length(ai)
  thetai <- rep(parms[1], k)

  # likelihood logarithmieren!!!
  if(intercept == "fixed"){
    lli <- sapply(1:k, function(i){log(.llnchgi(thetai = thetai[i], ai = ai[i], bi = bi[i], ci = ci[i], di = di[i],intercept = "fixed"))})
  }

  if(intercept == "random"){
    tau2 <- parms[2]

    if(tau2 <= 0){
      tau2 <- 0.00001
    }

    lli <- sapply(1:k, function(i){

      inti <- try(integrate(.llnchgi, lower = -Inf, upper = Inf,
                            ai = ai[i], bi = bi[i], ci = ci[i], di = di[i],
                            theta = thetai[i], tau2 = tau2,
                            intercept = "random"))

      if(inherits(inti, "try-error")){
        paste("Could not obtain integral for study", i, ".")
      }else{
        if(inti$value > 0){
          out <- log(inti$value)
        }else{
          out <- -Inf
        }
      }

      return(out)
    })

  }

  res <- -sum(lli)

  return(res)

}
