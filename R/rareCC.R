# Function for Continuity correction (and only continuity correction)


rareCC <- function(x, cc = "constant", ccval = 0.5, tccval, cccval, ccsum = 1,
                   ccto = "only0", drop00 = TRUE, measure, method = "FE"){


  # check if x is an object of class rareData
  if(!inherits(x,"rareData")){
    stop("x must be an object of class 'rareData'. See ?rareDescribe for more details.")
  }

  # extract counts and sample sizes
  ai <- x$ai
  bi <- x$bi
  ci <- x$ci
  di <- x$di
  n1i <- x$n1i
  n2i <- x$n2i


  ## check if cc is specified
  #if(missing(cc)){
  #  stop("Some studies have zero events. \n
  #       You must specify the 'cc' argument to determine how they are handled.\n
  #       In case you want to exclude all zero-studies, set 'cc' equal to 'none'.")
  #}

  # check if cc argument is valid
  if(!is.element(cc, c("none", "constant", "tacc", "empirical"))){
    stop("'cc' must be either 'none', 'constant', 'tacc', or 'empirical'.")
  }

  # check if drop00 argument is valid
  if(!is.logical(drop00)){
    stop("'drop00' must be a logical.")
  }

  # check if ccto argument is valid (if cc shall be applied)
  if(cc != "none" &&  !is.element(ccto, c("only0", "all", "if0all"))){
    stop("'ccto' must be either 'only0', 'all', or 'if0all'.")
  }

  # check ccval, tccval, cccval arguments (if cc == "constant")
  if(cc == "constant"){

    if(drop00 == FALSE && !is.element(length(ccval), c(1,x$k))){
      stop("'ccval' must have length 1 or length equal to the number of studies.")
    }

    if(drop00 == TRUE && !is.element(length(ccval), c(1,x$k,x$k-x$kdz))){
      stop("'ccval' must have length 1 or length equal to the number of studies (in- or excluding double-zero studies).")
    }

    if(any(ccval < 0)){
      stop("'ccval' must be non-negative.")
    }

    if(!missing(tccval) || !missing(cccval)){

      if((!missing(tccval) && missing(cccval)) || (missing(tccval) && !missing(cccval))){
        stop("Please specify both 'tccval' and 'cccval'.")
      }

      if(!equalLength(tccval, cccval)){
        stop("'tccval' and 'cccval' must have equal length.")
      }

      if(drop00 == FALSE && !is.element(length(tccval), c(1,x$k))){
        stop("'tccval' must have length 1 or length equal to the number of studies.")
      }

      #if(drop00 == TRUE && !is.element(length(tccval), c(1,x$k,x$k-x$kdz))){
      #  stop("'tccval' must have length 1 or length equal to the number of studies (in- or excluding double-zero studies).")
      #}

      if(any(c(tccval, cccval) < 0)){
        stop("All values in 'tccval' and 'cccval' must be non-negative.")
      }
    }
  }


  # check if measure argument is specified (if needed)
  if(cc == "empirical" && missing(measure)){
    stop("To apply the the empirical continuity correction, the 'measure' argument must be specified.")
  }

  # check if measure argument is valid (if needed)
  if(cc == "empirical" && !is.element(measure, c("logOR", "logRR", "RD"))){
    stop("To apply the empirical continuity correction, 'measure' must be either 'logOR', 'logRR' or 'RD'.")
  }

  # check if method argument is valid (if needed)
  if(cc == "empirical"
     && !is.element(method, c("FE", "EE", "DL", "HE", "HS", "HSk", "SJ", "ML",
                               "REML", "EB", "PM", "GENQ", "PMM", "GENQM"))){
    stop("'method' argument must be one of the following: 'FE','EE','DL','HE','HS',
         'HSK','SJ','ML','REML','EB','PM','GENQ','PMM','GENQM'.")
  }

  # tacc currently not supported for RD
  if(cc == "tacc" && measure == "RD"){
    stop("continuity correction of type 'tacc' is currently not supported for measure 'RD'.")
  }

  # empirical cc currently not supported for RD
  if(cc == "empirical" && measure == "RD"){
    stop("continuity correction of type 'empirical' is currently not supported for measure 'RD'.")
  }

  # check if there are non-zero studies when the empirical cc shall be applied
  if(cc == "empirical" && all(ai == 0 | bi == 0 | ci == 0 | di == 0)){
    stop("Continuity correction of type 'empirical' can only be applied if there is at least one non-zero study.")
  }

  # check ccsum argument:
  if(!is.numeric(ccsum) || length(ccsum) > 1 || ccsum < 0){
    stop("ccsum must be a scalar larger than 0.")
  }

  # check whether tccval and cccval are specified; if not, set them to ccval:
  if(missing(tccval)){
    tccval <- ccval
    cccval <- ccval
  }

  remove <- rep(FALSE, length(ai))

  # remove double-zero studies if desired:
  if(drop00 == TRUE){
    remove <- (ai == 0 & ci == 0) | (bi == 0 & di == 0)
    ai <- ai[!remove]
    bi <- bi[!remove]
    ci <- ci[!remove]
    di <- di[!remove]
    n1i <- n1i[!remove]
    n2i <- n2i[!remove]

    if(length(tccval) > length(ai)){
      tccval <- tccval[!remove]
      cccval <- cccval[!remove]
    }
  }

  ccstudies <- rep(FALSE, length(ai))

  # specify studies to be continuity corrected:
  if(ccto == "only0"){
    ccstudies <- (ai == 0 | ci == 0 | bi == 0 | di == 0)
  }

  if(ccto == "all" || (ccto == "if0all" && any(ai == 0 | ci == 0 | bi == 0 | di == 0) )){
    ccstudies <- rep(TRUE, length(ai))
  }

  tcc <- rep(0,length(ai))
  ccc <- rep(0,length(ci))

  if(cc == "constant"){

    if(length(tccval) == 1){
      tcc[ccstudies] <- tccval
      ccc[ccstudies] <- cccval
    }else{
      tcc <- tccval
      ccc <- cccval
    }


  }


  # continuity corrections for logOR and logRR:
  if(cc == "tacc" && is.element(measure, c("logOR", "logRR"))){
    rinv <- n2i[ccstudies]/n1i[ccstudies]

    tcc[ccstudies] <- ccsum*1/(rinv+1)
    ccc[ccstudies] <- ccsum*rinv/(rinv+1)
  }

  # convert measure to the metafor notation
  if(cc == "empirical")
    metafor_measure <- sub("log", "", measure)

  if(cc == "empirical" && is.element(measure, c("logOR", "logRR"))){
    rinv <- n2i[ccstudies]/n1i[ccstudies]
    nozero <- which((ai != 0) & (bi != 0) & (ci != 0) & (di != 0))
    fit_nozero <- metafor::rma(ai = ai[nozero], bi = bi[nozero],
                               ci = ci[nozero], di = di[nozero],
                               measure = metafor_measure, method = method)
    prior <- exp(as.numeric(fit_nozero$beta))

    tcc[ccstudies] <- ccsum*prior/(rinv+prior)
    ccc[ccstudies] <- ccsum*rinv/(rinv+prior)

  }

  ai.cc <- ai; bi.cc <- bi; ci.cc <- ci; di.cc <- di

  # apply continuity correction:
  ai.cc <- ai+tcc
  bi.cc <- bi+tcc
  ci.cc <- ci+ccc
  di.cc <- di+ccc

  #adding description of continuity corrected data
  data.cc <- data.frame(ai.cc = ai.cc, bi.cc = bi.cc, ci.cc = ci.cc, di.cc = di.cc)
  x.cc <- rareDescribe(ai.cc, bi.cc, ci.cc, di.cc, data = data.cc)

  out <- append(list(ai.cc = x.cc$ai, bi.cc = x.cc$bi, ci.cc = x.cc$ci,
                     di.cc = x.cc$di, n1i.cc = x.cc$n1i, n2i.cc = x.cc$n2i,
                     ni.cc = x.cc$ni, nratioi.cc = x.cc$nratioi,
                     nratio.cc = x.cc$nratio, k.cc = x.cc$k, kdz.cc = x.cc$kdz,
                     ksz.cc = x.cc$ksz, k1sz.cc = x.cc$k1sz, k2sz.cc = x.cc$k2sz,
                     n1.cc = x.cc$n1, n2.cc = x.cc$n2, n.cc = x.cc$n,
                     rf1.cc = x.cc$rf1, rf2.cc = x.cc$rf2, rf.cc = x.cc$rf,
                     krare.cc = x.cc$krare, kveryrare.cc = x.cc$kveryrare,
                     k1rare.cc = x.cc$k1rare, k1veryrare.cc = x.cc$k1veryrare,
                     k2rare.cc = x.cc$k2rare, k2veryrare.cc = x.cc$k2veryrare,
                     cc = cc, ccto = ccto, drop00 = drop00, ccstudies = ccstudies,
                     ccc = ccc, tcc = tcc, remove = remove), x)

  #report measure and method argument used in 'tacc'- and 'empirical'- continuity correction
  if(cc == "tacc" || cc == "empirical"){
    out <- append(out, list(method.cc = method, measure.cc = measure))
  }

  out <- rareData(out)
  return(out)

}



