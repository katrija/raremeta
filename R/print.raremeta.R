#' Print Method for raremeta objects.
#'
#' Print method for objects of class "raremeta".
#'
#' @param x an object of class raremeta
#' @param digits number of digits to be shown in the print output
#' @param ... other arguments.
#'
#' @return The print functions do not return an object
#' @export
print.raremeta <- function(x, digits, ...){

  if (!is.element("raremeta", class(x))){
    stop("x must be an object of class 'raremeta'.")
  }

  if(missing(digits)){
    digits <- x$digits
  }

  # rareIV ---------------------------------------------------------------------
  if(x$model == "rareIV"){

    cc <- ifelse(x$cc != "tacc", x$cc, "treatment-arm")
    drop00 <- ifelse(x$drop00 == TRUE, "excluded from", "included in")

    if(is.element(x$method, c("FE", "EE", "CE"))){
      cat("Fixed-effects meta-analysis using the inverse variance model:", "\n")
    }else{
      cat("Random-effects meta-analysis using the inverse variance model",  paste0("(tau^2 estimator: ", x$method, "):"), "\n")
    }

    cat("\nNumber of studies:", x$k)
    cat(paste0("\nContinuity correction: ", cc, ", applied to ", sum(x$cc.studies), " studies"))
    cat("\nDouble-zero studies were", drop00, "the analysis.", "\n")

    # Heterogeneity:
    qpval <- ifelse(x$QEp < .0001, "< .0001", format(round(x$QEp, digits), nsmall = digits))

    cat("\nHeterogeneity: \n")
    cat("\n", paste0("Q(df = ", x$k-1, ") = ", format(round(x$QE, digits), nsmall = digits), ", p-val ", qpval), "\n")

    if(is.element(x$method, c("CE", "EE", "FE"))){
      cat("\n",
          "I^2 (total heterogeneity / total variability): ", format(round(x$I2, 2), nsmall = 2), " %", "\n"
      )
    }else{
      cat("\n",
          "tau^2 (estimated amount of total heterogeneity):", format(round(x$tau2, digits), nsmall = digits,
                                                                     scientific = FALSE), "\n",
          "tau (square root of estimated tau^2 value):     ", format(round(sqrt(x$tau2), digits), nsmall = digits,
                                                                     scientific = FALSE), "\n",
          "I^2 (total heterogeneity / total variability):  ", format(round(x$I2, 2), nsmall = 2), " %", "\n"
      )
    }

    # Model results:
    signif <- stats::symnum(x$pval, corr=FALSE, na=FALSE,
                            cutpoints=c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))

    res.table <- format(data.frame(round(as.numeric(x$beta), digits), round(x$se, digits),
                            round(x$zval, digits), ifelse(x$pval < 0.001, "< .001", round(x$pval, digits)),
                            round(x$ci.lb, digits), round(x$ci.ub, digits)), nsmall = digits, scientific = FALSE)
    names(res.table) <- c(as.character(x$measure), "se", "zval", "pval",  "ci.lb", "ci.ub")


    res.table <- cbind(res.table, signif)
    colnames(res.table)[ncol(res.table)] <- ""

    cat("\nModel results: \n")
    cat("\n")
    print(res.table, row.names = FALSE)
    cat("\n")

    cat("---")
    cat("\nSignif. codes: ", "'***': < .001 '**': < .01 '*': < .05 '.': < .1" )
  }

  # rareMH ---------------------------------------------------------------------
  if(x$model == "rareMH"){

    #cc <- ifelse(x$cc != "tacc", x$cc, "treatment-arm")
    #drop00 <- ifelse(x$drop00 == TRUE, "excluded from", "included in")

    cat("Fixed-effects meta-analysis using the Mantel-Haenszel method:", "\n")

    cat("\nNumber of studies:", x$k, "\n")
    #cat(paste0("\nContinuity correction: ", cc, ", applied to ", sum(x$cc.studies), " studies"))
    #cat("\nDouble-zero studies were", drop00, "the analysis.", "\n")

    # Model results:
    signif <- stats::symnum(x$pval, corr=FALSE, na=FALSE,
                            cutpoints=c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))

    res.table <- format(data.frame(round(as.numeric(x$beta), digits), round(x$se, digits),
                            round(x$zval, digits), round(x$pval, digits),
                            round(x$ci.lb, digits), round(x$ci.ub, digits)), nsmall = digits, scientific = FALSE)
    names(res.table) <- c(as.character(x$measure), "se", "zval", "pval",  "ci.lb", "ci.ub")


    res.table <- cbind(res.table, signif)
    colnames(res.table)[ncol(res.table)] <- ""

    cat("\nModel results: \n")
    cat("\n")
    print(res.table, row.names = FALSE)
    cat("\n")

    cat("---")
    cat("\nSignif. codes: ", "'***': < .001 '**': < .01 '*': < .05 '.': < .1" )
  }

  # rarePeto ---------------------------------------------------------------------
  if(x$model == "rarePeto"){

    #cc <- ifelse(x$cc != "tacc", x$cc, "treatment-arm")
    #drop00 <- ifelse(x$drop00 == TRUE, "excluded from", "included in")

    cat("Fixed-effects meta-analysis using Peto's method:", "\n")

    cat("\nNumber of studies:", x$k, "\n")
    #cat(paste0("\nContinuity correction: ", cc, ", applied to ", sum(x$cc.studies), " studies"))
    #cat("\nDouble-zero studies were", drop00, "the analysis.", "\n")

    # Model results:
    signif <- stats::symnum(x$pval, corr=FALSE, na=FALSE,
                            cutpoints=c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))

    res.table <- format(data.frame(round(as.numeric(x$beta), digits), round(x$se, digits),
                            round(x$zval, digits), round(x$pval, digits),
                            round(x$ci.lb, digits), round(x$ci.ub, digits)), nsmall = digits, scientific = FALSE)
    names(res.table) <- c(as.character(x$measure), "se", "zval", "pval",  "ci.lb", "ci.ub")


    res.table <- cbind(res.table, signif)
    colnames(res.table)[ncol(res.table)] <- ""

    cat("\nModel results: \n")
    cat("\n")
    print(res.table, row.names = FALSE)
    cat("\n")

    cat("---")
    cat("\nSignif. codes: ", "'***': < .001 '**': < .01 '*': < .05 '.': < .1" )
  }

  # rareGLMM -------------------------------------------------------------------
  if(x$model == "rareGLMM"){

    drop00 <- ifelse(x$drop00 == TRUE, "excluded from", "included in")

    if(x$conditional == FALSE){
      if(x$slope == "fixed" & x$intercept == "fixed"){
        cat("Fixed-effects meta-analysis using the Generalised Linear Model (GLM):", "\n")

        cat("\nNumber of studies:", x$k, "\n")
        cat("\nDouble-zero studies were", drop00, "the analysis.", "\n")
      }

      if(x$slope == "fixed" & x$intercept == "random"){
        cat("Fixed-effects meta-analysis using the Generalised Linear Mixed Model (GLMM):", "\n")

        cat("\nNumber of studies:", x$k, "\n")
        cat("\nDouble-zero studies were", drop00, "the analysis.", "\n")

        cat("\nHeterogeneity: \n")

        cat("\n",
            "sigma^2 (estimated variance of the intercepts):        ", format(round(x$sigma2[[1]][1], digits), nsmall = digits,
                                                                              scientific = FALSE), "\n",
            "sigma (estimated standard deviation of the intercepts):", format(round(sqrt(x$sigma2[[1]][1]), digits), nsmall = digits,
                                                                              scientific = FALSE), "\n"
        )

      }

      if(x$slope == "random"){
        cat("Random-effects meta-analysis using the Generalised Linear Mixed Model (GLMM):", "\n")

        cat("\nNumber of studies:", x$k, "\n")
        cat("\nDouble-zero studies were", drop00, "the analysis.", "\n")

        cat("\nHeterogeneity: \n")

        # Heterogeneity:
        LRTp <- ifelse(x$LRT.pval < .0001, "< .0001", format(round(x$LRT.pval, digits), nsmall = digits, scientitific = FALSE))
        cat("\n", paste0("LRT(df = ", x$LRT.df, ") = ", format(round(x$LRT.Chisq, digits), nsmall = digits, scientific = FALSE), ", p-val ", LRTp), "\n")


        if(x$intercept == "random"){
          cat("\n",
              "sigma^2 (estimated variance of the intercepts):        ", format(round(x$sigma2[[1]][1], digits), nsmall = digits, scientific = FALSE), "\n",
              "sigma (estimated standard deviation of the intercepts):", format(round(sqrt(x$sigma2[[1]][1]), digits), nsmall = digits, scientific = FALSE), "\n"
          )
        }

        cat("\n",
            "tau^2 (estimated variance of the effect sizes):        ", format(round(x$tau2, digits), nsmall = digits, scientific = FALSE), "\n",
            "tau (estimated standard deviation of the effect sizes):", format(round(sqrt(x$tau2), digits), nsmall = digits, scientific = FALSE), "\n")

        if(x$cor){
          cat("\n",
              "Random-effects correlation: ", format(round(x$sigma2[[1]][2]/sqrt(x$sigma2[[1]][1]*x$sigma2[[1]][4]), digits), nsmalL = digits, scientific = FALSE), "\n")
        }

      }
    }

    if(x$conditional == TRUE){

      if(x$intercept == "fixed"){
        cat("Fixed-effects meta-analysis using the Conditional Generalised Linear Mixed Model (GLMM):", "\n")

        cat("\nNumber of studies:", x$k, "\n")
        cat("\nDouble-zero studies were", drop00, "the analysis.", "\n")
      }

      if(x$intercept == "random"){
        cat("Random-effects meta-analysis using the Conditional Generalised Linear Mixed Model (GLMM):", "\n")

        cat("\nNumber of studies:", x$k, "\n")
        cat("\nDouble-zero studies were", drop00, "the analysis.", "\n")

        cat("\nHeterogeneity: \n")

        # Heterogeneity:
        LRTp <- ifelse(x$LRT.pval < .0001, "< .0001", format(round(x$LRT.pval, digits), nsmall = digits, scientitific = FALSE))
        cat("\n", paste0("LRT(df = ", x$LRT.df, ") = ", format(round(x$LRT.Chisq, digits), nsmall = digits, scientific = FALSE), ", p-val ", LRTp), "\n")

        cat("\n",
            "tau^2 (estimated variance of the effect sizes):        ", format(round(x$tau2, digits), nsmall = digits, scientific = FALSE), "\n",
            "tau (estimated standard deviation of the effect sizes):", format(round(sqrt(x$tau2), digits), nsmall = digits, scientific = FALSE), "\n")

      }


    }


    # Model results:
    signif <- stats::symnum(x$pval, corr=FALSE, na=FALSE,
                            cutpoints=c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))

    res.table <- format(data.frame(round(as.numeric(x$beta), digits), round(x$se, digits),
                            round(x$zval, digits), round(x$pval, digits),
                            round(x$ci.lb, digits), round(x$ci.ub, digits)), nsmall = digits, scientific = FALSE)
    names(res.table) <- c("estimate", "se", "zval", "pval",  "ci.lb", "ci.ub")
    rownames(res.table) <- names(x$beta)

    groupInd <- which(rownames(res.table) == "group")
    rownames(res.table)[groupInd] <- as.character(x$measure)

    res.table <- cbind(res.table, signif)
    colnames(res.table)[ncol(res.table)] <- ""

    cat("\nModel results: \n")
    cat("\n")
    print(res.table, row.names = TRUE)
    cat("\n")

    cat("---")
    cat("\nSignif. codes: ", "'***': < .001 '**': < .01 '*': < .05 '.': < .1" )

  }

  # rareBetabin -------------------------------------------------------------
  if(x$model == "rareBetabin"){

    drop00 <- ifelse(x$drop00 == TRUE, "excluded from", "included in")

    cat("Random-effects meta-analysis using the beta-binomial model:", "\n")

    cat("\nNumber of studies:", x$k, "\n")
    cat("\nDouble-zero studies were", drop00, "the analysis.", "\n")

      cat("\nHeterogeneity: \n")

      # Heterogeneity:
      # LRTp <- ifelse(x$LRT.pval < .0001, "< .0001", format(round(x$LRT.pval, digits), nsmall = digits, scientific = FALSE))
      # cat("\n", paste0("LRT(df = ", x$LRT.df, ") = ", format(round(x$LRT.Chisq, digits), nsmall = digits, scientific = FALSE), ", p-val ", LRTp), "\n")

      if(x$common_rho == TRUE){
        cat("\n",
            "rho (estimated ICC): ", format(round(x$rho, digits), nsmall = digits, scientific = FALSE), "\n"
        )
      }else{
        cat("\n",
            "rho (estimated ICC), group 1: ", format(round(x$rho[1], digits), nsmall = digits, scientific = FALSE), "\n",
            "rho (estimated ICC), group 2: ", format(round(x$rho[2], digits), nsmall = digits, scientific = FALSE), "\n"
        )
      }

    # Model results:
    signif <- stats::symnum(x$pval, corr=FALSE, na=FALSE,
                            cutpoints=c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))

    res.table <- format(data.frame(round(as.numeric(x$beta), digits), round(x$se, digits),
                            round(x$zval, digits), round(x$pval, digits),
                            round(x$ci.lb, digits), round(x$ci.ub, digits)), nsmall = digits, scientific = FALSE)
    names(res.table) <- c("estimate", "se", "zval", "pval",  "ci.lb", "ci.ub")
    rownames(res.table) <- c("Intercept", x$measure)

    groupInd <- which(rownames(res.table) == "group")
    rownames(res.table)[groupInd] <- as.character(x$measure)

    res.table <- cbind(res.table, signif)
    colnames(res.table)[ncol(res.table)] <- ""

    cat("\nModel results: \n")
    cat("\n")
    print(res.table, row.names = TRUE)
    cat("\n")

    cat("---")
    cat("\nSignif. codes: ", "'***': < .001 '**': < .01 '*': < .05 '.': < .1" )

  }



}
