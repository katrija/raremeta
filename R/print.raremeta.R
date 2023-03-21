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
    qpval <- ifelse(x$QEp < .001, "< .001", round(x$QEp, 3))

    cat("\nHeterogeneity: \n")
    cat("\n", paste0("Q(df = ", x$k-1, ") = ", round(x$QE, digits), ", p-val ", qpval), "\n")

    if(is.element(x$method, c("CE", "EE", "FE"))){
      cat("\n",
          "I^2 (total heterogeneity / total variability): ", round(x$I2, 2), " %", "\n"
      )
    }else{
      cat("\n",
          "tau^2 (estimated amount of total heterogeneity):", round(x$tau2, digits), "\n",
          "tau (square root of estimated tau^2 value):     ", round(sqrt(x$tau2), digits), "\n",
          "I^2 (total heterogeneity / total variability):  ", round(x$I2, 2), " %", "\n"
      )
    }

    # Model results:
    signif <- stats::symnum(x$pval, corr=FALSE, na=FALSE,
                            cutpoints=c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))

    res.table <- data.frame(round(as.numeric(x$beta), digits), round(x$se, digits),
                            round(x$zval, digits), ifelse(x$pval < 0.001, "< .001", round(x$pval, 3)),
                            round(x$ci.lb, digits), round(x$ci.ub, digits))
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

    res.table <- data.frame(round(as.numeric(x$beta), digits), round(x$se, digits),
                            round(x$zval, digits), ifelse(x$pval < 0.001, "< .001", round(x$pval, 3)),
                            round(x$ci.lb, digits), round(x$ci.ub, digits))
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
          "sigma^2 (estimated variance of the intercepts):        ", round(x$sigma2[[1]][1], digits), "\n",
          "sigma (estimated standard deviation of the intercepts):", round(sqrt(x$sigma2[[1]][1]), digits), "\n"
      )

    }

    if(x$slope == "random"){
      cat("Random-effects meta-analysis using the Generalised Linear Mixed Model (GLMM):", "\n")

      cat("\nNumber of studies:", x$k, "\n")
      cat("\nDouble-zero studies were", drop00, "the analysis.", "\n")

      cat("\nHeterogeneity: \n")

      # Heterogeneity:
      LRTp <- ifelse(x$LRT.pval < .001, "< .001", round(x$LRT.pval, 3))
      cat("\n", paste0("LRT(df = ", x$LRT.df, ") = ", round(x$LRT.Chisq, digits), ", p-val ", LRTp), "\n")


      if(x$intercept == "random"){
        cat("\n",
            "sigma^2 (estimated variance of the intercepts):        ", round(x$sigma2[[1]][1], digits), "\n",
            "sigma (estimated standard deviation of the intercepts):", round(sqrt(x$sigma2[[1]][1]), digits), "\n"
        )
      }

      cat("\n",
          "tau^2 (estimated variance of the effect sizes):        ", round(x$tau2, digits), "\n",
          "tau (estimated standard deviation of the effect sizes):", round(sqrt(x$tau2), digits), "\n")

      if(x$cor){
        cat("\n",
            "Random-effects correlation: ", round(x$sigma2[[1]][2]/sqrt(x$sigma2[[1]][1]*x$sigma2[[1]][4]), digits), "\n")
      }

    }

    # Model results:
    signif <- stats::symnum(x$pval, corr=FALSE, na=FALSE,
                            cutpoints=c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))

    res.table <- data.frame(round(as.numeric(x$beta), digits), round(x$se, digits),
                            round(x$zval, digits), ifelse(x$pval < 0.001, "< .001", round(x$pval, 3)),
                            round(x$ci.lb, digits), round(x$ci.ub, digits))
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

    cat("Random-effects meta-analysis using the Beta-binomial Model:", "\n")

    cat("\nNumber of studies:", x$k, "\n")
    cat("\nDouble-zero studies were", drop00, "the analysis.", "\n")

      cat("\nHeterogeneity: \n")

      # Heterogeneity:
      LRTp <- ifelse(x$LRT.pval < .001, "< .001", round(x$LRT.pval, 3))
      cat("\n", paste0("LRT(df = ", x$LRT.df, ") = ", round(x$LRT.Chisq, digits), ", p-val ", LRTp), "\n")

      if(x$common_rho == TRUE){
        cat("\n",
            "rho (estimated ICC): ", round(x$rho, digits), "\n"
        )
      }else{
        cat("\n",
            "rho (estimated ICC), group 1: ", round(x$rho[1], digits), "\n",
            "rho (estimated ICC), group 2: ", round(x$rho[2], digits), "\n"
        )
      }

    # Model results:
    signif <- stats::symnum(x$pval, corr=FALSE, na=FALSE,
                            cutpoints=c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))

    res.table <- data.frame(round(as.numeric(x$beta), digits), round(x$se, digits),
                            round(x$zval, digits), ifelse(x$pval < 0.001, "< .001", round(x$pval, 3)),
                            round(x$ci.lb, digits), round(x$ci.ub, digits))
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
