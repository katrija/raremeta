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


}
