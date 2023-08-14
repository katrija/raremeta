#' Print Method for rareData objects.
#'
#' Print method for objects of class "rareData".
#'
#' @param x an object of class rareData
#' @param ... other arguments.
#'
#' @return The print functions do not return an object
#' @export
print.rareData <- function(x, ...) {
  sampleSize <- round(as.matrix(rbind(x$n1, x$n2, x$n)), 0)
  sampleSize <- sampleSize[, c("mean", "sd", "min", "q25", "median", "q75", "max", "IQR")]
  row.names(sampleSize) <- c("group 1", "group 2", "total")

  sampleRatio <- round(x$nratio, 2)
  sampleRatio <- t(sampleRatio[c("mean", "sd", "min", "q25", "median", "q75", "max", "IQR")])
  row.names(sampleRatio) <- "group1:group2"

  relFreq <- round(as.matrix(rbind(x$rf1, x$rf2, x$rf)), 4)
  relFreq <- relFreq[, c("mean", "sd", "min", "q25", "median", "q75", "max", "IQR")]
  row.names(relFreq) <- c("group 1", "group 2", "total")

  cat("\nStudies: \n",
    "Total number of studies:       ", x$k, "\n",
    "Number of double-zero studies: ", x$kdz, "\n",
    "Number of single-zero studies: ", x$ksz, " (", x$k1sz, " zeros in group 1 and ", x$k2sz, " zeros in group 2) \n",
    sep = ""
  )
  cat("\n")

  cat("\nSample sizes: \n")
  stats::printCoefmat(sampleSize)
  cat("\n")

  cat("\nSample size ratios: \n")
  stats::printCoefmat(sampleRatio)
  cat("\n")

  cat("\nRelative frequencies of the event: \n")
  stats::printCoefmat(relFreq)
  cat("\n")

  # cat("\nNumber of studies in which the event is rare (0.01 < rel. freq. < 0.05):", x$krare, "out of", x$k, "studies.")
  # cat("\nNumber of studies in which the event is very rare (rel. freq. < 0.01):  ", x$kveryrare, "out of", x$k, "studies.")



  #print summary of continuity corrected data
  #(commented out for now)

  #if(all(c("ai.cc","bi.cc","ci.cc","di.cc") %in% names(x))){
  if(FALSE){
    cat("\n")

    cat("\n####################################\n")

    cat("\n")

    cat("\nSummary of continuity corrected data \n" )

    sampleSize <- round(as.matrix(rbind(x$n1.cc, x$n2.cc, x$n.cc)), 2)
    sampleSize <- sampleSize[, c("mean", "min", "q25", "median", "q75", "max")]
    row.names(sampleSize) <- c("group 1", "group 2", "total")

    sampleRatio <- round(x$nratio, 2)
    sampleRatio <- t(sampleRatio[c("mean", "min", "q25", "median", "q75", "max")])
    row.names(sampleRatio) <- "group1:group2"

    relFreq <- round(as.matrix(rbind(x$rf1.cc, x$rf2.cc, x$rf.cc)), 4)
    relFreq <- relFreq[, c("mean", "min", "median", "max")]
    row.names(relFreq) <- c("group 1", "group 2", "total")

    cat("\nStudies: \n",
        "Total number of studies:       ", x$k.cc, "\n",
        "Number of double-zero studies: ", x$kdz.cc, "\n",
        "Number of single-zero studies: ", x$ksz.cc, ", thereof ", x$k1sz.cc, " zeros in group 1 and ", x$k2sz.cc, " zeros in group 2.\n",
        sep = ""
    )
    cat("\n")

    cat("\nSample sizes: \n")
    stats::printCoefmat(sampleSize)
    cat("\n")

    cat("\nSample size ratios: \n")
    stats::printCoefmat(sampleRatio)
    cat("\n")

    cat("\nRelative frequencies of the event: \n")
    stats::printCoefmat(relFreq)
    cat("\n")

    cat("\nNumber of studies in which the event is rare (0.01 < rel. freq. < 0.05):", x$krare.cc, "out of", x$k.cc, "studies.")
    cat("\nNumber of studies in which the event is very rare (rel. freq. < 0.01):  ", x$kveryrare.cc, "out of", x$k.cc, "studies.")

  }

  #print estimated effect sizes
  #(commented out for now)

  #if(all(c("yi", "vi") %in% names(x))){
  if(FALSE){
    cat("\n")

    cat("\n Estimated effect sizes \n")

    cat("\n")

    print(x$measure)
    print(x$yi)

    cat("\n")

    cat("\n variance \n")
    print(x$vi)
  }
}

# print summary of effect size and variance estimates
if(FALSE){
  if(all(c("yi","vi") %in% names(x))){
    EffectSize <- round(as.matrix(rbind(x$yi, x$vi)), 2)
    EffectSize <- EffectSize[, c("mean", "sd", "min", "q25", "median", "q75", "max", "IQR")]
    row.names(EffectSize) <- c(x$measure, "vi")

    cat("\nEffect Size Estimation: \n")
    stats::printCoefmat(EffectSize)
    cat("\n")

  }
}
