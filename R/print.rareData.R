#' Print Method for rareData objects.
#'
#' Print method for objects of class "rareData".
#'
#' @param x an object of class rareData
#' @param ... other arguments.
#'
#' @return The print functions do not return an object
#' @export
print.rareData <- function(x, ...){

  sampleSize <- as.matrix(rbind(x$n1, x$n2,x$n))
  sampleSize <- sampleSize[,c("mean", "min", "q25", "median", "q75", "max")]
  row.names(sampleSize) <- c("group 1", "group 2", "total")

  sampleRatio <- round(x$nratio,2)
  sampleRatio <- t(sampleRatio[c("mean", "min", "q25", "median", "q75", "max")])
  row.names(sampleRatio) <- "group1:group2"

  relFreq <- round(as.matrix(rbind(x$rf1, x$rf2, x$rf)),4)
  relFreq <- relFreq[,c("mean", "min", "median", "max")]
  row.names(relFreq) <- c("group 1", "group 2", "total")

  cat("\nStudies: \n",
      "Total number of studies:       ", x$k, "\n",
      "Number of double-zero studies: ", x$kdz, "\n",
      "Number of single-zero studies: ", x$ksz, ", thereof ", x$k1sz, " zeros in group 1 and ", x$k2sz, " zeros in group 2.\n", sep = "")
  cat("\n")

  cat("\nSample sizes: \n")
  printCoefmat(sampleSize)
  cat("\n")

  cat("\nSample size ratios: \n")
  printCoefmat(sampleRatio)
  cat("\n")

  cat("\nRelative frequencies of the event: \n")
  printCoefmat(relFreq)
  cat("\n")

  cat("\nNumber of studies in which the event is rare (0.01 < rel. freq. < 0.05):", x$krare, "out of", x$k, "studies.")
  cat("\nNumber of studies in which the event is very rare (rel. freq. < 0.01):  ", x$kveryrare, "out of", x$k, "studies.")

}
