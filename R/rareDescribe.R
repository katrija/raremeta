rareDescribe <- function(ai, bi, ci, di, n1i, n2i,
                         data){

  mf <- match.call()

  mf.ai <- mf[[match("ai", names(mf))]]
  mf.bi <- mf[[match("bi", names(mf))]]
  mf.ci <- mf[[match("ci", names(mf))]]
  mf.di <- mf[[match("di", names(mf))]]
  mf.n1i <- mf[[match("n1i", names(mf))]]
  mf.n2i <- mf[[match("n2i", names(mf))]]

  ai <- eval(mf.ai, data)
  bi <- eval(mf.bi, data)
  ci <- eval(mf.ci, data)
  di <- eval(mf.di, data)
  n1i <- eval(mf.n1i, data)
  n2i <- eval(mf.n2i, data)

  if(!is.data.frame(data)){
    stop("First argument must be a data frame.")
  }

  nEventsT <- sum(ai, na.rm = TRUE)
  nEventsC <- sum(ci, na.rm = TRUE)

  descEvents <- paste("There are", nEventsT, "events in the treatment group and", nEventsC, "in the control group.")

  return(descEvents)
}
