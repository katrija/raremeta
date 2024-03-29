% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rareDescribe.R
\name{rareDescribe}
\alias{rareDescribe}
\title{Descriptives for a meta-analysis of a rare event}
\usage{
rareDescribe(ai, bi, ci, di, n1i, n2i, data)
}
\arguments{
\item{ai}{data frame column to specify the number of events in group 1 (i.e., the treatment group).}

\item{bi}{data frame column to specify the number of non-events in group 1 (i.e., the treatment group).}

\item{ci}{data frame column to specify the number of events in group 2 (i.e., the control group).}

\item{di}{data frame column to specify number of non-events in group 2 (i.e., the control group).}

\item{n1i}{data frame column to specify the sample sizes in group 1 (i.e., the treatment group).}

\item{n2i}{data frame column to specify the sample sizes in group 2 (i.e., the control group).}

\item{data}{data frame.}
}
\value{
An object of the class \code{"rareData"}, which contains descriptives of the meta-analytic data.
}
\description{
Function to compute descriptives based on the 2x2 tables of the individual studies
for meta-analytic rare event data.
}
\examples{
# introducing the data
data(dat.nissen2007)
d <- dat.nissen2007

# pre-processing the data for further usage
x <- rareDescribe(ai=miRosiglitazone, ci=miControl, n1i=nRosiglitazone, n2i=nControl, data=d)
x


}
