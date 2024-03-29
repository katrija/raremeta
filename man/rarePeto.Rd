% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rarePeto.R
\name{rarePeto}
\alias{rarePeto}
\title{Conduct a meta-analysis using Peto's method}
\usage{
rarePeto(x, ai, bi, ci, di, n1i, n2i, data, level = 95, digits = 4)
}
\arguments{
\item{x}{an object of class \code{"rareData"}.}

\item{ai}{data frame column to specify the number of events in group 1 (i.e., the treatment group).}

\item{bi}{data frame column to specify the number of non-events in group 1 (i.e., the treatment group).}

\item{ci}{data frame column to specify the number of events in group 2 (i.e., the control group).}

\item{di}{data frame column to specify number of non-events in group 2 (i.e., the control group).}

\item{n1i}{data frame column to specify the sample sizes in group 1 (i.e., the treatment group).}

\item{n2i}{data frame column to specify the sample sizes in group 2 (i.e., the control group).}

\item{data}{data frame.}

\item{level}{numeric nbetween 0 and 100 specifying the confidence interval
level (the default is 95)}

\item{digits}{integer specifying the number of decimal places to which the printed results
should be rounded (if unspecified, the default is 4).}
}
\value{
an object of class "raremeta".
The object is a list containing the following elements:
\itemize{
\item \code{beta}, \code{b}: estimated effect size.
\item \code{se}: standard error of the  estimator.
\item \code{zval}: test statistics of the coefficients.
\item \code{pval}: p-values corresponding to the test statistics.
\item \code{ci.lb}: lower bound of the confidence intervals for the coefficients.
\item \code{ci.ub}: upper bound of the confidence intervals for the coefficients.
\item \code{k}: number of studies included in the analysis.
\item \code{kdz},\code{ksz}: number of double-zero and single-zero studies.
\item \code{k1sz}, \code{k2sz}: number of single-zero studies where the zero is in group 1 or group 2.
\item \code{ai}, \code{bi}, \code{ci}, \code{di}: original entries of the 2x2 tables for all studies.
\item \code{ni}, \code{n1i}, \code{n2i}: original total and group sample sizes.
\item ...
}
}
\description{
Function to conduct a meta-analysis of the log odds ratio using the Peto method.
}
\section{Details}{
\subsection{Data input}{

The data input can be specified either through the arguments \code{ai},\code{bi},\code{ci},\code{di},\code{n1i}, and \code{n2i} (columns of the data frame \code{data})
or through the argument \code{x}, which takes an object that results from applying the \code{rareDescribe()} function to the data
(i.e., the input for argument \code{x} must be an object of type \code{rareData}).
A \code{rareData} object can be produced from a data frame by applying the \code{rareDescribe()} function to it.
The \code{rareDescribe()} function pre-processes the data frame and stores the information required by the \code{rarePeto()} function in a list.
See \code{?rareDescribe} for more details.
}

\subsection{Effect size measure}{

The Peto method is a fixed-effects method for meta-analysis of the log odds ratio. For further details,
see Yusuf et al. (1985).
}
}

\examples{
# load a data set
data(dat.nissen2007)
d <- dat.nissen2007

# meta analysis using Peto's method
rarePeto(ai=miRosiglitazone, ci=miControl, n1i=nRosiglitazone, n2i=nControl, data=d)

# same analysis with pre-processed data
x <- rareDescribe(ai=miRosiglitazone, ci=miControl, n1i=nRosiglitazone, n2i=nControl, data=d)
rarePeto(x)

}
\references{
Yusuf, S., Peto, R., Lewis, J., Collins, R., & Sleight, P. (1985).
Beta blockade during and after myocardial infarction: an overview of
the randomized trials. Progress in cardiovascular diseases, 27(5), 335-371.
}
