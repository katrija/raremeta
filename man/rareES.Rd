% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rareES.R
\name{rareES}
\alias{rareES}
\title{Calculate individual effect sizes for meta-analyses of rare events.}
\usage{
rareES(
  x,
  ai,
  bi,
  ci,
  di,
  n1i,
  n2i,
  data,
  measure,
  cc = "none",
  ccval = 0.5,
  tccval,
  cccval,
  ccsum = 1,
  ccto = "only0",
  drop00 = TRUE,
  method = "FE",
  ...
)
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

\item{measure}{character string specifying the effect size or outcome measure to be used
(either \code{"logOR"} for the log odds ratio, \code{"logRR"} for the log relative risk,
or \code{"RD"} for the risk difference).}

\item{cc}{character string specifying the type of continuity correction to be used
(either \code{"constant"}, \code{"tacc"} or \code{"empirical"}). Default is "none". See 'Details'.}

\item{ccval}{scalar or numerical vector specifying the value of the continuity correction if
\code{cc = "constant"}. Must be a scalar or a vector of length equal to the number of studies.
Default is \code{ccval = 0.5}. If a scalar is specified, the value is added to all studies for
which the number of events is zero in at least one of the groups. This behavior can be changed
using the argument \code{ccto}. \code{ccval} is overwritten by tccval and cccval if both arguments are specified.}

\item{tccval}{scalar or numerical vector specifying the value of the continuity correction
applied to the observations from the treatment group (group 1) if \code{cc = "constant"}. Must be a scalar or a vector
of length equal to the number of studies. If \code{cc = "constant"} and \code{tccval} is not specified, \code{tccval} is
set to the value of \code{ccval} internally.}

\item{cccval}{scalar or numerical vector specifying the value of the continuity correction
applied to the observations from the control group (group 2) if \code{cc = "constant"}. Must be a scalar or a vector
of length equal to the number of studies. If \code{cc = "constant"} and \code{cccval} is not specified, \code{cccval} is
set to the value of \code{ccval} internally.}

\item{ccsum}{numeric value specifying the value of the sum of the continuity correction applied to the
observations from the treatment group and the continuity correction applied to the observations from
the control group. Default is \code{ccsum = 1}. Currently, setting this argument to a different number only has
an effect when \code{cc = "tacc"} or \code{cc = "empirical"}.}

\item{ccto}{character string indicating to which studies the continuity correction should
be applied. Either \code{"only0"}, for which the continuity correction is applied to all studies
for which the number of events is zero in at least one of the groups, \code{"all"}, for which the
continuity correction is applied to all studies, or \code{"if0all"}, for which the continuity
correction is applied to all studies if any of the individual studies has zero events in at
least one of the groups.}

\item{drop00}{logical indicating whether double-zero studies (i.e., studies with no events or only events in both groups)
should be excluded prior to calculating the studies' effect sizes and sampling variances.}

\item{method}{Needs to be specified when \code{cc = "empirical"}. Meta-analytic method to be used when
determining the prior used for the empirical continuity correction. Current default is \code{method = "FE"}.
See Sweeting et al. (2004) for details.}

\item{...}{additional arguments}
}
\value{
an object of class "raremeta". The object is a list containing the following elements:
\itemize{
\item \code{ai}, \code{bi}, \code{ci}, \code{di}: original entries of the 2x2 tables for all studies.
\item \code{measure}: effect size estimand
\item \code{yi}: vector which contains the estimated study effect sizes
\item \code{vi}: vector which contains the estimated sampling variances
\item \code{ai.cc}, \code{bi.cc}, \code{ci.cc}, \code{di.cc}: modified entries of the 2x2 tables for all studies after
application of the specified continuity correction
\item ...
}
}
\description{
Calculate individual effect sizes for meta-analyses of rare events.
}
\section{Details}{
\subsection{Data input}{

The data input can be specified either through the arguments \code{ai},\code{bi},\code{ci},\code{di},\code{n1i}, and \code{n2i} (columns of the data frame \code{data})
or through the argument \code{x}, which takes an object that results from applying the \code{rareDescribe()} function to the data
(i.e., the input for argument \code{x} must be an object of type \code{rareData}).
A \code{rareData} object can be produced from a data frame by applying the \code{rareDescribe()} function to it.
The \code{rareDescribe()} function pre-processes the data frame and stores the information required by the \code{rareES()} function in a list.
See \code{?rareDescribe} for more details.
}

\subsection{Effect size measures}{

The function includes methods for calculating log odds ratios, log risk ratios, and risk differences.
The effect size measure can be specified using the \code{measure} argument. The respective effect size,
along with an estimate of its sampling variance, is then calculated for each study based on the
entries of the study's 2x2 table:\tabular{lrr}{
    \tab event \tab no event \cr
   group1 (treatment) \tab ai \tab bi \cr
   group2 (control) \tab ci \tab di \cr
}

}

\subsection{Handling single-zero and double-zero studies}{

Single-zero studies are studies for which one entry of the 2x2 table is zero. Double-zero studies are
studies for which two entries of the same column of the 2x2 table are zero.
The function includes a variety of arguments to handle single-zero and double-zero studies. Per default,
double-zero studies are currently excluded from the analysis (this behavior might be changed in the future).
The inclusion of double-zero studies can be enforced by setting the argument \code{drop00} to \code{FALSE}.
If the data includes at least one single-zero study, the function throws an error if the user
did not specify whether or not a continuity correction shall be applied.
By setting \code{cc = "none"}, any zero-studies (studies with at least one zero-cell)
which remain in the data are excluded from the analysis. If it is desired that zero-studies are
included, the user needs to specify which continuity correction shall be used, and to which
studies it shall be applied. Per default, the continuity correction is only applied to zero-studies, while
studies for which all cells are larger than zero are left uncorrected. This behavior can be changed
using the argument \code{ccto}. Per default, the constant value 0.5 (\code{cc = "constant"}, \code{ccval = 0.5}) is added to
all cells of the studies specified by \code{ccto}. This continuity correcton was desribed by Gart and Zweifel (1967).
Alternative continuity corrections which were described by Sweeting et al. (2004) can be applied by setting \code{cc} to \code{"tacc"}
for the treatment-arm continuity correction, and to \code{"empirical"} for the empirical continuity correction.
Per default, the sum of the corrections for treatment and control groups is set to \code{1}, but this can be changed by setting the
the argument \code{ccsum} to a different value.
It is possible to set the continuity correction to a user-defined value (or a vector of user-defined values) using the
argument \code{ccval} (if the value). If the user wants to specify different values for the treatment and the control group,
this is possible via the arguments \code{tccval} and \code{cccval}.
}

\subsection{Differences between effect size measures in the application of continuity corrections}{

When either the log odds ratio or the log risk ratio is used as an effect size measure, both the effect sizes and
their sampling variances are calculated based on the continuity-corrected 2x2 table.
When the effect size measure is the risk difference, the continuity-corrected 2x2 table is only used in the
calculation of the sampling variances.
}
}

\examples{

# introducing a dataset
data(dat.nissen2007)
d <- dat.nissen2007


# estimating the log odds ratios with applicatoin of the default continuity correction of 0.5
d.RR <- rareES(ai=miRosiglitazone, ci=miControl, n1i=nRosiglitazone, n2i=nControl, data=d, measure="logRR", cc="constant")
d.RR$yi

# estimating log odds ratios with applicatoin of the empirical continuity correction
d.OR <- rareES(ai=miRosiglitazone, ci=miControl, n1i=nRosiglitazone, n2i=nControl, data=d, measure="logRR", cc="empirical")
d.OR$yi

# same analysis with pre-processed data
x <- rareDescribe(ai=miRosiglitazone, ci=miControl, n1i=nRosiglitazone, n2i=nControl, data=d)
x.OR <- rareES(x, measure="logOR", cc="empirical")
x.OR$yi
}
\references{
Gart, John J, and James R Zweifel. 1967. On the bias of various estimators of the logit and
its variance with application to quantal bioassay. Biometrika, 54, 181–187. doi:10.1093/BIOMET/54.1-2.181

Sweeting, M. J., Sutton, A. J., & Lambert, P. C. (2004). What to add to nothing? Use and avoidance of
continuity corrections in meta-analysis of sparse data. Statistics in Medicine, 23, 1351–1375. doi: 10.1002/sim.1761
}
