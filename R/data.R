#' Meta-analytic data from Yuan et al. (2018)
#'
#' A data frame with 8 rows and 5 columns.
#' * study: study identifier
#' * ai: number of events (dropout) in the treatment group
#' * ci: number of events (dropout) in the control group
#' * n1i: sample size of the treatment group
#' * n2i: sample size of the control group
#'
#' @docType data
#'
#' @usage dat.yuan2018
#'
#' @examples
#' # load data frame:
#'
#' data <- dat.yuan2018
#'
#' # calculate descriptives using rareDescribe():
#'
#' descr <- rareDescribe(ai = ai,
#'                       n1i = n1i,
#'                       ci = ci,
#'                       n2i = n2i,
#'                       data = data)
#'
#' print(descr)
#'
#' # conduct random-effects meta-analysis
#' # using the inverse variance model:
#'
#' ma <- rareIV(descr,
#'              measure = "logOR",
#'              method = "REML",
#'              cc = "constant",
#'              ccto = "only0",
#'              ccval = 0.5,
#'              drop00 = FALSE)
#'
#' summary(ma)
#'
#' @references Yuan, S., Zhou, X., Zhang, Y., Zhang, H., Pu, J., Yang, L., Liu, L., Jiang, X., & Xie, P. (2018). Comparative efficacy
#' and acceptability of bibliotherapy for depression and anxiety disorders in children and adolescents: A meta-analysis
#' of randomized clinical trials. Neuropsychiatric Disease and Treatment, 14, 353–365. \href{https://doi.org/10.2147/NDT.S152747}{NDT}.
"dat.yuan2018"


#' Meta-analytic data from Nissen and Wolski (2007)
#'
#' A data frame with 42 rows and 8 columns.
#' * study: study identifier
#' * nRosiglitazone: number of patients in the Rosiglitazone group
#' * miRosiglitazone: number of myocardial infarctions in the Rosiglitazone group
#' * cvRosiglitazone: number of deaths from cardiovascular causes in the Rosiglitazone group
#' * nControl: number of patients in the control group
#' * miControl: number of myocardial infarctions in the control group
#' * cvControl: number of deaths from cardiovascular causes in the control group
#' * methodControl: type of control group
#'
#' @docType data
#'
#' @usage dat.nissen2007
#'
#' @examples
#' # load data frame:
#'
#' data <- dat.nissen2007
#'
#' # calculate descriptives using rareDescribe():
#'
#' descr <- rareDescribe(ai = miRosiglitazone,
#'                       n1i = nRosiglitazone,
#'                       ci = miControl,
#'                       n2i = nControl,
#'                       data = data)
#'
#' print(descr)
#'
#' # conduct random-effects meta-analysis
#' # using the inverse variance model:
#'
#' ma <- rareIV(descr,
#'              measure = "logOR",
#'              method = "REML",
#'              cc = "constant",
#'              ccto = "only0",
#'              ccval = 0.5,
#'              drop00 = FALSE)
#'
#' summary(ma)
#'
#' @references Nissen and Wolski (2007) (\href{https://www.nejm.org/doi/full/10.1056/NEJMoa072761}{NEJM})
"dat.nissen2007"


#' Meta-analytic data from Nissen and Wolski (2010)
#'
#' A data frame with 56 rows and 10 columns.
#' * study: study identifier
#' * durationWeeks: study duration in weeks
#' * nRosiglitazone: number of patients in the Rosiglitazone group
#' * miRosiglitazone: number of myocardial infarctions in the Rosiglitazone group
#' * cvRosiglitazone: number of deaths from cardiovascular causes in the Rosiglitazone group
#' * nControl: number of patients in the control group
#' * miControl: number of myocardial infarctions in the control group
#' * cvControl: number of deaths from cardiovascular causes in the control group
#' * methodControl: type of control group
#' * durationRosiglitazone: person times of the Rosiglitazone group (durationWeeks &ast; nRosiglitazone)
#' * durationcontrol: person times of the control group (durationWeeks &ast; nControl)
#'
#' @docType data
#'
#' @usage dat.nissen2010
#'
#' @examples
#' # load data frame:
#'
#' data <- dat.nissen2010
#'
#' # calculate descriptives using rareDescribe():
#'
#' descr <- rareDescribe(ai = miRosiglitazone,
#'                       n1i = nRosiglitazone,
#'                       ci = miControl,
#'                       n2i = nControl,
#'                       data = data)
#'
#' print(descr)
#'
#' # conduct random-effects meta-analysis
#' # using the inverse variance model:
#'
#' ma <- rareIV(descr,
#'              measure = "logRR",
#'              method = "REML",
#'              cc = "constant",
#'              ccto = "only0",
#'              ccval = 0.5,
#'              drop00 = FALSE)
#'
#' summary(ma)
#'
#' @format A data frame with 56 rows and 10 columns.
#'
#' @references Nissen and Wolski (2010) (\href{https://jamanetwork.com/journals/jamainternalmedicine/article-abstract/225844}{JAMA})
"dat.nissen2010"
