#' Meta-analytic data from Hoppen et al. (2023): EMDR vs. passive control conditions
#'
#' A data frame with 7 rows and 5 columns.
#' * study: study identifier
#' * ai: number of events (dropout) in the treatment group (EMDR)
#' * ci: number of events (dropout) in the control group (passive control conditions)
#' * n1i: sample size of the treatment group
#' * n2i: sample size of the control group
#'
#' @docType data
#'
#' @usage dat.hoppen2023
#'
#' @examples
#' # load data frame:
#'
#' data <- dat.hoppen2023
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
#' @references Hoppen, T. H., Jehn, M., Holling, H., Mutz, J., Kip, A., & Morina, N. (2023).
#' The efficacy and acceptability of psychological interventions for adult PTSD:
#' A network and pairwise meta-analysis of randomized controlled trials.
#' Journal of Consulting and Clinical Psychology.
#' \href{https://doi.org/10.1037/ccp0000809}{J CONSULT CLIN PSYCH}.
"dat.hoppen2023"

#' Meta-analytic data from Yang et al. (2019): Psychological interventions for social anxiety disorder vs. control
#'
#' A data frame with 16 rows and 5 columns.
#' * study: study identifier
#' * ai: number of events (dropout) in the treatment group (psychological intervention group)
#' * ci: number of events (dropout) in the control group (control group)
#' * n1i: sample size of the treatment group
#' * n2i: sample size of the control group
#'
#' @docType data
#'
#' @usage dat.yang2019
#'
#' @examples
#' # load data frame:
#'
#' data <- dat.yang2019
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
#' @references Yang, L., Zhou, X., Pu, J., Liu, L., Cuijpers, P., Zhang, Y., ... & Xie, P. (2019).
#' Efficacy and acceptability of psychological interventions for social anxiety disorder in
#' children and adolescents: a meta-analysis of randomized controlled trials.
#' European Child & Adolescent Psychiatry, 28, 79-89.
#' \href{https://doi.org/10.1007/s00787-018-1189-x}{EUR CHILD ADOLES PSY}.
"dat.yang2019"

#' Meta-analytic data from Okumura et al. (2014): Group CBT for depression vs. non-active control
#'
#' A data frame with 27 rows and 5 columns.
#' * study: study identifier
#' * ai: number of events (dropout) in the treatment group (group CBT group)
#' * ci: number of events (dropout) in the control group (non-active control group)
#' * n1i: sample size of the treatment group
#' * n2i: sample size of the control group
#'
#' @docType data
#'
#' @usage dat.okumura2014
#'
#' @examples
#' # load data frame:
#'
#' data <- dat.okumura2014
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
#' @references Okumura, Y., & Ichikura, K. (2014). Efficacy and acceptability of
#' group cognitive behavioral therapy for depression: a systematic review and
#' meta-analysis. Journal of Affective Disorders, 164, 155-164.
#' \href{https://doi.org/10.1016/j.jad.2014.04.023}{JAD}.
"dat.okumura2014"

#' Meta-analytic data from Fodor et al. (2018): VR therapy vs. passive control groups
#'
#' A data frame with 26 rows and 5 columns.
#' * study: study identifier
#' * ai: number of events (dropout) in the treatment group (VR therapy group)
#' * ci: number of events (dropout) in the control group (passive control group)
#' * n1i: sample size of the treatment group
#' * n2i: sample size of the control group
#'
#' @docType data
#'
#' @usage dat.fodor2018.pass
#'
#' @examples
#' # load data frame:
#'
#' data <- dat.fodor2018.pass
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
#' @references Fodor, L. A., Coteț, C. D., Cuijpers, P., Szamoskozi, Ș., David, D., & Cristea, I. A. (2018).
#' The effectiveness of virtual reality based interventions for symptoms of anxiety and depression: A meta-analysis.
#' Scientific reports, 8(1), 10323. \href{https://doi.org/10.1038/s41598-018-28113-6}{SciRep}.
"dat.fodor2018.pass"

#' Meta-analytic data from Fodor et al. (2018): VR therapy vs. active control groups
#'
#' A data frame with 27 rows and 5 columns.
#' * study: study identifier
#' * ai: number of events (dropout) in the treatment group (VR therapy group)
#' * ci: number of events (dropout) in the control group (active control group)
#' * n1i: sample size of the treatment group
#' * n2i: sample size of the control group
#'
#' @docType data
#'
#' @usage dat.fodor2018.act
#'
#' @examples
#' # load data frame:
#'
#' data <- dat.fodor2018.act
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
#' @references Fodor, L. A., Coteț, C. D., Cuijpers, P., Szamoskozi, Ș., David, D., & Cristea, I. A. (2018).
#' The effectiveness of virtual reality based interventions for symptoms of anxiety and depression: A meta-analysis.
#' Scientific reports, 8(1), 10323. \href{https://doi.org/10.1038/s41598-018-28113-6}{SciRep}.
"dat.fodor2018.act"


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
