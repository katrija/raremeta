# usual data set:
data <- data.frame(
  ai = c(0, 3, 2, 0),
  bi = c(20, 18, 15, 19),
  ci = c(1, 4, 0, 0),
  di = c(19, 17, 16, 20),
  n1i = c(20, 21, 17, 19),
  n2i = c(20, 21, 16, 20)
)

x <- rareDescribe(
  ai = ai,
  bi = bi,
  ci = ci,
  di = di,
  n1i = n1i,
  n2i = n2i,
  data = data
)

# same data set but double zero studies removed
data_rareIV00 <- data.frame(
  ai = c(0, 3, 2),
  bi = c(20, 18, 15),
  ci = c(1, 4, 0),
  di = c(19, 17, 16),
  n1i = c(20, 21, 17),
  n2i = c(20, 21, 16)
)

x_rareIV00 <- rareDescribe(
  ai = ai,
  bi = bi,
  ci = ci,
  di = di,
  n1i = n1i,
  n2i = n2i,
  data = data_rareIV00
)

# data set with only zero-studies
data_only0 <- data.frame(
  ai = c(1, 0, 0, 1),
  bi = c(0, 2, 0, 0),
  ci = c(0, 4, 1, 0),
  di = c(19, 0, 0, 20),
  n1i = c(1, 2, 0, 1),
  n2i = c(19, 4, 1, 20)
)

x_only0 <- rareDescribe(
  ai = ai,
  bi = bi,
  ci = ci,
  di = di,
  n1i = n1i,
  n2i = n2i,
  data = data_only0
)

test_that("rareIV returns errors and warning messages", {

  # error if x is not object of class rareData
  expect_error(
    rareIV(x = data, measure = "logOR", method = "FE", cc = "constant"),
    "x must be an object of class 'rareData'."
  )


  # error if measure argument is not defined
  expect_error(
    rareIV(x, method = "FE", cc = "constant"),
    "'measure' argument must be specified."
  )

  # error if measure argument is not valid
  expect_error(
    rareIV(x, measure = "logoR", method = "FE", cc = "constant"),
    "'measure' must be either 'logOR', 'logRR', or 'RD'."
  )

  # error if method argument is not defined
  expect_error(
    rareIV(x, measure = "logOR", cc = "constant"),
    "'method' argument must be specified."
  )

  # error if method argument is not valid
  expect_error(
    rareIV(x, measure = "logOR", method = "fe", cc = "constant"),
    "unknown 'method' specified."
  )

  # error if cc is not specified
  expect_error(
    rareIV(x, measure = "logOR", method = "FE"),
    "Some studies have zero events."
  )

  # error if drop00 is not logical
  expect_error(
    rareIV(x, measure = "logOR", method = "FE", cc = "constant", drop00 = "true"),
    "'drop00' must be a logical"
  )

  # error if cc argument is not valid
  expect_error(
    rareIV(x, measure = "logOR", method = "FE", cc = "None"),
    "'cc' must be either 'none', 'constant', 'tacc', or 'empirical'"
  )

  # error if ccto argument is not valid
  expect_error(
    rareIV(x, measure = "logOR", method = "FE", cc = "constant", ccto = "Only0"),
    "'ccto' must be either 'only0', 'all', or 'if0all'"
  )

  # error if ccval argument is not valid
  ccval <- c(1, 1)

  # for drop00 = FALSE
  expect_error(
    rareIV(x, measure = "logOR", method = "FE", cc = "constant", drop00 = FALSE, ccval = ccval),
    "'ccval' must have length 1 or length equal to the number of studies."
  )

  # for drop00 = TRUE
  expect_error(
    rareIV(x, measure = "logOR", method = "FE", cc = "constant", drop00 = TRUE, ccval = ccval),
    "'ccval' must have length 1 or length equal to the number of studies"
  )

  # error if ccval argument is not non-negative
  expect_error(
    rareIV(x, measure = "logOR", method = "FE", cc = "constant", ccval = -1),
    "'ccval' must be non-negative."
  )

  # error if only tccval is defined
  expect_error(
    rareIV(x, measure = "logOR", method = "FE", cc = "constant", tccval = 0.5),
    "Please specify both 'tccval' and 'cccval'."
  )

  # error if only cccval is defined
  expect_error(
    rareIV(x, measure = "logOR", method = "FE", cc = "constant", cccval = 0.5),
    "Please specify both 'tccval' and 'cccval'."
  )

  # error if tccval and cccval have different lengths
  expect_error(
    rareIV(x, measure = "logOR", method = "FE", cc = "constant", cccval = c(0.5, 0.5, 0.5, 0.1), tccval = c(0.5, 0.5, 0.1)),
    "'tccval' and 'cccval' must have equal length."
  )

  # error if tccval is not length 1 or length equal to number of studies for drop00 = FALSE
  expect_error(
    rareIV(x, measure = "logOR", method = "FE", cc = "constant", drop00 = FALSE, cccval = c(0.5, 0.5, 0.5), tccval = c(0.5, 0.5, 0.1)),
    "'tccval' must have length 1 or length equal to the number of studies."
  )

  # # error if tccval is not length 1 or length equal to number of studies for drop00 = TRUE
  #
  expect_error(
    rareIV(x,
      measure = "logOR", method = "FE", cc = "constant", drop00 = TRUE, cccval = c(0.5, 0.5),
      tccval = c(0.5, 0.5)
    ),
    "'tccval' must have length 1 or length equal to the number of studies"
  )

  # error if tccval is not non-negative
  expect_error(
    rareIV(x, measure = "logOR", method = "FE", cc = "constant", cccval = 0.1, tccval = -0.5),
    "All values in 'tccval' and 'cccval' must be non-negative."
  )

  # error if cccval is not non-negative
  expect_error(
    rareIV(x, measure = "logOR", method = "FE", cc = "constant", tccval = 0.1, cccval = -0.5),
    "All values in 'tccval' and 'cccval' must be non-negative."
  )

  # error if test argument is not valid
  expect_error(
    rareIV(x, measure = "logOR", method = "FE", cc = "constant", test = "Z"),
    "'test' must be either 'z', 'knha', or 'hksj'"
  )

  # error if digits argument is not valid
  expect_error(
    rareIV(x, measure = "logOR", method = "FE", cc = "constant", digits = 0.5),
    "'digits' must be an integer of length 1."
  )

  # error if cc = tacc & measure = RD
  expect_error(
    rareIV(x, measure = "RD", method = "FE", cc = "tacc"),
    "continuity correction of type 'tacc' is currently not supported for measure 'RD'."
  )

  # error if cc = empirical & measure = RD
  expect_error(
    rareIV(x, measure = "RD", method = "FE", cc = "empirical"),
    "continuity correction of type 'empirical' is currently not supported for measure 'RD'."
  )

  # error if there are no non-zero studies and cc = empirical
  expect_error(
    rareIV(x_only0, measure = "logOR", method = "FE", cc = "empirical"),
    "continuity correction of type 'empirical' can only be applied if there is at least one non-zero study."
  )

  # error if cc = empirical & method = IPM
  expect_error(
    rareIV(x, measure = "logOR", method = "IPM", cc = "empirical"),
    "method = 'IPM' is currently not defined for cc = 'empirical'."
  )

  # error if ccsum is a scalar less than 0
  expect_error(
    rareIV(x, measure = "logOR", method = "FE", cc = "tacc", ccsum = -1),
    "ccsum must be a scalar larger than 0."
  )

  # error if level is a scalar larger than 100
  expect_error(
    rareIV(x, measure = "logOR", method = "FE", cc = "constant", level = 101),
    "level must be a scalar between 0 and 100."
  )

  # error if method = IPM & measure != logOR
  expect_error(
    rareIV(x, measure = "logRR", method = "IPM", cc = "constant"),
    "method = 'IPM' is only defined for measure = 'logOR'."
  )

  # # warning if method = IPM & cc = constant & ccto = only0
  expect_warning(
    rareIV(x, measure = "logOR", method = "IPM", cc = "constant", ccto = "only0"),
    "method = 'IPM' was developed for cc = 'constant' and ccto = 'all'."
  )
})

test_that("results from rareIV equal results from metafor package", {

  #############################################################################
  ##################### method = FE ###########################################
  ##################### measure = logOR #######################################
  #############################################################################
  # cc = constant, add 0.5
  # for method = FE
  # for measure = logOR
  # add to "only0"
  # drop00 = FALSE
  rare <- rareIV(x, measure = "logOR", method = "FE", cc = "constant", ccval = 0.5, ccto = "only0", drop00 = FALSE, level = 95)
  fit <- metafor::rma(ai = ai, bi = bi, ci = ci, di = di, n1i = n1i, n2i = n2i, data = data, measure = "OR", method = "FE", add = 0.5, to = "only0", drop00 = FALSE, level = 95)
  expect_equal(rare$beta, fit$beta)
  expect_equal(rare$tau2, fit$tau2)
  expect_equal(rare$se.tau2, fit$se.tau2)
  expect_equal(rare$se, fit$se)
  expect_equal(rare$pval, fit$pval)
  expect_equal(rare$ci.lb, fit$ci.lb)
  expect_equal(rare$ci.ub, fit$ci.ub)

  # cc = constant, add 0.5
  # for method = FE
  # for measure = logOR
  # add to "only0"
  # drop00 = TRUE
  rare <- rareIV(x, measure = "logOR", method = "FE", cc = "constant", ccval = 0.5, ccto = "only0", drop00 = TRUE, level = 95)
  fit <- suppressWarnings(metafor::rma(ai = ai, bi = bi, ci = ci, di = di, n1i = n1i, n2i = n2i, data = data, measure = "OR", method = "FE", add = 0.5, to = "only0", drop00 = TRUE, level = 95))
  expect_equal(rare$beta, fit$beta)
  expect_equal(rare$tau2, fit$tau2)
  expect_equal(rare$se.tau2, fit$se.tau2)
  expect_equal(rare$se, fit$se)
  expect_equal(rare$pval, fit$pval)
  expect_equal(rare$ci.lb, fit$ci.lb)
  expect_equal(rare$ci.ub, fit$ci.ub)

  # cc = constant, add 0.5
  # for method = FE
  # for measure = logOR
  # add to "all"
  # drop00 = FALSE
  rare <- rareIV(x, measure = "logOR", method = "FE", cc = "constant", ccval = 0.5, ccto = "all", drop00 = FALSE, level = 95)
  fit <- metafor::rma(ai = ai, bi = bi, ci = ci, di = di, n1i = n1i, n2i = n2i, data = data, measure = "OR", method = "FE", add = 0.5, to = "all", drop00 = FALSE, level = 95)
  expect_equal(rare$beta, fit$beta)
  expect_equal(rare$tau2, fit$tau2)
  expect_equal(rare$se.tau2, fit$se.tau2)
  expect_equal(rare$se, fit$se)
  expect_equal(rare$pval, fit$pval)
  expect_equal(rare$ci.lb, fit$ci.lb)
  expect_equal(rare$ci.ub, fit$ci.ub)
  # warning message: NAs omitted from model fitting for metafor::rma

  # cc = constant, add 0.5
  # for method = FE
  # for measure = logOR
  # add to "all"
  # drop00 = TRUE
  rare <- rareIV(x, measure = "logOR", method = "FE", cc = "constant", ccval = 0.5, ccto = "all", drop00 = TRUE, level = 95)
  fit <- suppressWarnings(metafor::rma(ai = ai, bi = bi, ci = ci, di = di, n1i = n1i, n2i = n2i, data = data, measure = "OR", method = "FE", add = 0.5, to = "all", drop00 = TRUE, level = 95))
  expect_equal(rare$beta, fit$beta)
  expect_equal(rare$tau2, fit$tau2)
  expect_equal(rare$se.tau2, fit$se.tau2)
  expect_equal(rare$se, fit$se)
  expect_equal(rare$pval, fit$pval)
  expect_equal(rare$ci.lb, fit$ci.lb)
  expect_equal(rare$ci.ub, fit$ci.ub)
  # warning message: NAs omitted from model fitting for metafor::rma

  # cc = constant, add 0.5
  # for method = FE
  # for measure = logOR
  # add to "if0all"
  # drop00 = FALSE
  rare <- rareIV(x, measure = "logOR", method = "FE", cc = "constant", ccval = 0.5, ccto = "if0all", drop00 = FALSE, level = 95)
  fit <- metafor::rma(ai = ai, bi = bi, ci = ci, di = di, n1i = n1i, n2i = n2i, data = data, measure = "OR", method = "FE", add = 0.5, to = "if0all", drop00 = FALSE, level = 95)
  expect_equal(rare$beta, fit$beta)
  expect_equal(rare$tau2, fit$tau2)
  expect_equal(rare$se.tau2, fit$se.tau2)
  expect_equal(rare$se, fit$se)
  expect_equal(rare$pval, fit$pval)
  expect_equal(rare$ci.lb, fit$ci.lb)
  expect_equal(rare$ci.ub, fit$ci.ub)

  # cc = constant, add 0.5
  # for method = FE
  # for measure = logOR
  # add to "if0all"
  # drop00 = TRUE
  rare <- rareIV(x, measure = "logOR", method = "FE", cc = "constant", ccval = 0.5, ccto = "if0all", drop00 = TRUE, level = 95)
  fit <- suppressWarnings(metafor::rma(ai = ai, bi = bi, ci = ci, di = di, n1i = n1i, n2i = n2i, data = data, measure = "OR", method = "FE", add = 0.5, to = "if0all", drop00 = TRUE, level = 95))
  expect_equal(rare$beta, fit$beta)
  expect_equal(rare$tau2, fit$tau2)
  expect_equal(rare$se.tau2, fit$se.tau2)
  expect_equal(rare$se, fit$se)
  expect_equal(rare$pval, fit$pval)
  expect_equal(rare$ci.lb, fit$ci.lb)
  expect_equal(rare$ci.ub, fit$ci.ub)
  # warning message: NAs omitted from model fitting for metafor::rma


  ##############################################################################
  ##################### measure = logRR #######################################

  # cc = constant, add 0.5
  # for method = FE
  # for measure = logRR
  # add to "only0"
  # drop00 = FALSE
  rare <- rareIV(x, measure = "logRR", method = "FE", cc = "constant", ccval = 0.5, ccto = "only0", drop00 = FALSE, level = 95)
  fit <- metafor::rma(ai = ai, bi = bi, ci = ci, di = di, n1i = n1i, n2i = n2i, data = data, measure = "RR", method = "FE", add = 0.5, to = "only0", drop00 = FALSE, level = 95)
  expect_equal(rare$beta, fit$beta)
  expect_equal(rare$tau2, fit$tau2)
  expect_equal(rare$se.tau2, fit$se.tau2)
  expect_equal(rare$se, fit$se)
  expect_equal(rare$pval, fit$pval)
  expect_equal(rare$ci.lb, fit$ci.lb)
  expect_equal(rare$ci.ub, fit$ci.ub)

  # cc = constant, add 0.5
  # for method = FE
  # for measure = logRR
  # add to "only0"
  # drop00 = TRUE
  rare <- rareIV(x, measure = "logRR", method = "FE", cc = "constant", ccval = 0.5, ccto = "only0", drop00 = TRUE, level = 95)
  fit <- suppressWarnings(metafor::rma(ai = ai, bi = bi, ci = ci, di = di, n1i = n1i, n2i = n2i, data = data, measure = "RR", method = "FE", add = 0.5, to = "only0", drop00 = TRUE, level = 95))
  expect_equal(rare$beta, fit$beta)
  expect_equal(rare$tau2, fit$tau2)
  expect_equal(rare$se.tau2, fit$se.tau2)
  expect_equal(rare$se, fit$se)
  expect_equal(rare$pval, fit$pval)
  expect_equal(rare$ci.lb, fit$ci.lb)
  expect_equal(rare$ci.ub, fit$ci.ub)
  # warning message: NAs omitted from model fitting for metafor::rma


  # cc = constant, add 0.5
  # for method = FE
  # for measure = logRR
  # add to "all"
  # drop00 = FALSE
  rare <- rareIV(x, measure = "logRR", method = "FE", cc = "constant", ccval = 0.5, ccto = "all", drop00 = FALSE, level = 95)
  fit <- metafor::rma(ai = ai, bi = bi, ci = ci, di = di, n1i = n1i, n2i = n2i, data = data, measure = "RR", method = "FE", add = 0.5, to = "all", drop00 = FALSE, level = 95)
  expect_equal(rare$beta, fit$beta)
  expect_equal(rare$tau2, fit$tau2)
  expect_equal(rare$se.tau2, fit$se.tau2)
  expect_equal(rare$se, fit$se)
  expect_equal(rare$pval, fit$pval)
  expect_equal(rare$ci.lb, fit$ci.lb)
  expect_equal(rare$ci.ub, fit$ci.ub)

  # cc = constant, add 0.5
  # for method = FE
  # for measure = logRR
  # add to "all"
  # drop00 = TRUE
  rare <- rareIV(x, measure = "logRR", method = "FE", cc = "constant", ccval = 0.5, ccto = "all", drop00 = TRUE, level = 95)
  fit <- suppressWarnings(metafor::rma(ai = ai, bi = bi, ci = ci, di = di, n1i = n1i, n2i = n2i, data = data, measure = "RR", method = "FE", add = 0.5, to = "all", drop00 = TRUE, level = 95))
  expect_equal(rare$beta, fit$beta)
  expect_equal(rare$tau2, fit$tau2)
  expect_equal(rare$se.tau2, fit$se.tau2)
  expect_equal(rare$se, fit$se)
  expect_equal(rare$pval, fit$pval)
  expect_equal(rare$ci.lb, fit$ci.lb)
  expect_equal(rare$ci.ub, fit$ci.ub)
  # warning message: NAs omitted from model fitting for metafor::rma


  # cc = constant, add 0.5
  # for method = FE
  # for measure = logRR
  # add to "if0all"
  # drop00 = FALSE
  rare <- rareIV(x, measure = "logRR", method = "FE", cc = "constant", ccval = 0.5, ccto = "if0all", drop00 = FALSE, level = 95)
  fit <- metafor::rma(ai = ai, bi = bi, ci = ci, di = di, n1i = n1i, n2i = n2i, data = data, measure = "RR", method = "FE", add = 0.5, to = "if0all", drop00 = FALSE, level = 95)
  expect_equal(rare$beta, fit$beta)
  expect_equal(rare$tau2, fit$tau2)
  expect_equal(rare$se.tau2, fit$se.tau2)
  expect_equal(rare$se, fit$se)
  expect_equal(rare$pval, fit$pval)
  expect_equal(rare$ci.lb, fit$ci.lb)
  expect_equal(rare$ci.ub, fit$ci.ub)

  # cc = constant, add 0.5
  # for method = FE
  # for measure = logRR
  # add to "if0all"
  # drop00 = TRUE
  rare <- rareIV(x, measure = "logRR", method = "FE", cc = "constant", ccval = 0.5, ccto = "if0all", drop00 = TRUE, level = 95)
  fit <- suppressWarnings(metafor::rma(ai = ai, bi = bi, ci = ci, di = di, n1i = n1i, n2i = n2i, data = data, measure = "RR", method = "FE", add = 0.5, to = "if0all", drop00 = TRUE, level = 95))
  expect_equal(rare$beta, fit$beta)
  expect_equal(rare$tau2, fit$tau2)
  expect_equal(rare$se.tau2, fit$se.tau2)
  expect_equal(rare$se, fit$se)
  expect_equal(rare$pval, fit$pval)
  expect_equal(rare$ci.lb, fit$ci.lb)
  expect_equal(rare$ci.ub, fit$ci.ub)
  # warning message: NAs omitted from model fitting for metafor::rma

  #############################################################################
  ##################### method = DL ###########################################
  ##################### measure = logOR #######################################
  #############################################################################
  # cc = constant, add 0.5
  # for method = DL
  # for measure = logOR
  # add to "only0"
  # drop00 = FALSE
  rare <- rareIV(x, measure = "logOR", method = "DL", cc = "constant", ccval = 0.5, ccto = "only0", drop00 = FALSE, level = 95)
  fit <- metafor::rma(ai = ai, bi = bi, ci = ci, di = di, n1i = n1i, n2i = n2i, data = data, measure = "OR", method = "DL", add = 0.5, to = "only0", drop00 = FALSE, level = 95)
  expect_equal(rare$beta, fit$beta)
  expect_equal(rare$tau2, fit$tau2)
  expect_equal(rare$se.tau2, fit$se.tau2)
  expect_equal(rare$se, fit$se)
  expect_equal(rare$pval, fit$pval)
  expect_equal(rare$ci.lb, fit$ci.lb)
  expect_equal(rare$ci.ub, fit$ci.ub)

  # cc = constant, add 0.5
  # for method = DL
  # for measure = logOR
  # add to "only0"
  # drop00 = TRUE
  rare <- rareIV(x, measure = "logOR", method = "DL", cc = "constant", ccval = 0.5, ccto = "only0", drop00 = TRUE, level = 95)
  fit <- suppressWarnings(metafor::rma(ai = ai, bi = bi, ci = ci, di = di, n1i = n1i, n2i = n2i, data = data, measure = "OR", method = "DL", add = 0.5, to = "only0", drop00 = TRUE, level = 95))
  expect_equal(rare$beta, fit$beta)
  expect_equal(rare$tau2, fit$tau2)
  expect_equal(rare$se.tau2, fit$se.tau2)
  expect_equal(rare$se, fit$se)
  expect_equal(rare$pval, fit$pval)
  expect_equal(rare$ci.lb, fit$ci.lb)
  expect_equal(rare$ci.ub, fit$ci.ub)
  # warning message: NAs omitted from model fitting for metafor::rma


  # cc = constant, add 0.5
  # for method = DL
  # for measure = logOR
  # add to "all"
  # drop00 = FALSE
  rare <- rareIV(x, measure = "logOR", method = "DL", cc = "constant", ccval = 0.5, ccto = "all", drop00 = FALSE, level = 95)
  fit <- metafor::rma(ai = ai, bi = bi, ci = ci, di = di, n1i = n1i, n2i = n2i, data = data, measure = "OR", method = "DL", add = 0.5, to = "all", drop00 = FALSE, level = 95)
  expect_equal(rare$beta, fit$beta)
  expect_equal(rare$tau2, fit$tau2)
  expect_equal(rare$se.tau2, fit$se.tau2)
  expect_equal(rare$se, fit$se)
  expect_equal(rare$pval, fit$pval)
  expect_equal(rare$ci.lb, fit$ci.lb)
  expect_equal(rare$ci.ub, fit$ci.ub)

  # cc = constant, add 0.5
  # for method = DL
  # for measure = logOR
  # add to "all"
  # drop00 = TRUE
  rare <- rareIV(x, measure = "logOR", method = "DL", cc = "constant", ccval = 0.5, ccto = "all", drop00 = TRUE, level = 95)
  fit <- suppressWarnings(metafor::rma(ai = ai, bi = bi, ci = ci, di = di, n1i = n1i, n2i = n2i, data = data, measure = "OR", method = "DL", add = 0.5, to = "all", drop00 = TRUE, level = 95))
  expect_equal(rare$beta, fit$beta)
  expect_equal(rare$tau2, fit$tau2)
  expect_equal(rare$se.tau2, fit$se.tau2)
  expect_equal(rare$se, fit$se)
  expect_equal(rare$pval, fit$pval)
  expect_equal(rare$ci.lb, fit$ci.lb)
  expect_equal(rare$ci.ub, fit$ci.ub)
  # warning message: NAs omitted from model fitting for metafor::rma

  # cc = constant, add 0.5
  # for method = DL
  # for measure = logOR
  # add to "if0all"
  # drop00 = FALSE
  rare <- rareIV(x, measure = "logOR", method = "DL", cc = "constant", ccval = 0.5, ccto = "if0all", drop00 = FALSE, level = 95)
  fit <- metafor::rma(ai = ai, bi = bi, ci = ci, di = di, n1i = n1i, n2i = n2i, data = data, measure = "OR", method = "DL", add = 0.5, to = "if0all", drop00 = FALSE, level = 95)
  expect_equal(rare$beta, fit$beta)
  expect_equal(rare$tau2, fit$tau2)
  expect_equal(rare$se.tau2, fit$se.tau2)
  expect_equal(rare$se, fit$se)
  expect_equal(rare$pval, fit$pval)
  expect_equal(rare$ci.lb, fit$ci.lb)
  expect_equal(rare$ci.ub, fit$ci.ub)

  # cc = constant, add 0.5
  # for method = DL
  # for measure = logOR
  # add to "if0all"
  # drop00 = TRUE
  rare <- rareIV(x, measure = "logOR", method = "DL", cc = "constant", ccval = 0.5, ccto = "if0all", drop00 = TRUE, level = 95)
  fit <- suppressWarnings(metafor::rma(ai = ai, bi = bi, ci = ci, di = di, n1i = n1i, n2i = n2i, data = data, measure = "OR", method = "DL", add = 0.5, to = "if0all", drop00 = TRUE, level = 95))
  expect_equal(rare$beta, fit$beta)
  expect_equal(rare$tau2, fit$tau2)
  expect_equal(rare$se.tau2, fit$se.tau2)
  expect_equal(rare$se, fit$se)
  expect_equal(rare$pval, fit$pval)
  expect_equal(rare$ci.lb, fit$ci.lb)
  expect_equal(rare$ci.ub, fit$ci.ub)
  # warning message: NAs omitted from model fitting for metafor::rma

  ##############################################################################
  ##################### measure = logRR #######################################

  # cc = constant, add 0.5
  # for method = DL
  # for measure = logRR
  # add to "only0"
  # drop00 = FALSE
  rare <- rareIV(x, measure = "logRR", method = "DL", cc = "constant", ccval = 0.5, ccto = "only0", drop00 = FALSE, level = 95)
  fit <- metafor::rma(ai = ai, bi = bi, ci = ci, di = di, n1i = n1i, n2i = n2i, data = data, measure = "RR", method = "DL", add = 0.5, to = "only0", drop00 = FALSE, level = 95)
  expect_equal(rare$beta, fit$beta)
  expect_equal(rare$tau2, fit$tau2)
  expect_equal(rare$se.tau2, fit$se.tau2)
  expect_equal(rare$se, fit$se)
  expect_equal(rare$pval, fit$pval)
  expect_equal(rare$ci.lb, fit$ci.lb)
  expect_equal(rare$ci.ub, fit$ci.ub)

  # cc = constant, add 0.5
  # for method = DL
  # for measure = logRR
  # add to "only0"
  # drop00 = TRUE
  rare <- rareIV(x, measure = "logRR", method = "DL", cc = "constant", ccval = 0.5, ccto = "only0", drop00 = TRUE, level = 95)
  fit <- suppressWarnings(metafor::rma(ai = ai, bi = bi, ci = ci, di = di, n1i = n1i, n2i = n2i, data = data, measure = "RR", method = "DL", add = 0.5, to = "only0", drop00 = TRUE, level = 95))
  expect_equal(rare$beta, fit$beta)
  expect_equal(rare$tau2, fit$tau2)
  expect_equal(rare$se.tau2, fit$se.tau2)
  expect_equal(rare$se, fit$se)
  expect_equal(rare$pval, fit$pval)
  expect_equal(rare$ci.lb, fit$ci.lb)
  expect_equal(rare$ci.ub, fit$ci.ub)
  # warning message: NAs omitted from model fitting for metafor::rma


  # cc = constant, add 0.5
  # for method = DL
  # for measure = logRR
  # add to "all"
  # drop00 = FALSE
  rare <- rareIV(x, measure = "logRR", method = "DL", cc = "constant", ccval = 0.5, ccto = "all", drop00 = FALSE, level = 95)
  fit <- suppressWarnings(metafor::rma(ai = ai, bi = bi, ci = ci, di = di, n1i = n1i, n2i = n2i, data = data, measure = "RR", method = "DL", add = 0.5, to = "all", drop00 = FALSE, level = 95))
  expect_equal(rare$beta, fit$beta)
  expect_equal(rare$tau2, fit$tau2)
  expect_equal(rare$se.tau2, fit$se.tau2)
  expect_equal(rare$se, fit$se)
  expect_equal(rare$pval, fit$pval)
  expect_equal(rare$ci.lb, fit$ci.lb)
  expect_equal(rare$ci.ub, fit$ci.ub)

  # cc = constant, add 0.5
  # for method = DL
  # for measure = logRR
  # add to "all"
  # drop00 = TRUE
  rare <- rareIV(x, measure = "logRR", method = "DL", cc = "constant", ccval = 0.5, ccto = "all", drop00 = TRUE, level = 95)
  fit <- suppressWarnings(metafor::rma(ai = ai, bi = bi, ci = ci, di = di, n1i = n1i, n2i = n2i, data = data, measure = "RR", method = "DL", add = 0.5, to = "all", drop00 = TRUE, level = 95))
  expect_equal(rare$beta, fit$beta)
  expect_equal(rare$tau2, fit$tau2)
  expect_equal(rare$se.tau2, fit$se.tau2)
  expect_equal(rare$se, fit$se)
  expect_equal(rare$pval, fit$pval)
  expect_equal(rare$ci.lb, fit$ci.lb)
  expect_equal(rare$ci.ub, fit$ci.ub)
  # warning message: NAs omitted from model fitting for metafor::rma

  # cc = constant, add 0.5
  # for method = DL
  # for measure = logRR
  # add to "if0all"
  # drop00 = FALSE
  rare <- rareIV(x, measure = "logRR", method = "DL", cc = "constant", ccval = 0.5, ccto = "if0all", drop00 = FALSE, level = 95)
  fit <- metafor::rma(ai = ai, bi = bi, ci = ci, di = di, n1i = n1i, n2i = n2i, data = data, measure = "RR", method = "DL", add = 0.5, to = "if0all", drop00 = FALSE, level = 95)
  expect_equal(rare$beta, fit$beta)
  expect_equal(rare$tau2, fit$tau2)
  expect_equal(rare$se.tau2, fit$se.tau2)
  expect_equal(rare$se, fit$se)
  expect_equal(rare$pval, fit$pval)
  expect_equal(rare$ci.lb, fit$ci.lb)
  expect_equal(rare$ci.ub, fit$ci.ub)

  # cc = constant, add 0.5
  # for method = DL
  # for measure = logRR
  # add to "if0all"
  # drop00 = TRUE
  rare <- rareIV(x, measure = "logRR", method = "DL", cc = "constant", ccval = 0.5, ccto = "if0all", drop00 = TRUE, level = 95)
  fit <- suppressWarnings(metafor::rma(ai = ai, bi = bi, ci = ci, di = di, n1i = n1i, n2i = n2i, data = data, measure = "RR", method = "DL", add = 0.5, to = "if0all", drop00 = TRUE, level = 95))
  expect_equal(rare$beta, fit$beta)
  expect_equal(rare$tau2, fit$tau2)
  expect_equal(rare$se.tau2, fit$se.tau2)
  expect_equal(rare$se, fit$se)
  expect_equal(rare$pval, fit$pval)
  expect_equal(rare$ci.lb, fit$ci.lb)
  expect_equal(rare$ci.ub, fit$ci.ub)
  # warning message: NAs omitted from model fitting for metafor::rma
})



# test if studies to be continuity corrected are specified correctly
# and test if continuity correction is applied correctly
# test_that("rareIV correctly specifies studies to be continuity corrected and applies cc correctly", {
# })

# test_that("Implementation of method = 'IPM' works as it should", {
# })

# test_that("Some basic results from rareIV align with those from metafor.", {
# })
