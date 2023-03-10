##Tests for the raremeta function rareMH()
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

test_that("rareMH runs with valid inputs", {

  expect_error(
    rareMH(x, measure = "logOR"),
    NA
  )

  expect_error(
    rareMH(x, measure = "logRR"),
    NA
  )

  expect_error(
    rareMH(x, measure = "RD"),
    NA
  )

  expect_error(
    rareMH(x, measure = "logOR", level = 0),
    NA
  )

  expect_error(
    rareMH(x, measure = "logOR", level = 50),
    NA
  )

  expect_error(
    rareMH(x, measure = "logOR", level = 100),
    NA
  )

  expect_error(
    rareMH(x, measure = "logRR", level = 0, digits = 1),
    NA
  )

  expect_error(
    rareMH(x, measure = "logRR", level = 50, digits = 2),
    NA
  )

  expect_error(
    rareMH(x, measure = "logRR", level = 100, digits = 3),
    NA
  )
})

test_that("rareMH returns errors and warning messages", {

  # testing the data input
  expect_error(
    rareMH(x = data, measure = "logOR"),
    "x must be an object of class 'rareData'."
    )

  #testing the measure input
  expect_error(
    rareMH(x),
    "'measure' argument must be specified."
  )

  expect_error(
    rareMH(x, measure = "OR"),
    "'measure' must be either 'logOR', 'logRR', or 'RD'."
  )

  expect_error(
    rareMH(x, measure = "RR"),
    "'measure' must be either 'logOR', 'logRR', or 'RD'."
  )

  expect_error(
    rareMH(x, measure = "logRD"),
    "'measure' must be either 'logOR', 'logRR', or 'RD'."
  )

  #testing the level input
  expect_error(
    rareMH(x, measure = "logOR", level = -1),
    "level must be a scalar between 0 and 100."
  )

  expect_error(
    rareMH(x, measure = "logRR", level = c(1,2)),
    "level must be a scalar between 0 and 100."
  )

  expect_error(
    rareMH(x, measure = "RD", level = "ninetyfive"),
    "level must be a scalar between 0 and 100."
  )

  expect_error(
    rareMH(x, measure = "logOR", level = 101),
    "level must be a scalar between 0 and 100."
  )

  expect_error(
    rareMH(x, measure = "logOR", digits = -1),
    "'digits' must be an integer of length 1."
  )

  expect_error(
    rareMH(x, measure = "logRR", digits = 1.6),
    "'digits' must be an integer of length 1."
  )

  expect_error(
    rareMH(x, measure = "RD", digits = c(1,2)),
    "'digits' must be an integer of length 1."
  )

})

test_that("comparing output of rareMH to output of metafor::rma.mh",{

  #log(OR) with defaults (level=95, digits=4)
  rare    <- rareMH(x, measure = "logOR")
  fit_rma <- metafor::rma.mh(ai=ai, bi=bi, ci=ci, di=di, data=data, measure="OR", drop00 = FALSE)

  expect_equal(rare$b, fit_rma$b, ignore_attr = TRUE)
  expect_equal(rare$beta, fit_rma$beta, ignore_attr = TRUE)
  expect_equal(rare$se, fit_rma$se, ignore_attr = TRUE)
  expect_equal(rare$zval, fit_rma$zval, ignore_attr = TRUE)
  expect_equal(rare$pval, fit_rma$pval, ignore_attr = TRUE)
  expect_equal(rare$ci.lb, fit_rma$ci.lb, ignore_attr = TRUE)
  expect_equal(rare$ci.ub, fit_rma$ci.ub, ignore_attr = TRUE)
  expect_equal(rare$k, fit_rma$k, ignore_attr = TRUE)

  #log(OR) with level=50, digits = 5
  rare    <- rareMH(x, measure = "logOR", level = 50, digits = 5)
  fit_rma <- rma.mh(ai=ai, bi=bi, ci=ci, di=di, data=data, measure="OR",
                    level = 50, digits = 5, drop00 = FALSE)

  expect_equal(rare$b, fit_rma$b, ignore_attr = TRUE)
  expect_equal(rare$beta, fit_rma$beta, ignore_attr = TRUE)
  expect_equal(rare$se, fit_rma$se, ignore_attr = TRUE)
  expect_equal(rare$zval, fit_rma$zval, ignore_attr = TRUE)
  expect_equal(rare$pval, fit_rma$pval, ignore_attr = TRUE)
  expect_equal(rare$ci.lb, fit_rma$ci.lb, ignore_attr = TRUE)
  expect_equal(rare$ci.ub, fit_rma$ci.ub, ignore_attr = TRUE)
  expect_equal(rare$k, fit_rma$k, ignore_attr = TRUE)

  #log(RR) with defaults (level=95, digits=4)
  rare    <- rareMH(x, measure = "logRR")
  fit_rma <- metafor::rma.mh(ai=ai, bi=bi, ci=ci, di=di, data=data, measure="RR", drop00 = FALSE)

  expect_equal(rare$b, fit_rma$b, ignore_attr = TRUE)
  expect_equal(rare$beta, fit_rma$beta, ignore_attr = TRUE)
  expect_equal(rare$se, fit_rma$se, ignore_attr = TRUE)
  expect_equal(rare$zval, fit_rma$zval, ignore_attr = TRUE)
  expect_equal(rare$pval, fit_rma$pval, ignore_attr = TRUE)
  expect_equal(rare$ci.lb, fit_rma$ci.lb, ignore_attr = TRUE)
  expect_equal(rare$ci.ub, fit_rma$ci.ub, ignore_attr = TRUE)
  expect_equal(rare$k, fit_rma$k, ignore_attr = TRUE)

  #log(OR) with level=50, digits = 5
  rare    <- rareMH(x, measure = "logRR", level = 50, digits = 5)
  fit_rma <- metafor::rma.mh(ai=ai, bi=bi, ci=ci, di=di, data=data, measure="RR",
                    level = 50, digits = 5, drop00 = FALSE)

  expect_equal(rare$b, fit_rma$b, ignore_attr = TRUE)
  expect_equal(rare$beta, fit_rma$beta, ignore_attr = TRUE)
  expect_equal(rare$se, fit_rma$se, ignore_attr = TRUE)
  expect_equal(rare$zval, fit_rma$zval, ignore_attr = TRUE)
  expect_equal(rare$pval, fit_rma$pval, ignore_attr = TRUE)
  expect_equal(rare$ci.lb, fit_rma$ci.lb, ignore_attr = TRUE)
  expect_equal(rare$ci.ub, fit_rma$ci.ub, ignore_attr = TRUE)
  expect_equal(rare$k, fit_rma$k, ignore_attr = TRUE)

  #RD with defaults (level=95, digits=4)
  rare    <- rareMH(x, measure = "RD")
  fit_rma <- metafor::rma.mh(ai=ai, bi=bi, ci=ci, di=di, data=data, measure="RD", drop00 = FALSE)

  expect_equal(rare$b, fit_rma$b, ignore_attr = TRUE)
  expect_equal(rare$beta, fit_rma$beta, ignore_attr = TRUE)
  expect_equal(rare$se, fit_rma$se, ignore_attr = TRUE)
  expect_equal(rare$zval, fit_rma$zval, ignore_attr = TRUE)
  expect_equal(rare$pval, fit_rma$pval, ignore_attr = TRUE)
  expect_equal(rare$ci.lb, fit_rma$ci.lb, ignore_attr = TRUE)
  expect_equal(rare$ci.ub, fit_rma$ci.ub, ignore_attr = TRUE)
  expect_equal(rare$k, fit_rma$k, ignore_attr = TRUE)

  #log(OR) with level=50, digits = 5
  rare    <- rareMH(x, measure = "RD", level = 50, digits = 5)
  fit_rma <- metafor::rma.mh(ai=ai, bi=bi, ci=ci, di=di, data=data, measure="RD",
                        level = 50, digits = 5, drop00 = FALSE)

  expect_equal(rare$b, fit_rma$b, ignore_attr = TRUE)
  expect_equal(rare$beta, fit_rma$beta, ignore_attr = TRUE)
  expect_equal(rare$se, fit_rma$se, ignore_attr = TRUE)
  expect_equal(rare$zval, fit_rma$zval, ignore_attr = TRUE)
  expect_equal(rare$pval, fit_rma$pval, ignore_attr = TRUE)
  expect_equal(rare$ci.lb, fit_rma$ci.lb, ignore_attr = TRUE)
  expect_equal(rare$ci.ub, fit_rma$ci.ub, ignore_attr = TRUE)
  expect_equal(rare$k, fit_rma$k, ignore_attr = TRUE)






})
