# usual data set:
data <- data.frame(
  ai = c(0,3,2,0),
  bi = c(20,18,15,19),
  ci = c(1,4,0,0),
  di = c(19,17,16,20),
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
  data = data)

# data set with only zero-studies
data_only0 <- data.frame(
  ai = c(1,0,0,1),
  bi = c(0,2,0,0),
  ci = c(0,4,1,0),
  di = c(19,0,0,20),
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
  data = data_only0)

test_that("rareIV returns errors and warning messages", {

  # error if x is not object of class rareData
  expect_error(rareIV(x = data, measure = "logOR", method = "FE", cc = "constant"),
               "x must be an object of class 'rareData'.")

  # error if measure argument is not defined
  expect_error(rareIV(x, method = "FE", cc = "constant"),
               "'measure' argument must be specified.")

  # error if measure argument is not valid
  expect_error(rareIV(x, measure = "logoR", method = "FE", cc = "constant"),
               "'measure' must be either 'logOR', 'logRR', or 'RD'.")

  # error if method argument is not defined
  expect_error(rareIV(x, measure = "logOR", cc = "constant"),
               "'method' argument must be specified.")

  # error if method argument is not valid
  expect_error(rareIV(x, measure = "logOR", method = "fe", cc = "constant"),
               "unknown 'method' specified.")

  # error if cc is not specified
  expect_error(rareIV(x, measure = "logOR", method = "FE"),
               "Some studies have zero events.")

  # error if drop00 is not logical
  expect_error(rareIV(x, measure = "logOR", method = "FE", cc = "constant", drop00 = "true"),
               "'drop00' must be a logical")

  # error if cc argument is not valid
  expect_error(rareIV(x, measure = "logOR", method = "FE", cc = "None"),
               "'cc' must be either 'none', 'constant', 'reciprocal', or 'empirical'")

  # error if ccto argument is not valid
  expect_error(rareIV(x, measure = "logOR", method = "FE", cc = "constant", ccto = "Only0"),
               "'ccto' must be either 'only0', 'all', or 'if0all'")

  # error if ccval argument is not valid
  ccval = c(1,1)

  #for drop00 = FALSE
  expect_error(rareIV(x, measure = "logOR", method = "FE", cc = "constant", drop00 = FALSE, ccval = ccval),
               "'ccval' must have length 1 or length equal to the number of studies.")

  # for drop00 = TRUE
  expect_error(rareIV(x, measure = "logOR", method = "FE", cc = "constant", drop00 = TRUE, ccval = ccval),
               "'ccval' must have length 1 or length equal to the number of studies")

  # error if test argument is not valid
  expect_error(rareIV(x, measure = "logOR", method = "FE", cc = "constant", test = "Z"),
               "'test' must be either 'z', 'knha', or 'hksj'")

  # error if digits argument is not valid
  expect_error(rareIV(x, measure = "logOR", method = "FE", cc = "constant", digits = 0.5),
               "'digits' must be an integer of length 1.")

  # error if ccval argument is negative
  expect_error(rareIV(x, measure = "logOR", method = "FE", cc = "constant", ccval = -1),
               "'ccval' must be non-negative.")

  # error if cc = empirical & method = RD
  expect_error(rareIV(x, measure = "RD", method = "FE", cc = "empirical"),
               "continuity correction of type 'empirical' is currently not supported for measure 'RD'.")

  # error if there are no non-zero studies and cc = empirical
  expect_error(rareIV(x_only0, measure = "logOR", method = "FE", cc = "empirical"),
               "continuity correction of type 'empirical' can only be applied if there is at least one non-zero study.")

  # error if cc = empirical & method = IPM
  expect_error(rareIV(x, measure = "logOR", method = "IPM", cc = "empirical"),
               "method = 'IPM' is currently not defined for cc = 'empirical'.")

  # error if method = IPM & measure != logOR
  expect_error(rareIV(x, measure = "logRR", method = "IPM", cc = "constant"),
               "method = 'IPM' is only defined for measure = 'logOR'.")

  # warning if method = IPM & cc = constant & ccto = only0
  expect_warning(rareIV(x, measure = "logOR", method = "IPM", cc = "constant", ccto = "only0"),
                  "method = 'IPM' was developed for cc = 'constant' and ccto = 'all'.")

  }
)

# test if studies to be continuity corrected are specified correctly
# and test if continuity correction is applied correctly
# test_that("rareIV correctly specifies studies to be continuity corrected and applies cc correctly", {
# })

# test_that("Implementation of method = 'IPM' works as it should", {
# })

# test_that("Some basic results from rareIV align with those from metafor.", {
# })
