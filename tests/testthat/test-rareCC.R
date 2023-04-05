##Tests for the raremeta- function rareCC()

#rareCC <- function(x, cc = "constant", ccval = 0.5, tccval, cccval, ccsum = 1,
#                   ccto = "only0", drop00 = TRUE, measure, method = "FE"){

#usual data frame including double/single/no-zero studies
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


test_that("rareCC runs with valid inputs",{
  # defaults
  expect_error(
    rareCC(x),
    NA
  )

  ##constant continuity correction

  #ccval/cccval+tccval
  expect_error(
    rareCC(x, ccval = 1.5),
    NA
  )

  expect_error(
    rareCC(x, ccval = c(0.5, 1, 1.5, 2)),
    NA
  )

  expect_error(
    rareCC(x, cccval = 1, tccval = 1.5),
    NA
  )

  expect_error(
    rareCC(x, cccval = c(0.5, 1, 1.5 , 2), tccval = c(1, 1.5 , 2, 2.5)),
    NA
  )

  #ccto
  expect_error(
    rareCC(x, ccto = "if0all"),
    NA
  )

  expect_error(
    rareCC(x, ccto = "all"),
    NA
  )

  expect_error(
    rareCC(x, ccto = "if0all", drop00 = FALSE),
    NA
  )


  ## treatment arm continuity correction
  expect_error(
    rareCC(x, cc = "tacc", measure = "logOR"),
    NA
  )

  expect_error(
    rareCC(x, cc = "tacc", measure = "logRR"),
    NA
  )

  expect_error(
    rareCC(x, cc = "tacc", measure = "logOR", ccsum = 2),
    NA
  )

  expect_error(
    rareCC(x, method = "tacc", measure = "logRR", ccsum = 2),
    NA
  )

  expect_error(
    rareCC(x, cc = "tacc", measure = "logOR", drop00 = FALSE),
    NA
  )

  expect_error(
    rareCC(x, cc = "tacc", measure = "logRR", drop00 = FALSE),
    NA
  )

  expect_error(
    rareCC(x, cc = "tacc", measure = "logOR", method = "DL"),
    NA
  )

  expect_error(
    rareCC(x, cc = "tacc", measure = "logRR", method = "DL"),
    NA
  )

  expect_error(
    rareCC(x, cc = "tacc", measure = "logOR", method = "DL", drop00 = FALSE),
    NA
  )

  expect_error(
    rareCC(x, cc = "tacc", measure = "logRR", method = "DL", drop00 = FALSE),
    NA
  )


  ## empirical continuity correction
  expect_error(
    rareCC(x, cc = "empirical", measure = "logOR"),
    NA
  )

  expect_error(
    rareCC(x, cc = "empirical", measure = "logRR"),
    NA
  )

  expect_error(
    rareCC(x, cc = "empirical", measure = "logOR", ccsum = 2),
    NA
  )

  expect_error(
    rareCC(x, cc = "empirical", measure = "logRR", ccsum = 2),
    NA
  )

  expect_error(
    rareCC(x, cc = "empirical", measure = "logOR", method = "REML"),
    NA
  )

  expect_error(
    rareCC(x, cc = "empirical", measure = "logRR", method = "REML"),
    NA
  )

  expect_error(
    rareCC(x, cc = "empirical", measure = "logOR", method = "SJ", drop00 = FALSE),
    NA
  )

  expect_error(
    rareCC(x, cc = "empirical", measure = "logRR", method = "SJ", drop00 = FALSE),
    NA
  )

})

## test warning messages

test_that("rareMH returns errors and warning messages", {

  expect_error(
    rareCC(x, cc="something"),
    "'cc' must be either 'none', 'constant', 'tacc', or 'empirical'."
  )

  expect_error(
    rareCC(x, cc = 0.5),
    "'cc' must be either 'none', 'constant', 'tacc', or 'empirical'."
  )

  expect_error(
    rareCC(x, drop00 = "yes"),
    "'drop00' must be a logical."
  )

  expect_error(
    rareCC(x, cc= "constant", ccto = "someofthem"),
    "'ccto' must be either 'only0', 'all', 'none', or 'if0all'."
  )

  #ccval, cccval, tccval

  expect_error(
    rareCC(x, cc="constant", drop00 = FALSE, ccval = c(1,2,3)),
    "'ccval' must have length 1 or length equal to the number of studies."
  )

  expect_error(
    rareCC(x, cc="constant", drop00 = TRUE, ccval = c(1,2)),
    "'ccval' must have length 1 or length equal to the number of studies"
  )

  expect_error(
    rareCC(x, cc="constant", ccval = -1),
    "'ccval' must be non-negative."
  )

  expect_error(
    rareCC(x, cc="constant", ccval = c(1,1,-1,1)),
    "'ccval' must be non-negative."
  )

  expect_error(
    rareCC(x, cc = "constant", tccval = 0.5),
    "Please specify both 'tccval' and 'cccval'."
  )

  expect_error(
    rareCC(x, cc = "constant", cccval = 0.5),
    "Please specify both 'tccval' and 'cccval'."
  )

  expect_error(
    rareCC(x, cc = "constant", cccval = c(1,1,1,1), tccval = c(1,1,1)),
    "'tccval' and 'cccval' must have equal length."
  )

  expect_error(
    rareCC(x, cc = "constant", drop00 = FALSE, cccval = c(1,1,1), tccval = c(1,1,1)),
    "'tccval' must have length 1 or length equal to the number of studies."
  )

  expect_error(
    rareCC(x, cc = "constant", drop00 = TRUE, cccval = c(1,1), tccval = c(1,1)),
    "'tccval' must have length 1 or length equal to the number of studies"
  )

  expect_error(
    rareCC(x, cc = "constant", drop00 = FALSE, cccval = c(1,1,1,-1), tccval = c(1,1,1,1)),
    "All values in 'tccval' and 'cccval' must be non-negative."
  )

  expect_error(
    rareCC(x, cc = "constant", drop00 = FALSE, cccval = c(1,1,1,1), tccval = c(1,1,1,-1)),
    "All values in 'tccval' and 'cccval' must be non-negative."
  )

  expect_error(
    rareCC(x, cc = "empirical"),
    "To apply the the empirical continuity correction, the 'measure' argument must be specified."
  )

  # 'empirical' and 'tacc not yet supported for RD
  expect_error(
    rareCC(x, cc = "tacc", measure = "RD"),
    "continuity correction of type 'tacc' is currently not supported for measure 'RD'."
  )

  expect_error(
    rareCC(x, cc = "empirical", measure = "RD"),
    "continuity correction of type 'empirical' is currently not supported for measure 'RD'."
  )

  expect_error(
    rareCC(x, cc = "empirical", measure = "Lebesgue"),
    "To apply the empirical continuity correction, 'measure' must be either 'logOR', 'logRR' or 'RD'."
  )

  expect_error(
    rareCC(x, cc = "empirical", measure = "logOR", method = "scientific")
  )

  expect_error(
    rareCC(x, cc = "tacc", measure = "logOR", ccsum = -1),
    "ccsum must be a scalar larger than 0."
  )

  expect_error(
    rareCC(x, cc = "tacc", measure = "logOR", ccsum = c(1,1)),
    "ccsum must be a scalar larger than 0."
  )

  expect_error(
    rareCC(x, cc = "tacc", measure = "logOR", ccsum = "one and a half"),
    "ccsum must be a scalar larger than 0."
  )


  })


#testing continuity correction vs correcting by hand for cc = "constant" with ccval = 1

#ccto  = "all"
dataCC_all <- data.frame(
  ai = c(1, 4, 3, 1),
  bi = c(21, 19, 16, 20),
  ci = c(2, 5, 1, 1),
  di = c(20, 18, 17, 21),
  n1i = c(22, 23, 19, 21),
  n2i = c(22, 23, 18, 22)
)

x_all <- rareDescribe(
  ai = ai,
  bi = bi,
  ci = ci,
  di = di,
  n1i = n1i,
  n2i = n2i,
  data = dataCC_all
)

#only 0
dataCC_only0 <- data.frame(
  ai = c(1, 3, 3, 1),
  bi = c(21, 18, 16, 20),
  ci = c(2, 4, 1, 1),
  di = c(20, 17, 17, 21),
  n1i = c(22, 21, 19, 21),
  n2i = c(22, 21, 18, 22)
)

x_only0 <- rareDescribe(
  ai = ai,
  bi = bi,
  ci = ci,
  di = di,
  n1i = n1i,
  n2i = n2i,
  data = dataCC_only0
)

test_that("comparing constant CC to result by hand",{

  #ccto = "all"
  x.rareCC  <- rareCC(x, cc = "constant", ccval = 1, ccto = "all", drop00 = FALSE)

  expect_equal(
    x.rareCC$ai.cc,
    x_all$ai
  )

  expect_equal(
    x.rareCC$bi.cc,
    x_all$bi
  )

  expect_equal(
    x.rareCC$ci.cc,
    x_all$ci
  )

  expect_equal(
    x.rareCC$di.cc,
    x_all$di
  )

  expect_equal(
    x.rareCC$n1i.cc,
    x_all$n1i
  )

  expect_equal(
    x.rareCC$n2i.cc,
    x_all$n2i
  )

  #ccto = "only0"
  x.rareCC  <- rareCC(x, cc = "constant", ccval = 1, ccto = "only0", drop00 = FALSE)

  expect_equal(
    x.rareCC$ai.cc,
    x_only0$ai
  )

  expect_equal(
    x.rareCC$bi.cc,
    x_only0$bi
  )

  expect_equal(
    x.rareCC$ci.cc,
    x_only0$ci
  )

  expect_equal(
    x.rareCC$di.cc,
    x_only0$di
  )

  expect_equal(
    x.rareCC$n1i.cc,
    x_only0$n1i
  )

  expect_equal(
    x.rareCC$n2i.cc,
    x_only0$n2i
  )
})

