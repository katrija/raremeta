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
    "x must be an object of class 'rareData'. See ?rareDescribe for more details."
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
  expect_equal(rareMH(x,"logOR"))





})
