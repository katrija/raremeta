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


  ##empirical continuity correction
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

#test_that("rareMH returns errors and warning messages", {

#  expect_error()

#  })
