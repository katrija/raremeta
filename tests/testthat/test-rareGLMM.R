# usual data set:
data <- data.frame(
  ai = c(0, 3, 2, 0),
  bi = c(200, 180, 150, 190),
  ci = c(1, 4, 0, 0),
  di = c(190, 170, 160, 200)
)

data$n1i <- data$ai+data$bi
data$n2i <- data$ci+data$di

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
data_no00 <- data[!(data$ai == 0 & data$ci == 0),]

x_no00 <- rareDescribe(
  ai = ai,
  bi = bi,
  ci = ci,
  di = di,
  n1i = n1i,
  n2i = n2i,
  data = data_no00
)

##### ERRORS AND WARNINGS ------------------------------------------------------

test_that("rareGLMM returns errors and warning messages", {

  # error if x is not object of class rareData
  expect_error(
    rareGLMM(x = data, measure = "logOR", intercept = "random", slope = "fixed"),
    "x must be an object of class 'rareData'."
  )

  expect_error(
    rareGLMM(x = x, intercept = "random", slope = "fixed"),
    "'measure' argument must be specified."
  )



})

##### EXCLUSION OF DOUBLE-ZERO STUDIES -----------------------------------------

fit1 <- rareGLMM(x, measure = "logOR", intercept = "random", slope = "random",
                 drop00 = TRUE)
fit1_no00 <- rareGLMM(x_no00, measure = "logOR", intercept = "random", slope = "random",
                      drop00 = TRUE)

test_that("rareGLMM handles the exclusion of double-zero studies correctly", {

  expect_equal(fit1$beta, fit1_no00$beta)
  expect_equal(fit1$ci.lb, fit1_no00$ci.lb)
  expect_equal(fit1$ci.ub, fit1_no00$ci.ub)
  expect_equal(fit1$tau2, fit1_no00$tau2)

})


##### OUTPUT STRUCTURE ---------------------------------------------------------

fit_rifs <- rareGLMM(x, measure = "logOR", intercept = "random", slope = "fixed",
                     drop00 = TRUE)
fit_rirs <- rareGLMM(x, measure = "logOR", intercept = "random", slope = "random",
                     drop00 = TRUE)

test_that("rareGLMM returns the appropriate outputs", {

  expect_equal(fit_rifs$tau2, 0)
  expect_equal(length(fit_rirs$tau2), 1)

})
