# usual data set:
data <- data.frame(
  ai = c(0, 3, 2, 0, 5, 4, 3),
  bi = c(200, 180, 150, 190, 170, 185, 205),
  ci = c(1, 4, 0, 0, 3, 4, 5),
  di = c(190, 170, 160, 200, 180, 180, 195)
)

data$n1i <- data$ai + data$bi
data$n2i <- data$ci + data$di

# long format (for aods3):
data_long <- data.frame(
  xi = c(data$ai, data$ci),
  ni = c(data$n1i, data$n2i),
  group = rep(1:0, each = nrow(data))
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
data_no00 <- data[!(data$ai == 0 & data$ci == 0), ]

x_no00 <- rareDescribe(
  ai = ai,
  bi = bi,
  ci = ci,
  di = di,
  n1i = n1i,
  n2i = n2i,
  data = data_no00
)

##### VALID INPUTS -------------------------------------------------------------
test_that("rareBetabin runs with valid inputs", {
  expect_error(
    rareBetabin(x, measure = "logOR", common_rho = TRUE, drop00 = FALSE),
    NA
  )

  expect_error(
    rareBetabin(x, measure = "logRR", common_rho = TRUE, drop00 = FALSE),
    NA
  )


})

##### ERRORS AND WARNINGS ------------------------------------------------------
test_that("rareBetabin returns errors and warning messages", {
  # error if x is not object of class rareData
  expect_error(
    rareBetabin(x = data, measure = "logOR"),
    "x must be an object of class 'rareData'."
  )

  expect_error(
    rareBetabin(x = x),
    "'measure' argument must be specified."
  )
})

##### EXCLUSION OF DOUBLE-ZERO STUDIES -----------------------------------------
fit1 <- rareBetabin(x, measure = "logOR", drop00 = TRUE)
fit1_no00 <- rareBetabin(x_no00, measure = "logOR", drop00 = FALSE)

test_that("rareBetabin handles the exclusion of double-zero studies correctly", {
  expect_equal(fit1$beta, fit1_no00$beta)
  expect_equal(fit1$ci.lb, fit1_no00$ci.lb)
  expect_equal(fit1$ci.ub, fit1_no00$ci.ub)
  expect_equal(fit1$tau2, fit1_no00$tau2)
})

##### COMPARISON TO ALTERNATIVE IMPLEMENTATIONS --------------------------------

# log OR with common rho:
fit <- rareBetabin(x, measure = "logOR", drop00 = FALSE)
fit_aods <- suppressWarnings(aods3::aodml(cbind(xi, ni-xi)~1+group, data = data_long,
                         link = "logit", family = "bb"))

test_that("rareBetabin yields similar results as aodml", {
  expect_equal(fit$beta, fit_aod$b, tolerance = 0.01, ignore_attr = TRUE)
  expect_equal(fit$se, sqrt(diag(fit_aod$varparam)[1:2]), tolerance = 0.01, ignore_attr = TRUE)
  expect_equal(fit$rho, fit_aod$phi/(1+fit_aod$phi), tolerance = 0.01, ignore_attr = TRUE)
})

# log OR with different rhos:
# currently not testable: the following code yields an error:
# fit_aods <- aods3::aodml(cbind(xi, ni-xi)~1+group, data = data_long,
# phi.formula = ~1+group, link = "logit", family = "bb")
# (also, choosing another optimizer does not seem to help)

# log RR:
# currently not testable (no alternative implementations as far as I am aware)

# look out for future developments!

##### OUTPUT STRUCTURE ---------------------------------------------------------
test_that("rareBetabin returns the appropriate outputs", {
  expect_equal(length(fit1$beta), 2)
  expect_equal(length(rareBetabin(x, measure = "logOR", common_rho = TRUE)$rho), 1)
  expect_equal(length(rareBetabin(x, measure = "logOR", common_rho = FALSE)$rho), 2)
})

##### CHECKING EQUIVALENT WAYS OF DATA INPUT -----------------------------------
test_that("does not matter if data is put in as data frame or rareData-object",{

  k <- which(names(rareBetabin(x=x, measure = "logOR")) == "call")

  expect_identical(
    rareBetabin(x=x, measure = "logOR")[-k],
    rareBetabin(ai=ai, bi=bi, ci=ci, di=di, n1i=n1i, n2i=n2i, data=data, measure = "logOR")[-k],
  )

  expect_identical(
    rareBetabin(x=x, measure = "logOR")[-k],
    rareBetabin(ai=ai, ci=ci, n1i=n1i, n2i=n2i, data=data, measure = "logOR")[-k]
  )

  expect_identical(
    rareBetabin(x=x, measure = "logOR")[-k],
    rareBetabin(ai=ai, ci=ci, di=di, n1i=n1i, data=data, measure = "logOR")[-k]
  )

  expect_identical(
    rareBetabin(x=x, measure = "logOR")[-k],
    rareBetabin(ai=ai, bi=bi, ci=ci, di=di, data=data, measure = "logOR")[-k]
  )

})

