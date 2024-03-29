# usual data set:
data <- data.frame(
  ai = c(0, 3, 2, 0, 5, 4, 3),
  bi = c(200, 180, 150, 190, 170, 185, 205),
  ci = c(1, 4, 0, 0, 3, 4, 5),
  di = c(190, 170, 160, 200, 180, 180, 195)
)

data$n1i <- data$ai + data$bi
data$n2i <- data$ci + data$di

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
test_that("rareGLMM runs with valid inputs", {
  expect_error(
    suppressMessages(rareGLMM(x,
      measure = "logOR",
      slope = "random", intercept = "random"
    )),
    NA
  )

  expect_error(
    suppressWarnings(suppressMessages(rareGLMM(x,
      measure = "logOR",
      slope = "random", intercept = "fixed"
    ))),
    NA
  )

  expect_error(
    rareGLMM(x,
      measure = "logOR",
      slope = "fixed", intercept = "fixed"
    ),
    NA
  )

  expect_error(
    suppressMessages(rareGLMM(x,
                              measure = "logRR",
                              slope = "random", intercept = "random"
    )),
    NA
  )

  expect_error(
    suppressWarnings(suppressMessages(rareGLMM(x,
                                               measure = "logRR",
                                               slope = "random", intercept = "fixed"
    ))),
    NA
  )

  expect_error(
    rareGLMM(x,
             measure = "logRR",
             slope = "fixed", intercept = "fixed"
    ),
    NA
  )

  expect_error(
    suppressMessages(rareGLMM(x,
                              measure = "logOR",
                              conditional = TRUE, approx = TRUE
    )),
    NA
  )

  expect_error(
    suppressMessages(rareGLMM(x,
                              measure = "logRR",
                              conditional = TRUE, approx = FALSE
    )),
    NA
  )
})

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
fit1 <- rareGLMM(x,
  measure = "logOR", intercept = "random", slope = "random",
  drop00 = TRUE
)
fit1_no00 <- rareGLMM(x_no00,
  measure = "logOR", intercept = "random", slope = "random"
)

test_that("rareGLMM handles the exclusion of double-zero studies correctly", {
  expect_equal(fit1$beta, fit1_no00$beta)
  expect_equal(fit1$ci.lb, fit1_no00$ci.lb)
  expect_equal(fit1$ci.ub, fit1_no00$ci.ub)
  expect_equal(fit1$tau2, fit1_no00$tau2)
})

### LOG OR ---------------------------------------------------------------------

##### OUTPUT STRUCTURE ---------------------------------------------------------
fit_rifs <- rareGLMM(x, drop00 = FALSE,
  measure = "logOR", intercept = "random", slope = "fixed"
)
fit_rirs <- rareGLMM(x, drop00 = FALSE,
  measure = "logOR", intercept = "random", slope = "random"
)
fit_firs <- suppressWarnings(rareGLMM(x, drop00 = FALSE,
                     measure = "logOR", intercept = "fixed", slope = "random"
))
fit_cond <- suppressWarnings(rareGLMM(x, measure = "logOR", intercept = "random",
                                      conditional = TRUE, approx = TRUE))

test_that("rareGLMM returns the appropriate outputs", {
  expect_equal(fit_rifs$tau2, 0)
  expect_equal(length(fit_rirs$tau2), 1)
  expect_equal(length(fit_rifs$beta), 2)
  expect_equal(length(fit_firs$beta), nrow(data)+1)
})

##### COMPARE TO METAFOR -------------------------------------------------------
fit_rirs_rma <- suppressWarnings(metafor::rma.glmm(ai = ai, bi = bi, ci = ci, di = di,
                                  data = data, measure = "OR", drop00 = FALSE,
                                  model = "UM.RS"))

fit_firs_rma <- suppressWarnings(metafor::rma.glmm(ai = ai, bi = bi, ci = ci, di = di,
                                  data = data, measure = "OR", drop00 = FALSE,
                                  model = "UM.FS"))

fit_cond_rma <- suppressWarnings(metafor::rma.glmm(ai = ai, bi = bi, ci = ci, di = di,
                                                   data = data, measure = "OR", drop00 = TRUE,
                                                   model = "CM.AL"))
# setting drop00 = TRUE ensures that the df of the LRT are the same as in rareGLMM

test_that("compare rareGLMM to metafor",{

  expect_equal(fit_rirs$beta[2], fit_rirs_rma$beta, ignore_attr = TRUE,
               tolerance = 1e-5)
  expect_equal(fit_rirs$se[2], fit_rirs_rma$se, ignore_attr = TRUE,
               tolerance = 1e-5)
  expect_equal(fit_rirs$tau2, fit_rirs_rma$tau2, ignore_attr = TRUE,
               tolerance = 1e-5)
  expect_equal(fit_rirs$LRT.Chisq, fit_rirs_rma$QE.LRT, ignore_attr = TRUE,
               tolerance = 1e-5)
  expect_equal(fit_rirs$LRT.df, fit_rirs_rma$QE.df, ignore_attr = TRUE)

  expect_equal(fit_firs$beta["group"], fit_firs_rma$beta, ignore_attr = TRUE,
               tolerance = 1e-5)
  expect_equal(fit_firs$se["group"], fit_firs_rma$se, ignore_attr = TRUE,
               tolerance = 1e-5)
  expect_equal(fit_firs$tau2, fit_firs_rma$tau2, ignore_attr = TRUE,
               tolerance = 1e-5)
  expect_equal(fit_firs$LRT.Chisq, fit_firs_rma$QE.LRT, ignore_attr = TRUE,
               tolerance = 1e-5)
  expect_equal(fit_rirs$LRT.df, fit_rirs_rma$QE.df, ignore_attr = TRUE)

  expect_equal(fit_cond$beta, fit_cond_rma$beta, ignore_attr = TRUE,
               tolerance = 1e-5)
  expect_equal(fit_cond$se, fit_cond_rma$se, ignore_attr = TRUE,
               tolerance = 1e-5)
  expect_equal(fit_cond$tau2, fit_cond_rma$tau2, ignore_attr = TRUE,
               tolerance = 1e-5)
  expect_equal(fit_cond$LRT.Chisq, fit_cond_rma$QE.LRT, ignore_attr = TRUE,
               tolerance = 1e-5)
  expect_equal(fit_cond$LRT.df, fit_cond_rma$QE.df, ignore_attr = TRUE)

})

### LOG RR ---------------------------------------------------------------------

##### OUTPUT STRUCTURE ---------------------------------------------------------
fit_rifs <- rareGLMM(x, drop00 = FALSE,
                     measure = "logRR", intercept = "random", slope = "fixed"
)
fit_rirs <- rareGLMM(x, drop00 = FALSE,
                     measure = "logRR", intercept = "random", slope = "random"
)
fit_firs <- suppressWarnings(rareGLMM(x, drop00 = FALSE,
                                      measure = "logRR", intercept = "fixed", slope = "random"
))

test_that("rareGLMM returns the appropriate outputs", {
  expect_equal(fit_rifs$tau2, 0)
  expect_equal(length(fit_rirs$tau2), 1)
  expect_equal(length(fit_rifs$beta), 2)
  expect_equal(length(fit_firs$beta), nrow(data)+1)
})

##### COMPARE TO METAFOR -------------------------------------------------------
fit_rirs_rma <- suppressWarnings(metafor::rma.glmm(x1i = ai, t1i = n1i, x2i = ci, t2i = n2i,
                                                   data = data, measure = "IRR", drop00 = FALSE,
                                                   model = "UM.RS"))

fit_firs_rma <- suppressWarnings(metafor::rma.glmm(x1i = ai, t1i = n1i, x2i = ci, t2i = n2i,
                                                   data = data, measure = "IRR", drop00 = FALSE,
                                                   model = "UM.FS"))

test_that("compare rareGLMM to metafor",{

  # rirs:
  expect_equal(fit_rirs$beta[2], fit_rirs_rma$beta, ignore_attr = TRUE,
               tolerance = 1e-5)
  expect_equal(fit_rirs$se[2], fit_rirs_rma$se, ignore_attr = TRUE,
               tolerance = 1e-5)
  expect_equal(fit_rirs$tau2, fit_rirs_rma$tau2, ignore_attr = TRUE,
               tolerance = 1e-5)
  expect_equal(fit_rirs$LRT.Chisq, fit_rirs_rma$QE.LRT, ignore_attr = TRUE,
               tolerance = 1e-5)
  expect_equal(fit_rirs$LRT.df, fit_rirs_rma$QE.df, ignore_attr = TRUE)

  # firs:
  expect_equal(fit_firs$beta["group"], fit_firs_rma$beta, ignore_attr = TRUE,
               tolerance = 1e-4)
  expect_equal(fit_firs$se["group"], fit_firs_rma$se, ignore_attr = TRUE,
               tolerance = 1e-5)
  expect_equal(fit_firs$tau2, fit_firs_rma$tau2, ignore_attr = TRUE,
               tolerance = 1e-5)
  expect_equal(fit_firs$LRT.Chisq, fit_firs_rma$QE.LRT, ignore_attr = TRUE,
               tolerance = 1e-5)
  expect_equal(fit_rirs$LRT.df, fit_rirs_rma$QE.df, ignore_attr = TRUE)




})

##### CHECKING EQUIVALENT WAYS OF DATA INPUT -----------------------------------

test_that("does not matter if data is put in as data frame or rareData-object",{

  k <- which(names(rareGLMM(x=x, measure = "logOR", intercept = "random", slope = "random")) == "call")

  expect_identical(
    rareGLMM(x=x, measure = "logOR", intercept = "random", slope = "random")[-k],
    rareGLMM(ai=ai, bi=bi, ci=ci, di=di, n1i=n1i, n2i=n2i, data=data, measure = "logOR", intercept = "random", slope = "random")[-k],
  )

  expect_identical(
    rareGLMM(x=x, measure = "logOR", intercept = "random", slope = "random")[-k],
    rareGLMM(ai=ai, ci=ci, n1i=n1i, n2i=n2i, data=data, measure = "logOR", intercept = "random", slope = "random")[-k]
  )

  expect_identical(
    rareGLMM(x=x, measure = "logOR", intercept = "random", slope = "random")[-k],
    rareGLMM(ai=ai, ci=ci, di=di, n1i=n1i, data=data, measure = "logOR", intercept = "random", slope = "random")[-k]
  )

  expect_identical(
    rareGLMM(x=x, measure = "logOR", intercept = "random", slope = "random")[-k],
    rareGLMM(ai=ai, bi=bi, ci=ci, di=di, data=data, measure = "logOR", intercept = "random", slope = "random")[-k]
  )
})

