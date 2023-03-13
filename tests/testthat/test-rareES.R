#Test file for raremeta-function 'rareES'

#since the input is already tested by 'rareCC', we do only write tests
#for computation of estimated effect sizes and estimated variances

#Compare rareIV computation to metafor's rma

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

test_that("compare to metafor::rma.uni() calculations",{
  #odds ratio
  x.rareES <- rareES(x, cc= "constant", ccval = 0.5, ccto = "only0", drop00 = FALSE,
                     measure = "logOR")
  x.escalc <- metafor::escalc(measure = "OR", ai=ai, bi=bi, ci=ci, di=di, data=data,
                              add = 0.5, to = "only0", drop00=FALSE)

  expect_equal(round(x.rareES$yi, digits = 4),
               as.numeric(round(x.escalc$yi, digits = 4)))

  expect_equal(round(x.rareES$vi, digits = 4),
               as.numeric(round(x.escalc$vi, digits = 4)))

  #relative risk
  x.rareES <- rareES(x, cc= "constant", ccval = 0.5, ccto = "only0", drop00 = FALSE,
                     measure = "logRR")
  x.escalc <- metafor::escalc(measure = "RR", ai=ai, bi=bi, ci=ci, di=di, data=data,
                              add = 0.5, to = "only0", drop00=FALSE)

  expect_equal(round(x.rareES$yi, digits = 4),
               as.numeric(round(x.escalc$yi, digits = 4)))

  expect_equal(round(x.rareES$vi, digits = 4),
               as.numeric(round(x.escalc$vi, digits = 4)))

  #risk difference
  x.rareES <- rareES(x, cc= "constant", ccval = 0.5, ccto = "only0", drop00 = FALSE,
                     measure = "RD")
  x.escalc <- metafor::escalc(measure = "RD", ai=ai, bi=bi, ci=ci, di=di, data=data,
                              add = 0.5, to = "only0", drop00=FALSE)

  expect_equal(round(x.rareES$yi, digits = 4),
               as.numeric(round(x.escalc$yi, digits = 4)))

  expect_equal(round(x.rareES$vi, digits = 4),
               as.numeric(round(x.escalc$vi, digits = 4)))
})
