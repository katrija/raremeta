#tests for the function 'rarePeto()'
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

#runs with legitimate inputs
test_that("rarePeto() runs with valid inputs",{
  expect_error(rarePeto(x),
               NA)

})



#Comparing to metafor::rma.peto

test_that("comparing to metafor::rma.peto()",{
  rareM     <- rarePeto(x)
  metaforM  <- suppressWarnings(metafor::rma.peto(ai=ai,bi=bi,ci=ci,di=di,data=data))

  expect_equal(as.numeric(rareM["beta"]),
               as.numeric(metaforM["beta"]))

  expect_equal(as.numeric(rareM["se"]),
               as.numeric(metaforM["se"]))

  expect_equal(as.numeric(rareM["zval"]),
               as.numeric(metaforM["zval"]))

  expect_equal(as.numeric(rareM["pval"]),
               as.numeric(metaforM["pval"]))

  expect_equal(as.numeric(rareM["ci.lb"]),
               as.numeric(metaforM["ci.lb"]))

  expect_equal(as.numeric(rareM["ci.ub"]),
               as.numeric(metaforM["ci.ub"]))




})

# checking equivalent ways of data input

test_that("does not matter if data is put in as data frame or rareData-object",{
  expect_equal(
    rarePeto(x=x),
    rarePeto(ai=ai, bi=bi, ci=ci, di=di, n1i=n1i, n2i=n2i, data=data)
  )

  expect_equal(
    rarePeto(x=x),
    rarePeto(ai=ai, ci=ci, n1i=n1i, n2i=n2i, data=data)
  )

  expect_equal(
    rarePeto(x=x),
    rarePeto(ai=ai, ci=ci, di=di, n1i=n1i, data=data)
  )

  expect_equal(
    rarePeto(x=x),
    rarePeto(ai=ai, bi=bi, ci=ci, di=di, data=data)
  )


})
