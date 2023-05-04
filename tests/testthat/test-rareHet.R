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


# checking equivalent ways of data input

test_that("does not matter if data is put in as data frame or rareData-object",{
  expect_equal(
    rareHet(x=x, measure="logOR", cc="constant"),
    rareHet(ai=ai, bi=bi, ci=ci, di=di, n1i=n1i, n2i=n2i, data=data,  measure="logOR", cc="constant")
  )

  expect_equal(
    rareHet(x=x, measure="logOR", cc="constant"),
    rareHet(ai=ai, ci=ci, n1i=n1i, n2i=n2i, data=data, measure="logOR", cc="constant")
  )

  expect_equal(
    rareHet(x=x, measure="logOR", cc="constant"),
    rareHet(ai=ai, ci=ci, di=di, n1i=n1i, data=data, measure="logOR", cc="constant")
  )

  expect_equal(
    rareHet(x=x, measure="logOR", cc="constant"),
    rareHet(ai=ai, bi=bi, ci=ci, di=di, data=data, measure="logOR", cc="constant")
  )


})
