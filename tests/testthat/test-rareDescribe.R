data <- data.frame(
  ai = c(0,3,2,0),
  bi = c(20,18,15,19),
  ci = c(1,4,0,0),
  di = c(19,17,16,20),
  n1i = c(20, 21, 17, 19),
  n2i = c(20, 21, 16, 20)
)

dataCC <- data.frame(
  ai = c(0.5,3,2,0.5),
  bi = c(20,18,15,19),
  ci = c(1,4,0.5,0.5),
  di = c(19,17,16,20),
  n1i = c(20.5, 21, 17, 19.5),
  n2i = c(20, 21, 16.5, 20.5)
)

dataDesc <- rareDescribe(ai = ai, bi = bi, ci = ci, di = di, data = data)

test_that("rareDescribe returns errors and warning messages", {
  # error if ai is missing
  expect_error(rareDescribe(bi = bi, ci = ci, di = di, data = data), "ai must be specified.")

  # error if bi and n1i are missing
  expect_error(rareDescribe(ai = ai, ci = ci, di = di, data = data), "bi or n1i must be specified.")

  # error if ci is missing
  expect_error(rareDescribe(ai = ai, bi = bi, di = di, data = data), "ci must be specified.")

  # error if di and n2i are missing
  expect_error(rareDescribe(ai = ai, bi = bi, ci = ci, data = data), "di or n2i must be specified.")

  # error if data is no data.frame
  list = as.list(data)
  expect_error(rareDescribe(ai = ai, bi = bi, ci = ci, di = di, n1i = n1i, n2i = n2i, data = list), "Argument 'data' must be a data frame.")

  # error if negative arguments occur
  data$negative.ai = c(-1,3,2,0)
  expect_error(rareDescribe(ai = negative.ai, bi = bi, ci = ci, di = di, n1i = n1i, n2i = n2i, data = data), "Arguments must not be negative.")

  # error if ai + bi do not add up to n1i
  data$add.ai <- c(1,3,2,0)
  expect_error(rareDescribe(ai = add.ai, bi = bi, ci = ci, di = di, n1i = n1i, n2i = n2i, data = data), "ai and bi must add up to n1i.")

  # error if ci + di do not add up to n2i
  data$add.ci <- c(1,4,1,0)
  expect_error(rareDescribe(ai = ai, bi = bi, ci = add.ci, di = di, n1i = n1i, n2i = n2i, data = data), "ci and di must add up to n2i.")

  # warning if there are missing values in group 1
  data$NA.ai <- c(0,3,2,NA)
  expect_warning(rareDescribe(ai = NA.ai, bi = bi, ci = ci, di = di, n1i = n1i, n2i = n2i, data = data), "There are missing values in group 1 in studies: ")

  # warning if there are missing values in group 2
  data$NA.ci <- c(1,4,0,NA)
  expect_warning(rareDescribe(ai = ai, bi = bi, ci = NA.ci, di = di, n1i = n1i, n2i = n2i, data = data), "There are missing values in group 2 in studies: ")

  # warning if arguments are not integer
  #expect_warning(rareDescribe(ai = ai, bi = bi, ci = ci, di = di, data = dataCC), "Some values in ai, bi, ci or di are not integers. Please check whether this is intended.")
})

test_that("summaries are calculated correctly", {
  ni = data$n1i + data$n2i
  n <- c(
    mean(ni, na.rm = TRUE), stats::median(ni, na.rm = TRUE),
    stats::quantile(ni, c(0.25, 0.75), na.rm = TRUE),
    min(ni, na.rm = TRUE), max(ni, na.rm = TRUE)
  )
  names(n) <- c("mean", "median", "q25", "q75", "min", "max")
  expect_equal(dataDesc$n, n)
})

test_that("rareDescribe correctly calculates bi and di", {
  # does rareDescribe correctly calculate bi?
  expect_equal(dataDesc$bi, data$bi)

  # does rareDescribe correctly calculate di?
  expect_equal(dataDesc$di, data$di)
})

test_that("rareDescribe correctly sums up numbers of studies", {
  # counting all studies
  expect_equal(dataDesc$k, 4)

  # counting double zero studies
  expect_equal(dataDesc$kdz, 1)

  # counting single zero studies
  expect_equal(dataDesc$ksz, 2)

  # counting single zero studies in group 1
  expect_equal(dataDesc$k1sz, 1)

  # counting single zero studies in group 2
  expect_equal(dataDesc$k2sz, 1)
})
