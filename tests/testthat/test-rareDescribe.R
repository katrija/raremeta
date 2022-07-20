test_that("rareDescribe returns the correct number of studies", {
  data <- data.frame(
    ai = c(0,3,2,0),
    bi = c(20,18,15,19),
    ci = c(1,4,0,0),
    di = c(19,17,16,20)
  )

  desc <- rareDescribe(ai = ai, bi = bi, ci = ci, di = di, data = data)

  expect_equal(desc$k, 4)
})
