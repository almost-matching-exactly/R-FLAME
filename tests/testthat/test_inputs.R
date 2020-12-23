p <- 4
n <- 250
data <- gen_data(n, p)
holdout <- gen_data(n, p)

test_that("binary numerics must be 0-1", {
  tmp_data <- data
  tmp_holdout <- holdout
  tmp_data$outcome <- sample(1:2, n, replace = T)
  tmp_holdout$outcome <- sample(0:1, n, replace = T)
  expect_error(FLAME(tmp_data, tmp_holdout))
})
