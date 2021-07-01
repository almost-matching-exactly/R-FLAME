p <- 4
n <- 250
weights <- runif(p)
data <- gen_data(n, p)
holdout <- gen_data(n, p)

test_that("binary numerics must be 0-1", {
  tmp_data <- data
  tmp_holdout <- holdout
  tmp_data$outcome <- sample(1:2, n, replace = T)
  tmp_holdout$outcome <- sample(0:1, n, replace = T)
  expect_error(FLAME(tmp_data, tmp_holdout, weights = weights))
})

test_that("data can be numeric", {
  tmp_data <- data
  tmp_holdout <- holdout
  tmp_data[, c(1:p + 2)] <- lapply(tmp_data[, c(1:p + 2)], as.numeric)
  tmp_holdout[, c(1:p + 2)] <- lapply(tmp_holdout[, c(1:p + 2)], as.numeric)
  flout <- FLAME(tmp_data, tmp_holdout, weights = weights)
  expect_true(TRUE)
})

test_that("covariates can be of mixed types", {
  tmp_data <- data
  tmp_holdout <- holdout
  tmp_data[1] <- lapply(tmp_data[1], as.numeric)
  tmp_data[2] <- lapply(tmp_data[2], as.character)
  tmp_data[3] <- lapply(tmp_data[3], as.factor)
  flout <- FLAME(tmp_data, tmp_holdout, weights = weights)
  expect_true(TRUE)
})

test_that("outlier units don't break anything", {
  tmp_data <- data
  tmp_holdout <- holdout
  max_val <- max(apply(data[, c(1:4)], 2, max))
  for (unit in sample(1:n, n %/% 20, replace=FALSE)) {
    tmp_data[unit, sample(1:p, 1)] <- max_val
    max_val <- max_val + 1
  }
  flout <- FLAME(tmp_data, tmp_holdout, weights = weights)
  expect_true(TRUE)
})
