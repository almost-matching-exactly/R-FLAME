# Tests ensuring functionality when only one covariate is supplied
n <- 250
p <- 1
df <- data.frame(X = sample(0:3, n, TRUE),
                 treated = sample(c(T, F), n, TRUE))
df$outcome <- df$X + 5 * df$treated + rnorm(n)

holdout <- data.frame(X = sample(0:3, n, TRUE),
                 treated = sample(c(T, F), n, TRUE))
holdout$outcome <- holdout$X + 5 * holdout$treated + rnorm(n)

flout <- FLAME(df, holdout)

test_that("runs without error", {
  expect_true(TRUE)
})

test_that("column names are correct", {
  expect_equal(colnames(flout$data), c('X', 'outcome', 'treated', 'matched', 'weight'))
})

df$X <- sample(c('Red', 'White', 'Green'), n, TRUE)
holdout$X <- sample(c('Red', 'White', 'Green'), n, TRUE)

test_that("runs with non-numeric covariate", {
  flout <- FLAME(df, holdout)
  expect_true(TRUE)
})

test_that("works without outcome", {
  df$outcome <- NULL
  flout <- FLAME(df, holdout)
  expect_true(TRUE)
})

df$outcome <- sample(0:3, n, TRUE) + 5 * df$treated + rnorm(n)

test_that("column order doesn't matter", {
  df <- df[, c(2, 3, 1)]
  holdout <- holdout[, c(2, 3, 1)]
  flout <- FLAME(df, holdout)
  expect_true(TRUE)
})

df$X[sample(1:n, 10)] <- NA
holdout$X[sample(1:n, 10)] <- NA

test_that("runs with 'drop' missingness", {
  flout <- FLAME(df, holdout, missing_data = 'drop', missing_holdout = 'drop')
  expect_true(TRUE)
})

test_that("runs with 'impute' missingness", {
  flout <- FLAME(df, holdout, missing_data = 'impute', missing_holdout = 'impute')
  expect_true(TRUE)
})

test_that("runs with 'ignore' missingness", {
  flout <- FLAME(df, holdout, missing_data = 'ignore', missing_holdout = 'drop')
  expect_true(TRUE)
})

