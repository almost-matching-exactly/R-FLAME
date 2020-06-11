n <- 500
data <- gen_data(n)
holdout <- gen_data(n)

test_that("ridge works with continuous outcome", {
  flout <- FLAME(data, holdout, PE_method = 'ridge')
  expect_true(TRUE)
})

test_that("XGBoost works with continuous outcome", {
  flout <- FLAME(data, holdout, PE_method = 'xgb')
  expect_true(TRUE)
})

data$outcome <-
  exp(scale(data$outcome, scale = F)) /
  (1 + exp(scale(data$outcome, scale = F)))

data$outcome <- round(data$outcome)

holdout$outcome <-
  exp(scale(holdout$outcome, scale = F)) /
  (1 + exp(scale(holdout$outcome, scale = F)))

holdout$outcome <- round(holdout$outcome)

test_that("ridge works with binary continuous outcome", {
  flout <- FLAME(data, holdout, PE_method = 'ridge')
  expect_true(TRUE)
})

test_that("XGBoost works with binary continuous outcome", {
  flout <- FLAME(data, holdout, PE_method = 'xgb')
  expect_true(TRUE)
})

data$outcome <- factor(data$outcome)
holdout$outcome <- factor(holdout$outcome)

test_that("ridge works with binary factor outcome", {
  flout <- FLAME(data, holdout, PE_method = 'ridge')
  expect_true(TRUE)
})

test_that("XGBoost works with binary factor outcome", {
  flout <- FLAME(data, holdout, PE_method = 'xgb')
  expect_true(TRUE)
})

data$outcome <- factor(sample(c('Green', 'White', 'Red'), n, TRUE))
holdout$outcome <- factor(sample(c('Green', 'White', 'Red'), n, TRUE))

test_that("ridge works with multiclass outcomes", {
  flout <- FLAME(data, holdout, PE_method = 'ridge')
  expect_true(TRUE)
})

test_that("XGBoost works with multiclass outcomes", {
  flout <- FLAME(data, holdout, PE_method = 'xgb')
  expect_true(TRUE)
})
