## Tests relating to column order
p <- 4
data <- gen_data(n = 250, p = p)
holdout <- gen_data(n = 250, p = p)
flout <- FLAME(data = data, holdout = holdout)

test_that("column order doesn't matter", {
  scrambling <- order(sample(1:(p + 2)))
  scrambled_data <- data[, scrambling]
  scrambled_holdout <- holdout[, scrambling]
  scrambled_flout <- FLAME(data = scrambled_data, holdout = scrambled_holdout)
  if (!identical(flout$MGs, scrambled_flout))
  expect_identical(flout$MGs, scrambled_flout$MGs)
})

test_that("breaks with extra covariates", {
  data$extra <- factor(sample(0:5, nrow(data), T))
  expect_error(flout(data, holdout))

  holdout$extra <- data$extra
  data$extra <- NULL
  expect_error(flout(data, holdout))

  holdout$extra <- NULL
})

test_that("excluding outcome doesn't change matches", {
  data$outcome <- NULL

  flout_wo_outcome <- FLAME(data = data, holdout = holdout)

  expect_identical(flout$MGs, flout_wo_outcome$MGs)
})

test_that("column order doesn't matter when no outcome", {
  tmp <- holdout$outcome
  holdout$outcome <- NULL
  holdout <- cbind(outcome = tmp, holdout)
  flout <- FLAME(data, holdout)
  expect_true(TRUE)
})

p <- 4
data <- gen_data(n = 250, p = p)
holdout <- gen_data(n = 250, p = p)

## Tests relating to column naming
test_that("outcome/treatment name doesn't matter", {
  renamed_data <- data
  renamed_holdout <- holdout

  colnames(renamed_data)    <- c('X1', 'X2', 'X3', 'X4', 'myout', 'myt')
  colnames(renamed_holdout) <- c('X1', 'X2', 'X3', 'X4', 'myout', 'myt')

  flout <- FLAME(data = data, holdout = holdout)
  renamed_flout <- FLAME(data = renamed_data, holdout = renamed_holdout,
                         treated_column_name = 'myt',
                         outcome_column_name = 'myout')

  expect_identical(flout$MGs, renamed_flout$MGs)
})

names_w_spaces <- c('Var 1', 'Var 2', 'Var 3', 'Var 4', 'My Outcome', 'The Treatment')
colnames(data) <- names_w_spaces
colnames(holdout) <- names_w_spaces

test_that("column names with spaces are ok", {
  flout <- FLAME(data, holdout,
                 treated_column_name = 'The Treatment', outcome_column_name = 'My Outcome')
  expect_true(TRUE)
})
