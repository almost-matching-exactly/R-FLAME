test_that("column order doesn't matter", {
  p <- 4
  data <- gen_data(n = 250, p = p)
  holdout <- gen_data(n = 250, p = p)
  flout <- FLAME(data = data, holdout = holdout)

  scrambling <- order(sample(1:(p + 2)))
  scrambled_data <- data[, scrambling]
  scrambled_holdout <- holdout[, scrambling]
  scrambled_flout <- FLAME(data = scrambled_data, holdout = scrambled_holdout)
  if (!identical(flout$MGs, scrambled_flout))
  expect_identical(flout$MGs, scrambled_flout$MGs)
})


test_that("outcome/treatment name doesn't matter", {
  p <- 4

  data <- gen_data(n = 250, p = p)
  holdout <- gen_data(n = 250, p = p)

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

