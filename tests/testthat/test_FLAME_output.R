test_that('all correct drops', {
  data <- gen_data(500, 8)
  holdout <- gen_data(500, 8)
  flout <- FLAME(data = data, holdout = holdout)
  expect_true(check_all_matches(flout, data))
})

test_that('dropped_covs and matching_covs in sync', {
  flout <- FLAME(data = gen_data(200, 10), holdout = gen_data(200, 10))
  covs <- c('X1', 'X2', 'X3', 'X4', 'X5', 'X6', 'X7', 'X8', 'X9', 'X10')
  dropped_covs <- flout$dropped
  matching_covs <- flout$matching_covs

  if (!isTRUE(all.equal(matching_covs[[1]], covs))) {
    matching_covs <- c(list(covs), matching_covs)
  }

  expect_equal(length(setdiff(covs, dropped_covs)), 1)
  expect_equal(length(dropped_covs) + 1, length(matching_covs))
  for (i in 1:length(dropped_covs)) {
    expect_equal(dropped_covs[i], setdiff(matching_covs[[i]], matching_covs[[i + 1]]))
  }
})

test_that('covariates dropped in right order', {
  flout <-
    FLAME(data = gen_data(200, 10), holdout = gen_data(200, 10), C = 0)
  n_dropped <- length(flout$dropped)
  expect_equal(flout$dropped[n_dropped - 1], 'X3')
  expect_equal(flout$dropped[n_dropped], 'X2')
})
