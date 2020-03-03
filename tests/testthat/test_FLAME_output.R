# test_that('all correct drops', {
#   data <- gen_data(500, 8)
#   holdout <- gen_data(500, 8)
#   flout <- FLAME(data = data, holdout = holdout)
#   expect_true(check_all_matches(flout, data))
# })

test_that('dropped_covs and matching_covs in sync', {
  flout <- FLAME(data = gen_data(200, 10), holdout = gen_data(200, 10))
  covs <- c('X1', 'X2', 'X3', 'X4', 'X5', 'X6', 'X7', 'X8', 'X9', 'X10')
  dropped_covs <- flout$dropped
  matching_covs <- flout$matching_covs

  if (!isTRUE(all.equal(matching_covs[[1]], covs))) {
    matching_covs <- c(list(covs), matching_covs)
  }

  expect_equal(length(dropped_covs) + 1, length(matching_covs))
  for (i in 1:length(dropped_covs)) {
    expect_equal(dropped_covs[i], setdiff(matching_covs[[i]], matching_covs[[i + 1]]))
  }
})

is_nested <- function(lst) {
  if (length(lst) == 1) {
    return(TRUE)
  }
  for (i in 1:(length(lst) - 1)) {
    if (!(all(lst[[i]] %in% lst[[i + 1]]))) {
      return(FALSE)
    }
  }
  return(TRUE)
}

test_that('matches propagate down if replace = TRUE', {
  n <- 200
  p <- 8
  flout <- FLAME(data = gen_data(n, p), holdout = gen_data(n, p), replace = TRUE)
  min_matched_on <- p - length(flout$dropped)
  for (i in 1:n) {
    if (!flout$data$matched[i]) {
      next
    }
    MGs <- which(sapply(flout$MGs, function(x) i %in% x))
    expect_true(is_nested(flout$MGs[MGs]))
    n_matched <- vapply(MGs, function(x) length(flout$matched_on[[x]]), numeric(1))
    if (length(n_matched) == 1) {
      expect_equal(n_matched, min_matched_on)
    }
    else {
      expect_equal(n_matched, n_matched[1]:min_matched_on)
    }
  }
})

# test_that('covariates dropped in right order', {
#   flout <-
#     FLAME(data = gen_data(200, 10), holdout = gen_data(200, 10), C = 0)
#   n_dropped <- length(flout$dropped)
#   expect_equal(flout$dropped[n_dropped - 1], 'X3')
#   expect_equal(flout$dropped[n_dropped], 'X2')
# })
