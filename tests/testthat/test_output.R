if (F) {
  skip('Ignoring output tests')
}

# test_that('all correct drops', {
#   data <- gen_data(500, 8)
#   holdout <- gen_data(500, 8)
#   flout <- FLAME(data = data, holdout = holdout)
#   expect_true(check_all_matches(flout, data))
# })

p <- 5
weights <- abs(rnorm(p))

test_that('FLAME MGs contain self', {
  dat <- gen_data(p = p)
  holdout <- gen_data(p = p)
  flout <- FLAME(dat, holdout, weights = weights)
  for (i in seq_along(flout$MGs)) {
    if (is.null(flout$MGs[[i]])) {
      next
    }
    expect_true(i %in% flout$MGs[[i]])
  }

  flout <- FLAME(dat, holdout, replace = TRUE, weights = weights)
  for (i in seq_along(flout$MGs)) {
    if (is.null(flout$MGs[[i]])) {
      next
    }
    expect_true(i %in% flout$MGs[[i]])
  }
})

test_that('DAME MGs contain self', {
  dat <- gen_data(p = p)
  holdout <- gen_data(p = p)
  dout <- DAME(dat, holdout, weights = weights)
  for (i in seq_along(dout$MGs)) {
    if (is.null(dout$MGs[[i]])) {
      next
    }
    expect_true(i %in% dout$MGs[[i]])
  }

  dout <- DAME(dat, holdout, replace = TRUE, weights = weights)
  for (i in seq_along(dout$MGs)) {
    if (is.null(dout$MGs[[i]])) {
      next
    }
    expect_true(i %in% dout$MGs[[i]])
  }
})
