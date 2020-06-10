if (F) {
  skip('Ignoring missingness tests')
}
test_that("replaced values don't show", {
  p <- 3
  data <- gen_missing_data(n = 250, p = p)
  holdout <- gen_data(n = 250, p = p)
  flout <- FLAME(data = data, holdout = holdout, missing_data = 3)
  no_extra_vals <-
    sapply(1:p, function(cov) {
      length(setdiff(setdiff(unique(flout$data[[cov]]), '*'), unique(data[[cov]]))) == 0
    })
  expect_true(all(no_extra_vals))
})

p <- 3
n <- 250
data <- gen_data(n = n, p = p)
holdout <- gen_data(n = n, p = p)
flout <- FLAME(data = data, holdout = holdout)
replace_inds_data <- c(sample(1:n, 1), sample(1:p, 1))
replace_inds_holdout <- c(sample(1:n, 1), sample(1:p, 1))
data[replace_inds_data[1], replace_inds_data[2]] <- NA
holdout[replace_inds_holdout[1], replace_inds_holdout[2]] <- NA

test_that("missing option 1 works", {
  flout1 <- FLAME(data = data, holdout = holdout,
                  missing_data = 1, missing_holdout = 1)

  expect_identical(flout$data[-replace_inds_data[1], ],
                   flout1$data[-replace_inds_data[1], ])
})

test_that("missing option 2 works", {
  n_imps <- 2
  flout2 <- FLAME(data = data, holdout = holdout,
                  missing_data = 2, missing_holdout = 2,
                  missing_data_imputations = n_imps)

  for (i in 1:n_imps) {
    expect_identical(flout$data[-replace_inds_data[1], ],
                     flout2[[i]]$data[-replace_inds_data[1], ])
  }
})

test_that("missing option 3 works", {
  flout3 <- FLAME(data = data, holdout = holdout,
                  missing_data = 3, missing_holdout = 1)

  flout$data[[replace_inds_data[[2]]]] <-
    factor(flout$data[[replace_inds_data[[2]]]],
           levels = levels(flout3$data[[replace_inds_data[[2]]]]))

  expect_identical(flout$data[-replace_inds_data[1], ],
                   flout3$data[-replace_inds_data[1], ])
})
