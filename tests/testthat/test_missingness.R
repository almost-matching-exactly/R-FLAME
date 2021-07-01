strip_missing <- function(vals) {
  # Gets rid of all " (m)" substrings at the end of values,
  # then removes all "*"
  tmp <- strsplit(as.character(vals), ' ')
  tmp <- lapply(tmp, function(x) {
    if (x[length(x)] == '(m)') {
      x <- x[length(x) - 1]
    }
    return(x)
  })
  return(setdiff(sapply(tmp, paste, collapse = ' ', sep = ''), '*'))
}

test_that("missing_data `drop` works", {
  p <- 5
  n <- 2500
  weights <- runif(p)
  data <- gen_missing_data(n, p)
  holdout <- gen_data(n, p)

  n_miss <- 250

  data[arrayInd(sample(n * p, n_miss), c(n, p))] <- NA

  # Drop missingness before
  ATE_predrop <- ATE(FLAME(data[apply(data, 1, function(x) !any(is.na(x))), ],
                           holdout, weights = weights))

  # Drop missingness within algo
  ATE_algodrop <-
    ATE(FLAME(data, holdout, missing_data = 'drop', weights = weights))

  expect_equal(ATE_predrop, ATE_algodrop)

  ATE_predrop <- ATE(DAME(data[apply(data, 1, function(x) !any(is.na(x))), ],
                          holdout, weights = weights))

  # Drop missingness within algo
  ATE_algodrop <-
    ATE(DAME(data, holdout, missing_data = 'drop', weights = weights))

  expect_equal(ATE_predrop, ATE_algodrop)
})

test_that("replaced values don't show", {
  p <- 3
  weights <- runif(p)
  data <- gen_missing_data(n = 250, p = p)
  holdout <- gen_data(n = 250, p = p)
  flout <- FLAME(data = data, holdout = holdout,
                 missing_data = 'impute', weights = weights)
  no_extra_vals <-
    vapply(1:p, function(cov) {
      length(setdiff(unique(flout$data[[cov]]),
                     unique(data[[cov]]))) == 0
    }, logical(1))
  expect_true(all(no_extra_vals))
})

p <- 4
weights <- runif(p)
n <- 250
data <- gen_data(n = n, p = p)
holdout <- gen_data(n = n, p = p)
flout <- FLAME(data = data, holdout = holdout, weights = weights)
replace_inds_data <- c(sample(1:n, 1), sample(1:p, 1))
replace_inds_holdout <- c(sample(1:n, 1), sample(1:p, 1))
data[replace_inds_data[1], replace_inds_data[2]] <- NA
holdout[replace_inds_holdout[1], replace_inds_holdout[2]] <- NA

# Former matched group of now missing unit
MG_of_missing <- MG(replace_inds_data[1], flout, index_only = TRUE)

# Did the unit originally match on the value they're now missing
matched_on_missing <-
  flout$data[replace_inds_data[1], replace_inds_data[2]] != '*'

test_that("dropping missing data works", {
  flout1 <-
    FLAME(data = data, holdout = holdout,
          missing_data = 'drop', missing_holdout = 'drop', weights = weights)

  # Avoid case in which the unit made missing was the only match for another
  # unit equivalent, not identical, to ignore discrepancies in factor levels due
  # to ' (m)'
  if (length(MG_of_missing) > 2) {
    expect_equivalent(flout$data[-replace_inds_data[1], ],
                     flout1$data[-replace_inds_data[1], ])
  }
  else {
    expect_true(TRUE)
  }
})

test_that("not matching on missing data works", {
  n_imps <- 2
  flout2 <- FLAME(data = data, holdout = holdout,
                  missing_data = 'ignore', missing_holdout = 'drop',
                  missing_data_imputations = n_imps, weights = weights)

  for (i in 1:n_imps) {
    # the missingness may have made me eligible
    # for someone else that didn't otherwise get matched
    if (length(MG_of_missing) > 2 &
        length(MG(replace_inds_data[1], flout2, index_only = TRUE)) > 2) {
      expect_identical(flout$data[-replace_inds_data[1], ],
                       flout2[[i]]$data[-replace_inds_data[1], ])
    }
    else {
      expect_true(TRUE)
    }
  }
})

# Check if we changed the output format here
test_that("missing option 3 works", {
  flout3 <-
    FLAME(data = data, holdout = holdout,
          missing_data = 'impute', missing_holdout = 'drop', weights = weights)

  flout$data[[replace_inds_data[[2]]]] <-
    factor(flout$data[[replace_inds_data[[2]]]],
           levels = levels(flout3$data[[replace_inds_data[[2]]]]))

  if (matched_on_missing) {
    if (length(MG_of_missing) > 2 &
        length(MG(replace_inds_data[1], flout3, index_only = TRUE)) > 2) {
      expect_identical(flout$data[-replace_inds_data[1], ],
                       flout3$data[-replace_inds_data[1], ])
    }
    else {
      expect_true(TRUE)
    }
  }
  else {
    flout3$data[which(flout3$data == '* (m)', arr.ind = T)] <- '*'
    expect_identical(flout$data, flout3$data)
  }
})

test_that("imputation works with no outcome", {
  data$outcome <- NULL
  flout <- FLAME(data, holdout, missing_data = 'ignore', missing_holdout = 'drop', weights = weights)
  expect_true(TRUE)
})
