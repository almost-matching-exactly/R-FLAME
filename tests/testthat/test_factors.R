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
  return(unique(setdiff(sapply(tmp, paste, collapse = ' ', sep = ''), '*')))
}

n <- 250
p <- 5
weights <- runif(p)
data <- gen_data(n = n, p = p)
holdout <- gen_data(n = n, p = p)
flout_not_0_ind <- FLAME(data, holdout, weights = weights)

## Careful because these tests depend on the covariates generated
## in gen_data.

mapping <- as.character(c(1, 2, 3, 4))
names(mapping) <- as.character(c(0, 1, 2, 3))
mapping <- rep(list(mapping), p)

data_0_ind <- factor_remap(data, 'treated', 'outcome', mapping = mapping)$df
holdout_0_ind <-
  factor_remap(holdout, 'treated', 'outcome', mapping = mapping)$df

flout_0_ind <- FLAME(data_0_ind, holdout_0_ind, weights = weights)

test_that("non 0 indexed factors work", {
  expect_identical(flout_not_0_ind$MGs, flout_0_ind$MGs)
})

test_that("non consecutive-level factors work", {
  mapping <- as.character(c(1, 2, 3, 4))
  names(mapping) <- as.character(c(1, 3, 5, 8))
  mapping <- rep(list(mapping), p)

  data_non_consec <-
    factor_remap(data, 'treated', 'outcome', mapping = mapping)$df
  holdout_non_consec <-
    factor_remap(holdout, 'treated', 'outcome', mapping = mapping)$df

  flout_non_consec <-
    FLAME(data_non_consec, holdout_non_consec, weights = weights)
  expect_identical(flout_non_consec$MGs, flout_0_ind$MGs)
})

test_that("non numeric factors work", {
  mapping <- as.character(c(1, 2, 3, 4))
  names(mapping) <- c('white', 'american indian', 'black', 'asian')
  mapping <- rep(list(mapping), p)

  data_non_num <-
    factor_remap(data, 'treated', 'outcome',  mapping = mapping)$df
  holdout_non_num <-
    factor_remap(holdout, 'treated', 'outcome',  mapping = mapping)$df

  flout_non_num <- FLAME(data_non_num, holdout_non_num, weights = weights)
  expect_identical(flout_non_num$MGs, flout_0_ind$MGs)
})

test_that("missing data 'impute' doesn't leave new levels", {
  mapping <- as.character(c(1, 2, 3, 4))
  names(mapping) <- c('white', 'american indian', 'black', 'asian')
  mapping <- rep(list(mapping), p)

  data_non_num <-
    factor_remap(data, 'treated', 'outcome', mapping = mapping)$df
  holdout_non_num <-
    factor_remap(holdout, 'treated', 'outcome', mapping = mapping)$df

  levels_in <- lapply(data_non_num[, 1:p], function(x) levels(factor(x)))
  for (i in 1:p) {
    data_non_num[[i]][sample(1:n, 10)] <- NA
  }
  flout <- FLAME(data_non_num, holdout_non_num,
                 missing_data = 'impute', weights = weights)
  levels_out <- lapply(flout$data[, 1:p], levels)
  no_extra_vals <-
    vapply(1:p, function(cov) {
      length(setdiff(strip_missing(flout$data[[cov]]),
                     unique(data_non_num[[cov]]))) == 0
    }, logical(1))
  expect_true(all(no_extra_vals))
})

test_that("we don't reorder factor levels", {
  data[, 1:p] <- lapply(data[, 1:p], as.factor)
  holdout[, 1:p] <- lapply(holdout[, 1:p], as.factor)

  data_reordered <- data
  holdout_reordered <- holdout

  data_reordered$X1 <- relevel(data_reordered$X1, levels(data_reordered$X1)[2])
  holdout_reordered$X1 <-
    relevel(holdout_reordered$X1, levels(holdout_reordered$X1)[2])

  dout_reordered <- DAME(data_reordered, holdout_reordered, weights = weights)
  if (nlevels(data$X1) == nlevels(dout_reordered$data$X1)) {
    expect_true(all(levels(data_reordered$X1) ==
                      levels(dout_reordered$data$X1)))
  }
})
