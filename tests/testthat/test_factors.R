p <- 5
data <- gen_data(n = 250, p = p)
holdout <- gen_data(n = 250, p = p)
flout_not_0_ind <- FLAME(data, holdout, C = 1e5)

data_0_ind <- data
holdout_0_ind <- holdout
mapping <- as.character(c(1, 2, 3, 4))
names(mapping) <- as.character(c(0, 1, 2, 3))
for (j in 1:p) {
  data_0_ind[[j]] <- do.call(forcats::fct_recode, c(list(data_0_ind[[j]]), mapping))
  holdout_0_ind[[j]] <- do.call(forcats::fct_recode, c(list(holdout_0_ind[[j]]), mapping))
}

flout_0_ind <- FLAME(data_0_ind, holdout_0_ind, C = 1e5)

test_that("non 0 indexed factors work", {
  expect_identical(flout_not_0_ind$MGs, flout_0_ind$MGs)
})

test_that("non consecutive-level factors work", {
  mapping <- as.character(c(1, 2, 3, 4))
  names(mapping) <- as.character(c(1, 3, 5, 8))
  data_non_consec <- data
  holdout_non_consec <- holdout
  for (j in 1:p) {
    data_non_consec[[j]] <-
      do.call(forcats::fct_recode, c(list(data_non_consec[[j]]), mapping))
    holdout_non_consec[[j]] <-
      do.call(forcats::fct_recode, c(list(holdout_non_consec[[j]]), mapping))
  }
  flout_non_consec <- FLAME(data_non_consec, holdout_non_consec, C = 1e5)
  expect_identical(flout_non_consec$MGs, flout_0_ind$MGs)
})


test_that("non numeric factors work", {
  mapping <- as.character(c(1, 2, 3, 4))
  names(mapping) <- c('white', 'hispanic', 'black', 'asian')
  data_non_num <- data
  holdout_non_num <- holdout
  for (j in 1:p) {
    data_non_num[[j]] <-
      do.call(forcats::fct_recode, c(list(data_non_num[[j]]), mapping))
    holdout_non_num[[j]] <-
      do.call(forcats::fct_recode, c(list(holdout_non_num[[j]]), mapping))
  }
  flout_non_num <- FLAME(data_non_num, holdout_non_num, C = 1e5)
  expect_identical(flout_non_num$MGs, flout_0_ind$MGs)
})

test_that("factors play nice with missing data", {

})

test_that("factors play nice with continuous data", {

})
