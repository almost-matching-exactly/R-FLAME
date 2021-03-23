# test warnings

test_that("complains about continuous variables", {
  p <- 5
  weights <- runif(p)
  data <- gen_data(p = p)
  holdout <- gen_data(p = p)
  data$age <- sample(15:100, nrow(data), TRUE)
  holdout$age <- sample(15:100, nrow(holdout), TRUE)
  expect_warning(flout <- FLAME(data, holdout, weights = weights))
})
