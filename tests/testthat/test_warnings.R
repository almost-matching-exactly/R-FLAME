# test warnings

test_that("complains about continuous variables", {
  data <- gen_data()
  holdout <- gen_data()
  data$age <- sample(15:100, nrow(data), TRUE)
  holdout$age <- sample(15:100, nrow(holdout), TRUE)
  expect_warning(flout <- FLAME(data, holdout))
})
