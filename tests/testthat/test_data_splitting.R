n <- 250
p <- 5
data <- gen_data(5 * n, p)

test_that("data splitting doesn't affect post matching", {
  flout <- FLAME(data)
  MGs <- MG(1:nrow(flout$data), flout)
  expect_true(TRUE)

  holdout_inds <- sample(1:nrow(data), round(0.1 * nrow(data)), TRUE)
  holdout <- data[holdout_inds, ]
  data <- data[-holdout_inds, ]

  flout <- FLAME(data, holdout)
  MGs <- MG(1:nrow(flout$data), flout)
  expect_true(TRUE)
})

test_that("can't require data splitting without outcome", {
  data$outcome <- NULL
  expect_error(flout <- FLAME(data))
})
