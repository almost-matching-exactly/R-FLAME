n <- 250
p <- 5
data <- gen_data(n, p)
holdout <- gen_data(n, p)

colnames(data) <- c('X1', 'X2', 'X3', 'X4', 'X5', 'outcome', 'treated')
colnames(holdout) <- c('X1', 'X2', 'X3', 'X4', 'X5', 'outcome', 'treated')

data$outcome <- sample(c('Sick', 'Super Sick', 'A Little Sick', 'Not Sick'),
                       n, replace = TRUE)

data$outcome <- factor(data$outcome)

holdout$outcome <- sample(c('Sick', 'Super Sick', 'A Little Sick', 'Not Sick'),
                       n, replace = TRUE)

holdout$outcome <- factor(holdout$outcome)

test_that("works with multinomial outcome", {
  flout <- FLAME(data, holdout)
  expect_true(TRUE)

  expect_null(flout$CATE)
})

data$outcome <- sample(c('Sick', 'Not Sick'),
                       n, replace = TRUE)

data$outcome <- factor(data$outcome)

holdout$outcome <- sample(c('Sick', 'Not Sick'),
                          n, replace = TRUE)

holdout$outcome <- factor(holdout$outcome)

test_that("works with binary outcome", {
  flout <- FLAME(data, holdout)
  expect_true(TRUE)

  expect_null(flout$CATE)
})
