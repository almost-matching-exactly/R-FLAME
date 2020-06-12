data <- gen_data(500)
holdout <- gen_data(500)

test_that("can use BART to predict", {

  my_fit <- dbarts::bart
  my_fit_params <- list(ntree = 100, verbose = FALSE, keeptrees = TRUE)
  my_predict <- function(bart_fit, new_data) {
    return(colMeans(predict(bart_fit, new_data)))
  }
  FLAME_out <-
    FLAME(data = data, holdout = holdout,
          user_PE_fit = my_fit, user_PE_fit_params = my_fit_params,
          user_PE_predict = my_predict)

  expect_true(TRUE)
})

