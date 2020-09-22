data <- gen_data(500)
holdout <- gen_data(500)

test_that("user specified PE methods work", {

  my_fit <- glmnet::cv.glmnet
  my_fit_params <- list(nfolds = 5)
  my_predict <- predict
  my_predict_params <- list(s = "lambda.min")

  FLAME_out <-
    FLAME(data = data, holdout = holdout,
          user_PE_fit = my_fit, user_PE_fit_params = my_fit_params,
          user_PE_predict = my_predict)

  expect_true(TRUE)
})

