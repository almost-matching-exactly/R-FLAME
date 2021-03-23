# Four options:
# 1. OK: PE_method is ridge/xgb and user_PE* is NULL
# 2. PE_method is a function and user_PE* is NULL
# 3. OK: PE_method is a function and user_PE* is not NULL (previously caught)
# 4. OK: PE_method is ridge / xgb and user_PE* is not NULL


#' If a function, denotes a
#'   user-supplied function that should be used for computing PE. This function
#'   will be passed a data frame of covariates as its first argument and a
#'   vector of outcome values as its second argument. It must return a vector of
#'   in-sample predictions, which, if the outcome is binary or multi-class, must
#'   be maximum probability class labels. See below for examples.

my_PE_cont <- function(X, Y) {
  lm(Y ~ ., data = as.data.frame(cbind(X, Y = Y)))$fitted.values
}

my_PE_bin <- function(X, Y) {
  df <- as.data.frame(cbind(X, Y = Y))
  X <- model.matrix(Y ~ ., data = df)
  fit <- glmnet::glmnet(X, Y, family = 'binomial')
  return(predict(fit, X, s = 0.01, type = 'class'))
}

my_PE_multiclass <- function(X, Y) {
  df <- as.data.frame(cbind(X, Y = Y))
  X <- model.matrix(Y ~ ., data = df)
  fit <- glmnet::glmnet(X, Y, family = 'multinomial')
  return(predict(fit, X, s = 0.01, type = 'class'))
}

test_that("user continuous PE methods works", {
  data <- gen_data(500, 4)
  holdout <- gen_data(500, 4)
  flout <- FLAME(data, holdout, PE_method = my_PE_cont)
  expect_true(TRUE)
})

test_that("user binary factor PE methods works", {
  data <- gen_data(500, 4)
  holdout <- gen_data(500, 4)
  data$outcome <- cut(data$outcome, breaks = 2, labels = c('A', 'B'))
  holdout$outcome <- cut(holdout$outcome, breaks = 2, labels = c('A', 'B'))
  flout <- FLAME(data, holdout, PE_method = my_PE_bin)
  expect_true(TRUE)
})

test_that("user binary numeric PE methods works", {
  data <- gen_data(500, 4)
  holdout <- gen_data(500, 4)
  data$outcome <- ifelse(data$outcome > median(data$outcome), 1, 0)
  holdout$outcome <- ifelse(holdout$outcome > median(holdout$outcome), 1, 0)
  flout <- FLAME(data, holdout, PE_method = my_PE_bin)
  expect_true(TRUE)
})

test_that("user multiclass PE methods works", {
  data <- gen_data(500, 4)
  holdout <- gen_data(500, 4)
  data$outcome <- cut(data$outcome, breaks = 3, labels = c('A', 'B', 'C'))
  holdout$outcome <- cut(holdout$outcome, breaks = 3, labels = c('A', 'B', 'C'))
  flout <- FLAME(data, holdout, PE_method = my_PE_multiclass)
  expect_true(TRUE)
})

test_that("old user continuous PE methods work", {
  data <- gen_data(500, 4)
  holdout <- gen_data(500, 4)

  my_fit <- glmnet::cv.glmnet
  my_fit_params <- list(nfolds = 5)
  my_predict <- predict
  my_predict_params <- list(s = "lambda.min")

  FLAME_out <-
    FLAME(data = data, holdout = holdout,
          user_PE_fit = my_fit, user_PE_fit_params = my_fit_params,
          user_PE_predict = my_predict,
          user_PE_predict_params = my_predict_params)

  expect_true(TRUE)
})

test_that("old user binary numeric PE methods work", {
  data <- gen_data(500, 4)
  holdout <- gen_data(500, 4)
  data$outcome <- ifelse(data$outcome > median(data$outcome), 1, 0)
  holdout$outcome <- ifelse(holdout$outcome > median(holdout$outcome), 1, 0)
  my_fit <- glmnet::cv.glmnet
  my_fit_params <- list(nfolds = 5, family = 'binomial')
  my_predict <- predict
  my_predict_params <- list(s = "lambda.min", type = "class")

  FLAME_out <-
    FLAME(data = data, holdout = holdout,
          user_PE_fit = my_fit, user_PE_fit_params = my_fit_params,
          user_PE_predict = my_predict,
          user_PE_predict_params = my_predict_params)

  expect_true(TRUE)
})

test_that("old user binary factor PE methods work", {
  data <- gen_data(500, 4)
  holdout <- gen_data(500, 4)
  data$outcome <- cut(data$outcome, breaks = 2, labels = c('A', 'B'))
  holdout$outcome <- cut(holdout$outcome, breaks = 2, labels = c('A', 'B'))
  my_fit <- glmnet::cv.glmnet
  my_fit_params <- list(nfolds = 5, family = 'binomial')
  my_predict <- predict
  my_predict_params <- list(s = "lambda.min", type = "class")

  FLAME_out <-
    FLAME(data = data, holdout = holdout,
          user_PE_fit = my_fit, user_PE_fit_params = my_fit_params,
          user_PE_predict = my_predict,
          user_PE_predict_params = my_predict_params)

  expect_true(TRUE)
})

test_that("old user multiclass PE methods work", {
  data <- gen_data(500, 4)
  holdout <- gen_data(500, 4)
  data$outcome <- cut(data$outcome, breaks = 3, labels = c('A', 'B', 'C'))
  holdout$outcome <- cut(holdout$outcome, breaks = 3, labels = c('A', 'B', 'C'))
  my_fit <- glmnet::cv.glmnet
  my_fit_params <- list(nfolds = 5, family = 'multinomial')
  my_predict <- predict
  my_predict_params <- list(s = "lambda.min", type = "class")

  FLAME_out <-
    FLAME(data = data, holdout = holdout,
          user_PE_fit = my_fit, user_PE_fit_params = my_fit_params,
          user_PE_predict = my_predict,
          user_PE_predict_params = my_predict_params)

  expect_true(TRUE)
})

test_that("old-new PE method conflicts are caught", {
  data <- gen_data(500, 4)
  holdout <- gen_data(500, 4)
  my_fit <- glmnet::cv.glmnet
  my_fit_params <- list(nfolds = 5)
  my_predict <- predict
  my_predict_params <- list(s = "lambda.min")

  expect_error(FLAME(data = data, holdout = holdout,
                     PE_method = my_PE_cont,
                     user_PE_fit = my_fit, user_PE_fit_params = my_fit_params,
                     user_PE_predict = my_predict,
                     user_PE_predict_params = my_predict_params))
})
