# Some hard-coded parameter values to cross-validate XGBoost over
#   If the user cares about these they'll just input their own PE function.
cv_xgboost <- function(X, Y, obj) {
  # Return best 5-fold XGBoost fit for Y ~ X across various parameter
  #  configurations
  eta <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5)
  max_depth <- c(2, 3, 4, 6, 8)
  alpha <- c(0.01, 0.1, 0.5, 1, 5)
  nrounds <- c(5, 10, 50, 200, 500)
  subsample <- c(0.1, 0.3, 0.5, 0.75, 1)


  if (is.factor(Y)) {
    Y <- as.numeric(Y) - 1 # -1 for XGBoost compatibility
  }

  param_combs <- expand.grid(eta, max_depth, alpha, nrounds, subsample)

  colnames(param_combs) <-
    c('eta', 'max_depth', 'alpha', 'nrounds', 'subsample')

  error <- vector(mode = 'numeric', length = length(param_combs))

  for (i in seq_along(param_combs)) {
    params <- list(objective = obj,
                   eta = param_combs$eta[i],
                   max_depth = param_combs$max_depth[i],
                   alpha = param_combs$alpha[i],
                   subsample = param_combs$subsample[i])
    if (obj == 'multi:softmax') {
      params <- c(params, list(num_class = length(unique(Y))))
    }

    cv <-
      xgboost::xgb.cv(data = X,
                      label = Y,
                      params = params,
                      nrounds = param_combs$nrounds[i],
                      nfold = 5, verbose = 0)

    error[i] <- cv$evaluation_log[param_combs$nrounds[i], 4]
  }

  best_params <- param_combs[which.min(error), ]
  params <- c(best_params, list(objective = obj))
  if (obj == 'multi:softmax') {
    params <- c(params, list(num_class = length(unique(Y))))
  }
  params$nrounds <- NULL
  fit <- xgboost::xgboost(data = X,
                  label = Y,
                  params = params,
                  nrounds = best_params$nrounds,
                  verbose = 0)

  return(fit)
}

setup_preds <- function(holdout, covs, covs_to_drop) {

  # Covariates for prediction
  cov_set <- setdiff(covs, covs_to_drop)

  # Fit treated, control separately
  Y_treat <- holdout$outcome[holdout$treated == 1]
  Y_control <- holdout$outcome[holdout$treated == 0]

  X_treat <- holdout[holdout$treated == 1, cov_set, drop = FALSE]
  X_control <- holdout[holdout$treated == 0, cov_set, drop = FALSE]

  return(list(X_treat = X_treat,
              X_control = X_control,
              Y_treat = Y_treat,
              Y_control = Y_control))
}

get_error <- function(X, Y, PE_method,
                      fit_fun, predict_fun, fit_params, predict_params,
                      user_fit_predict, user_backpassed) {

  if (!is.factor(Y) & length(unique(Y)) > 2) {
    outcome_type <- 'continuous'
  }
  else {
    outcome_type <- 'discrete'
  }

  # Only in the case that the user specifies nothing
  if (is.null(user_fit_predict)) {
    X <- model.matrix(~ ., data = X)
    fit <- do.call(fit_fun, c(list(X, Y), fit_params))
    if (PE_method == 'ridge') {
      if (outcome_type == 'continuous') {
        preds <- predict(fit, X)
      }
      else {
        preds <- predict(fit, X, type = 'class')
      }
    }
    else if (PE_method == 'xgb') {
      if (is.factor(Y)) {
        if (length(unique(Y)) == 2) {
          preds <- levels(Y)[(predict(fit, X) > 0.5) + 1]
        }
        else {
          preds <- levels(Y)[predict(fit, X) + 1]
        }
      }
      else {
        if (length(unique(Y)) == 2) {
          preds <- predict(fit, X) > 0.5
        }
        else {
          preds <- predict(fit, X)
        }
      }
    }
    else {
      stop('Unimplemented PE_method')
    }
  }
  else {
    if (user_backpassed) {
      X <- model.matrix(~ ., data = X)
    }

    preds <- user_fit_predict(X, Y)
  }

  if (outcome_type != 'continuous') {
    if (length(unique(preds)) > length(unique(Y))) {
      warning('It looks like your function for computing PE ',
              'does not return predicted class labels. If so,
              the PE computation will not work as expected.', call. = FALSE)
    }
    error <- mean(preds != Y)
  }
  else {
    error <- mean((preds - Y) ^ 2)
  }

  return(error)
}

predict_master <-
  function(holdout, covs, covs_to_drop, PE_method,
           PE_fit, PE_predict, PE_fit_params, PE_predict_params,
           user_fit_predict, user_backpassed) {

  n_imputations <- length(holdout) # List of data frames

  PE <- vector('numeric', length = n_imputations)

  for (i in 1:n_imputations) {
    setup_out <- setup_preds(holdout[[i]], covs, covs_to_drop)
    X_treat <- setup_out[[1]]
    X_control <- setup_out[[2]]
    Y_treat <- setup_out[[3]]
    Y_control <- setup_out[[4]]

    error_treat <-
      get_error(X_treat, Y_treat, PE_method,
               PE_fit, PE_predict, PE_fit_params, PE_predict_params,
               user_fit_predict, user_backpassed)

    error_control <-
      get_error(X_control, Y_control, PE_method,
                PE_fit, PE_predict, PE_fit_params, PE_predict_params,
                user_fit_predict, user_backpassed)

    PE[i] <- error_treat + error_control
  }
  return(mean(PE))
  }

get_PE <- function(covs_to_drop, covs, holdout, PE_method,
                   user_PE_fit, user_PE_fit_params,
                   user_PE_predict, user_PE_predict_params) {

  user_fit_predict <- NULL
  user_backpassed <- FALSE # For back-compatibility; to be removed in future

  PE_predict_params <- list()

  if (!is.null(user_PE_fit)) {
    PE_fit <- user_PE_fit
    PE_fit_params <- user_PE_fit_params

    user_backpassed <- TRUE
  }
  else {
    if (is.function(PE_method)) {
      user_fit_predict <- PE_method
    }
    else {
      if (PE_method == 'ridge') {
        PE_fit <- glmnet::cv.glmnet
        if (length(unique(holdout[[1]]$outcome)) == 2) {
          family <- 'binomial'
          PE_predict_params <- list(type = 'class')
        }
        else if (is.factor(holdout[[1]]$outcome)) {
          family <- 'multinomial'
          PE_predict_params <- list(type = 'class')
        }
        else {
          family <- 'gaussian'
        }
        PE_fit_params <- list(family = family, nfolds = 5)
      }
      else if (PE_method == 'xgb') {
        PE_fit <- cv_xgboost
        if (length(unique(holdout[[1]]$outcome)) == 2) {
          obj <- 'binary:logistic'
        }
        else if (is.factor(holdout[[1]]$outcome)) {
          obj <- 'multi:softmax'
        }
        else {
          obj <- 'reg:squarederror'
        }
        PE_fit_params <- list(obj = obj)
      }
      # Else caught by arg_checker
    }
  }

  if (!is.null(user_PE_predict)) {
    user_backpassed <- TRUE
    PE_predict <- user_PE_predict
    PE_predict_params <- user_PE_predict_params
  }
  else {
    PE_predict <- predict
    PE_predict_params <- list()
  }

  # Purely for back-compatibility
  if (user_backpassed) {
    user_fit_predict <- function(X, Y) {
      fit <- do.call(PE_fit, c(list(X, Y), PE_fit_params))
      preds <- do.call(PE_predict, c(list(fit, X), PE_predict_params))
    }
  }

  PE <- predict_master(holdout, covs, covs_to_drop, PE_method,
                       PE_fit, PE_predict, PE_fit_params, PE_predict_params,
                       user_fit_predict, user_backpassed)
  return(PE)
}
