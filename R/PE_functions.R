setup_preds <- function(holdout, covs, cov_to_drop) {

  Y_treat <-
    holdout %>%
    dplyr::filter(treated == 1) %>%
    dplyr::pull(outcome)

  Y_control <-
    holdout %>%
    dplyr::filter(treated == 0) %>%
    dplyr::pull(outcome)

  covs_treat <-
    holdout %>%
    dplyr::filter(treated == 1) %>%
    dplyr::select(setdiff(covs, cov_to_drop), outcome)
    # dplyr::select(-c(cov_to_drop, treated))
  X_treat <- model.matrix(outcome ~ ., covs_treat)

  covs_control <-
    holdout %>%
    dplyr::filter(treated == 0) %>%
    dplyr::select(setdiff(covs, cov_to_drop), outcome)
    # dplyr::select(-c(cov_to_drop, treated))
  X_control <- model.matrix(outcome ~ ., covs_control)
  # browser()
  return(list(X_treat = X_treat,
              X_control = X_control,
              Y_treat = Y_treat,
              Y_control = Y_control))
}

predict_custom <- function(holdout, covs, cov_to_drop) {
  c(X_treat, X_control, Y_treat, Y_control) %<-%
    setup_preds(holdout, covs, cov_to_drop)

}

get_MSE <- function(X, Y, func, ...) {
  MSE <-
    tryCatch(
    error = function(cnd) 0,
    func(X, Y, ...) %>%
      predict(X) %>%
      magrittr::subtract(Y) %>%
      magrittr::raise_to_power(2) %>%
      mean()
    )
  return(MSE)
}

predict_master <- function(holdout, covs, cov_to_drop, PE_func, PE_func_params) {
  c(X_treat, X_control, Y_treat, Y_control) %<-%
    setup_preds(holdout, covs, cov_to_drop)
  # browser()
  MSE_treat <- do.call(get_MSE,
                       c(list(X = X_treat, Y = Y_treat, func = PE_func),
                         PE_func_params))
  # browser()
  MSE_control <- do.call(get_MSE,
                       c(list(X = X_control, Y = Y_control, func = PE_func),
                         PE_func_params))
  return(MSE_treat + MSE_control)
}

predict_xgb <- function(holdout, covs, cov_to_drop) {
  require(xgboost)
  c(X_treat, X_control, Y_treat, Y_control) %<-%
    setup_preds(holdout, covs, cov_to_drop)

  MSE_treated <-
    tryCatch(
      error = function(cnd) 0,
      xgboost::xgboost(data = X_treat, label = Y_treat, verbose = 0, nrounds = 100) %>%
      predict(X_treat) %>%
      magrittr::subtract(Y_treat) %>%
      magrittr::raise_to_power(2) %>%
      mean()
    )

  MSE_control <-
    tryCatch(
      error = function(cnd) 0,
      xgboost::xgboost(data = X_control, label = Y_control, verbose = 0, nrounds = 100) %>%
        predict(X_treat) %>%
        magrittr::subtract(Y_treat) %>%
        magrittr::raise_to_power(2) %>%
        mean()
    )
  return(MSE_treated + MSE_control)
}

predict_elasticnet <- function(holdout, covs, cov_to_drop, alpha) {
  c(X_treat, X_control, Y_treat, Y_control) %<-%
    setup_preds(holdout, covs, cov_to_drop)
  # browser()
  MSE_treated <-
    tryCatch(
      error = function(cnd) 0,
      glmnet::glmnet(X_treat,
             Y_treat,
             alpha = 0,
             lambda = 0.2 / nrow(X_treat)) %>%
      predict(X_treat) %>%
      magrittr::subtract(Y_treat) %>%
      magrittr::raise_to_power(2) %>%
      mean()
    )

  MSE_control <-
    tryCatch(
      error = function(cnd) 0,
      glmnet::glmnet(X_control,
             Y_control,
             alpha = 0,
             lambda = 0.2 / nrow(X_control)) %>%
        predict(X_control) %>%
        magrittr::subtract(Y_control) %>%
        magrittr::raise_to_power(2) %>%
        mean()
    )
 return(MSE_treated + MSE_control)
}
