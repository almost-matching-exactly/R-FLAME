setup_preds <- function(holdout, covs, cov_to_drop) {
  covs_to_test <- setdiff(covs, cov_to_drop)
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
    dplyr::select(covs_to_test, outcome)

  X_treat <- model.matrix(outcome ~ ., covs_treat)

  covs_control <-
    holdout %>%
    dplyr::filter(treated == 0) %>%
    dplyr::select(covs_to_test, outcome)

  X_control <- model.matrix(outcome ~ ., covs_control)

  return(list(X_treat = X_treat,
              X_control = X_control,
              Y_treat = Y_treat,
              Y_control = Y_control))
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

  n_imputations <- length(holdout) # It's a list of dataframes

  PE <- vector(mode = 'numeric', length = n_imputations)
  for (i in 1:n_imputations) {
    # c(X_treat, X_control, Y_treat, Y_control) %<-%
      # setup_preds(holdout[[i]], covs, cov_to_drop)

    setup_out <- setup_preds(holdout[[i]], covs, cov_to_drop)
    X_treat <- setup_out[[1]]
    X_control <- setup_out[[2]]
    Y_treat <- setup_out[[3]]
    Y_control <- setup_out[[4]]

    browser()
    MSE_treat <- do.call(get_MSE,
                         c(list(X = X_treat, Y = Y_treat, func = PE_func),
                           PE_func_params))
    # browser()
    MSE_control <- do.call(get_MSE,
                         c(list(X = X_control, Y = Y_control, func = PE_func),
                           PE_func_params))
    PE[i] <- (MSE_treat + MSE_control) / 2
  }
  return(mean(PE))
}
