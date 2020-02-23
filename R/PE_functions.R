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

predict_elasticnet <- function(holdout, covs, cov_to_drop, alpha) {
  c(X_treat, X_control, Y_treat, Y_control) %<-%
    setup_preds(holdout, covs, cov_to_drop)
  # browser()
  MSE_treated <-
    tryCatch(
      error = function(cnd) 0,
      glmnet(X_treat,
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
      glmnet(X_control,
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
