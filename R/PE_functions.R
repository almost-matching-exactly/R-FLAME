setup_preds <- function(holdout, cov_to_drop) {
  Y_treat <-
    holdout %>%
    filter(treated == 1) %>%
    pull(outcome)

  Y_control <-
    holdout %>%
    filter(treated == 0) %>%
    pull(outcome)

  covs_treat <-
    holdout %>%
    dplyr::filter(treated == 1) %>%
    dplyr::select(-c(cov_to_drop, treated))
  X_treat <- model.matrix(outcome ~ ., covs_treat)

  covs_control <-
    holdout %>%
    dplyr::filter(treated == 0) %>%
    dplyr::select(-c(cov_to_drop, treated))
  X_control <- model.matrix(outcome ~ ., covs_control)

  return(list(X_treat = X_treat,
              X_control = X_control,
              Y_treat = Y_treat,
              Y_control = Y_control))
}

predict_elasticnet <- function(holdout, cov_to_drop, alpha) {
  c(X_treat, X_control, Y_treat, Y_control) %<-%
    setup_preds(holdout, cov_to_drop)

  MSE_treated <-
    tryCatch(
      error = function(cnd) 0,
      glmnet(X_treat,
             Y_treat,
             alpha = alpha) %>%
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
             alpha = alpha) %>%
        predict(X_control) %>%
        magrittr::subtract(Y_control) %>%
        magrittr::raise_to_power(2) %>%
        mean()
    )
 return(MSE_treated + MSE_control)
}
