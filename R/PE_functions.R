predict_elasticnet <- function(holdout, cov_to_drop, alpha) {
  holdout_ctl <-
    holdout %>%
    filter(treated == 0) %>%
    dplyr::select(-c(cov_to_drop, treated))

  holdout_trt <-
    holdout %>%
    filter(treated == 1) %>%
    dplyr::select(-c(cov_to_drop, treated))

  MSE_treated <-
    tryCatch(
      error = function(cnd) 0,
      glmnet(dplyr::select(holdout_trt, -outcome),
             holdout_trt$outcome,
             alpha = alpha) %>%
      predict(dplyr::select(holdout_trt, -outcome)) %>%
      magrittr::subtract(holdout_trt$outcome) %>%
      magrittr::raise_to_power(2) %>%
      mean()
    )
  MSE_control <-
    tryCatch(
      error = function(cnd) 0,
      glmnet(dplyr::select(holdout_trt, -outcome),
             holdout_trt$outcome,
             alpha = alpha) %>%
      predict(dplyr::select(holdout_trt, -outcome)) %>%
      magrittr::subtract(holdout_trt$outcome) %>%
      magrittr::raise_to_power(2) %>%
      mean()
    )
 return(MSE_treated + MSE_control)
}
