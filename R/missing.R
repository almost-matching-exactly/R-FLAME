impute_missing <- function(data, n_imputations) {
  pred_mat <-
    matrix(1, nrow = ncol(data), ncol = ncol(data) - 1) %>%
    cbind(0)
  diag(pred_mat) <- 0
  mice::mice(data, m = n_imputations,
             predictorMatrix = pred_mat, printFlag = FALSE) %>%
    mice::complete(action = 'all') %>%
    return()
}

handle_missing_data <- function(data, holdout,
                                missing_data_replace, missing_holdout_replace,
                                missing_data_imputations, missing_holdout_imputations) {

  if (missing_data_replace == 0) {
    if (sum(is.na(data)) > 0) {
      stop("Found missingness in 'data' but was told to assume there was none.\n
           Please either change 'missing_data_replace' to 1, 2, or 3, or supply 'data'
           without missingness.")
    }
  }
  else if (missing_data_replace == 1) {
    data %<>%
      dplyr::drop_na()
  }
  else if (missing_data_replace == 2) {

  }
  else if (missing_data_replace == 3) {
    message('Imputing missingness in data using MICE\r')
    data <- impute_missing(data, missing_data_imputations)
    message('Finished imputation')
  }

  if (missing_holdout_replace == 0) {
    if (sum(is.na(holdout)) > 0) {
      stop("Found missingness in 'holdout' but was told to assume there was none.\n
           Please either change 'missing_holdout_replace' to 1 or 2, or supply 'holdout'
           without missingness.")
    }
  }
  else if (missing_holdout_replace == 1) {
    holdout %<>%
      dplyr::drop_na()
  }
  else if (missing_holdout_replace == 2) {
    message('Imputing missingness in holdout using MICE\r')
    holdout <- impute_missing(holdout, missing_holdout_imputations)
    message('Finished imputation')
  }

  # Change levels to allow for 'unmatched on this covariate' indicator: '*'
  if (is.data.frame(data)) {
    n_imputations <- 1
  }
  else {
    n_imputations <- length(data)
  }
  for (i in 1:n_imputations) {
    for (j in 1:(ncol(data[[i]]) - 2)) {
      levels(data[[i]][, j]) %<>% c('*')
    }
  }

  return(list(data = data,
              holdout = holdout))
}
