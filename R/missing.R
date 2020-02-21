impute_missing <- function(data, n_imputations) {
  # Ignore outcome when imputing missing data?
  mice::mice(data, m = n_imputations, printFlag = FALSE) %>%
    mice::complete() %>%
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
    data <- impute_missing(data, missing_data_imputations)
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
    holdout <- impute_missing(holdout, missing_holdout_imputations)
  }

  return(list(data = data,
              holdout = holdout))
}
