impute_missing <- function(data, n_imputations) {
  # browser()
  # pred_mat <-
  #   matrix(1, nrow = ncol(data), ncol = ncol(data)) %>%
  #   cbind(0)
  # diag(pred_mat) <- 0
  mice::mice(data, m = n_imputations, printFlag = FALSE) %>%
    mice::complete(action = 'all') %>%
    return()
}

handle_missing_data <- function(data, holdout,
                                missing_data, missing_holdout,
                                n_data_imputations, n_holdout_imputations) {

  if (missing_data == 0) {
    if (sum(is.na(data)) > 0) {
      stop("Found missingness in 'data' but was told to assume there was none.\n
           Please either change 'missing_data' to 1, 2, or 3, or supply 'data'
           without missingness.")
    }
  }
  else if (missing_data == 1) {
    data %<>%
      tidyr::drop_na()
  }
  else if (missing_data == 2) {
    # Replace the missing values with unique integers, each larger than all
    # observed covariate values in data
    replace_vals <-
      max(data[, 1:(ncol(data) - 2)], na.rm = TRUE) + seq_len(sum(is.na(data)))
    data[which(is.na(data), arr.ind = TRUE)] <- replace_vals

  }
  else if (missing_data == 3) {
    if (sum(is.na(data)) > 0) {
      data <- impute_missing(data, n_data_imputations)
    }
  }

  if (missing_holdout == 0) {
    if (sum(is.na(holdout)) > 0) {
      stop("Found missingness in 'holdout' but was told to assume there was none.\n
           Please either change 'missing_holdout' to 1 or 2, or supply 'holdout'
           without missingness.")
    }
  }
  else if (missing_holdout == 1) {
    holdout %<>%
      tidyr::drop_na()
  }
  else if (missing_holdout == 2) {
    if (sum(is.na(holdout)) > 0) {
      holdout <- impute_missing(holdout, n_holdout_imputations)
    }
  }

  if (is.data.frame(holdout)) {
    holdout %<>% list()
  }

  # Change levels to allow for 'unmatched on this covariate' indicator: '*'
  if (is.data.frame(data)) {
    n_imputations <- 1
    data %<>% list()
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
