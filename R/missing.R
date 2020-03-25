impute_missing <- function(data, outcome_in_data, n_imputations,
                           treated_column_name, outcome_column_name,
                           impute_with_treatment, impute_with_outcome) {

  treatment_ind <- which(colnames(data) == treated_column_name)
  outcome_ind <- which(colnames(data) == outcome_column_name)

  pred_mat <- matrix(1, nrow = ncol(data), ncol = ncol(data))
  diag(pred_mat) <- 0

  if (!impute_with_treatment) {
    pred_mat[, treatment_ind] <- 0
  }
  if (!impute_with_outcome) {
    pred_mat[, outcome_ind] <- 0
  }

  pred_mat[c(treatment_ind, outcome_ind), ] <- 0

  mice::mice(data, m = n_imputations,
             predictorMatrix = pred_mat, printFlag = FALSE) %>%
    mice::complete(action = 'all') %>%
    return()
}

handle_missing_data <-
  function(data, holdout, outcome_in_data,
           treated_column_name, outcome_column_name,
           missing_data, missing_holdout,
           missing_data_imputations, missing_holdout_imputations,
           impute_with_treatment, impute_with_outcome) {

  treatment_ind_data <- which(colnames(data) == treated_column_name)
  treatment_ind_holdout <- which(colnames(holdout) == treated_column_name)
  treatment_data <- data[[treatment_ind_data]]
  treatment_holdout <- holdout[[treatment_ind_holdout]]

  # outcome_ind <- which(colnames(data) == outcome_column_name)

  # treatment_data <- dplyr::pull(data, !!rlang::enquo(treated_column_name))
  # treatment_holdout <- dplyr::pull(holdout, !!rlang::enquo(treated_column_name))

  if (any(is.na(treatment_data)) | any(is.na(treatment_holdout))) {
    message('Found missingness in the treatment. Corresponding rows will automatically be dropped.')
    # Do we drop these?
  }

  outcome_ind_holdout <- which(colnames(holdout) == outcome_column_name)
  outcome_holdout <- holdout[[outcome_ind_holdout]]
  # outcome_holdout <- dplyr::pull(holdout, !!rlang::enquo(outcome_column_name))

  if (outcome_in_data) {
    outcome_ind_data <- which(colnames(data) == outcome_column_name)
    outcome_data <- data[[outcome_ind_data]]
    # outcome_data <- dplyr::pull(data, !!rlang::enquo(outcome_column_name))
    if (all(is.na(outcome_data))) {
      stop('Outcome in data is entirely missing. If you do not have an outcome, do not supply a corresponding column.')
    }
    if (all(is.na(outcome_holdout))) {
      stop('Outcome in holdout is entirely missing.')
    }
    if (any(is.na(outcome_data)) | any(is.na(outcome_holdout))) {
      message('Found missingness in the outcome. Corresponding rows will automatically be dropped.')
    }
    to_drop_data <- is.na(outcome_data) | is.na(treatment_data)
  }
  else {
    if (any(is.na(outcome_holdout))) {
      message('Found missingness in the outcome. Corresponding rows will automatically be dropped.')
    }
    to_drop_data <- is.na(treatment_data)
  }

  to_drop_holdout <- is.na(outcome_holdout) | is.na(treatment_holdout)

  data <- data[!to_drop_data, ]
  holdout <- holdout[!to_drop_holdout, ]

  if (missing_data == 0) {
    is_missing <- FALSE
    if (sum(is.na(data)) > 0) {
      stop("Found missingness in 'data' but was told to assume there was none.\n
           Please either change 'missing_data' to 1, 2, or 3, or supply 'data'
           without missingness.")
    }
  }
  else if (missing_data == 1) {
    is_missing <- apply(data, 1, function(row) any(is.na(row)))
    if (all(is_missing)) {
      stop('All rows in data contain missingness. In this case, matches may only be made if missing_data = 2 or missing_data = 3')
    }
  }
  else if (missing_data == 2) {
    is_missing <- FALSE
    if (sum(is.na(data)) > 0) {
      data <- impute_missing(data, outcome_in_data, missing_data_imputations,
                             treated_column_name, outcome_column_name,
                             impute_with_treatment, impute_with_outcome)

    }
    else {
      message('No missing data found; skipping imputation.')
    }
  }
  else if (missing_data == 3) {
    is_missing <- FALSE
    if (sum(is.na(data)) > 0) {
      if (outcome_in_data) {
        cov_inds <- setdiff(1:ncol(data), c(treatment_ind_data, outcome_ind_data))
      }
      else {
        cov_inds <- setdiff(1:ncol(data), treatment_ind_data)
      }

      tmp_data <- data

      for (cov in cov_inds) {
        tmp_data[[cov]] %<>% as.numeric()
        which_missing <- is.na(tmp_data[[cov]])
        if (sum(which_missing) > 0) {
          tmp_data[[cov]][which_missing] <-
            max(as.numeric(tmp_data[[cov]]), na.rm = TRUE) %>%
            magrittr::add(seq_len(sum(which_missing)))
        }
        tmp_data[[cov]] %<>% as.factor()
      }

      data <- tmp_data
    }
    else {
      warning('Was directed to skip matches on missing values, but no missing values found.')
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
      holdout <- impute_missing(holdout, outcome_in_data, missing_holdout_imputations,
                                treated_column_name, outcome_column_name,
                                impute_with_treatment, impute_with_outcome)
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

  if (outcome_in_data) {
    covariates <-
      which(!(1:ncol(data[[1]]) %in% c(treatment_ind_data, outcome_ind_data)))
  }
  else {
    covariates <-
      which(!(1:ncol(data[[1]]) %in% c(treatment_ind_data)))
  }

  for (i in 1:n_imputations) {
    # for (j in 1:(ncol(data[[i]]) - 1 - outcome_in_data)) {
    #   levels(data[[i]][, j]) %<>% c('*')
    # }
    for (j in covariates) {
      levels(data[[i]][, j]) %<>% c('*')
    }
  }

  return(list(data = data,
              holdout = holdout,
              is_missing = is_missing))
}
