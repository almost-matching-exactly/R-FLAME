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

  mice_out <-
    mice::mice(data, m = n_imputations,
               predictorMatrix = pred_mat, printFlag = FALSE)

  imputed_data <- mice::complete(mice_out, action = 'all')

  return(imputed_data)
}

handle_missing_data <-
  function(data, holdout,
           treated_column_name, outcome_column_name,
           missing_data, missing_holdout,
           missing_holdout_imputations,
           impute_with_treatment, impute_with_outcome) {

    outcome_in_data <- !is.null(data[[outcome_column_name]])

    # Corresponds to data only
    cov_inds <- which(!(colnames(data) %in%
                        c(treated_column_name, outcome_column_name)))

    if (outcome_in_data) {
      to_drop_data <- is.na(data[[outcome_column_name]]) |
                      is.na(data[[treated_column_name]])
    }
    else {
      to_drop_data <- is.na(data[[treated_column_name]])
    }

    to_drop_holdout <- is.na(holdout[[outcome_column_name]]) |
                       is.na(holdout[[treated_column_name]])

    if (any(to_drop_data)) {
      message('Found missingness in `data` in treatment and/or outcome; ',
              'corresponding rows will be dropped.')
    }
    if (any(to_drop_holdout)) {
      message('Found missingness in `holdout` in treatment and/or outcome; ',
              'corresponding rows will be dropped.')
    }

    if (all(to_drop_data)) {
      stop('Dropping all rows in `data` due to missingness ',
           'in treatment and/or outcome.')
    }
    if (all(to_drop_holdout)) {
      stop('Dropping all rows in `holdout` due to missingness ',
           'in treatment and/or outcome.')
    }

    data <- data[!to_drop_data, ]
    holdout <- holdout[!to_drop_holdout, ]

    orig_missing <- which(is.na(data), arr.ind = TRUE)

    if (missing_data == 'none') {
      is_missing <- FALSE
      if (sum(is.na(data)) > 0) {
        stop('Found missingness in `data` but was told to assume there ',
             'was none. Please either change `missing_data` or ',
             'supply `data` without missingness.')
      }
    }
    else if (missing_data == 'drop') {
      is_missing <- apply(data, 1, function(row) any(is.na(row)))
      if (all(is_missing)) {
        stop('All rows in `data` contain missingness. ',
             'In this case, matches may only be made if `missing_data` ',
             " = 'keep' or `missing_data` = 'impute'.")
      }
    }
    else if (missing_data == 'impute') {
      is_missing <- FALSE
      if (sum(is.na(data)) > 0) {
        message('Starting imputation of `data`\r', appendLF = FALSE)
        data <- impute_missing(data, outcome_in_data, 1,
                               treated_column_name, outcome_column_name,
                               impute_with_treatment, impute_with_outcome)
        message('Finished imputation of `data`\r', appendLF = FALSE)
      }
      else {
        message('No missing data found; skipping imputation.')
      }
    }
    else if (missing_data == 'keep') {
      is_missing <- FALSE
      if (sum(is.na(data)) > 0) {
        # Substitute large, unique values for missing values, preventing
        # matches from being made on them
        tmp_data <- data
        for (cov in cov_inds) {
          # -1 for conversion from factor to numeric
          max_val <- max(as.numeric(tmp_data[[cov]]), na.rm = TRUE) - 1
          which_missing <- is.na(tmp_data[[cov]])
          n_missing <- sum(which_missing)
          if (sum(which_missing) > 0) {
            new_levels <-
              c(levels(tmp_data[[cov]]), max_val + seq_len(n_missing))

            tmp_data[[cov]] <- as.numeric(tmp_data[[cov]]) - 1

            tmp_data[[cov]][which_missing] <-
              max_val + seq_len(sum(which_missing))

            tmp_data[[cov]] <- factor(tmp_data[[cov]], levels = new_levels)
          }
        }

        data <- tmp_data
      }
      else {
        warning('Was directed to skip matches on missing values, ',
                'but no missing values found.')
      }
    }

    if (missing_holdout == 'none') {
      if (sum(is.na(holdout)) > 0) {
        stop('Found missingness in `holdout` but was told to assume ',
             'there was none. Please either change `missing_holdout` or ',
             'supply `holdout` without missingness.')
      }
    }
    else if (missing_holdout == 'drop') {
      holdout <- holdout[complete.cases(holdout), ]
    }
    else if (missing_holdout == 'impute') {
      if (sum(is.na(holdout)) > 0) {
        message('Starting imputation of `holdout`\r', appendLF = FALSE)
        holdout <-
          impute_missing(holdout, outcome_in_data, missing_holdout_imputations,
                         treated_column_name, outcome_column_name,
                         impute_with_treatment, impute_with_outcome)
        message('Finished imputation of `holdout`\r', appendLF = FALSE)
      }
    }

    if (is.data.frame(data)) {
      data <- list(data)
    }

    if (is.data.frame(holdout)) {
      holdout <- list(holdout)
    }

    return(list(data = data,
                holdout = holdout,
                is_missing = is_missing,
                orig_missing = orig_missing))
  }
