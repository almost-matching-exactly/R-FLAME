

# NULL default for holdout set
sort_cols <-
  function(df, treated_column_name, outcome_column_name,
           type, is_missing = NULL) {

  outcome_in_data <- !is.null(df[[1]][[outcome_column_name]])

  n_covs <- ncol(df[[1]]) - 1 - outcome_in_data # Ignore treatment, outcome
  n_df <- length(df) # Always pass in a list of data frames

  # Treatment and outcome will be constant across imputations
  treatment_col <- df[[1]][treated_column_name]

  if (type == 'holdout' | (type == 'data' & outcome_in_data)) {
    outcome_col <- df[[1]][outcome_column_name]
  }

  n <- nrow(df[[1]])

  # easier way to do this
  treatment_col_ind <- which(colnames(df[[1]]) == treated_column_name)
  outcome_col_ind <- which(colnames(df[[1]]) == outcome_column_name)
  covariates <-
    which(!(1:ncol(df[[1]]) %in% c(treatment_col_ind, outcome_col_ind)))

  # For all imputed data sets
  for (i in 1:n_df) {
    tmp_df <- df[[i]]

    if (type == 'holdout' | (type == 'data' & outcome_in_data)) {
      tmp_df <- cbind(tmp_df[, covariates, drop = FALSE], outcome_col, treatment_col)
    }
    else {
      tmp_df <- cbind(tmp_df[, covariates, drop = FALSE], treatment_col)
    }

    cov_names <- colnames(tmp_df)[1:n_covs]

    # Sorting data column names
    if (type == 'holdout' | (type == 'data' & outcome_in_data)) {
      colnames(tmp_df) <- c(cov_names, 'outcome', 'treated')
    }
    else {
      colnames(tmp_df) <- c(cov_names, 'treated')
    }

    # Don't I do this twice?
    if (type == 'data') {
      for (j in 1:n_covs) {
        levels(tmp_df[, j]) <- c(levels(tmp_df[, j]), '*')
      }
    }

    # covs denotes the covariates currently being matched on
    covs <- 1:n_covs

    if (type == 'data') {
      tmp_df$matched <- rep(FALSE, n)
      tmp_df$weight <- rep(0, n)
      tmp_df$missing <- is_missing
    }

    df[[i]] <- tmp_df
  }
  return(list(df = df,
              covs = covs,
              n_covs = n_covs,
              cov_names = cov_names))
}
