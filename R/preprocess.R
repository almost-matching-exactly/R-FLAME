preprocess <- function(data, holdout, C,
                       treated_column_name, outcome_column_name, n_flame_iters,
                       PE_method, user_PE_fit, user_PE_fit_params,
                       user_PE_predict, user_PE_predict_params,
                       replace, verbose, return_pe, return_bf,
                       early_stop_params,
                       missing_data, missing_holdout,
                       missing_holdout_imputations,
                       impute_with_outcome, impute_with_treatment) {
  # Get matching and holdout data, from input .csv files, if necessary
  read_data_out <-
    read_data(data, holdout, treated_column_name, outcome_column_name)

  data <- read_data_out[[1]]
  holdout <- read_data_out[[2]]

  # Was outcome supplied by user (in matching data)?
  # outcome_in_data <- !is.null(data[[outcome_column_name]])
  # Make sure the user didn't do anything funny
  check_args(data, holdout, C,
             treated_column_name, outcome_column_name, n_flame_iters,
             PE_method, user_PE_fit, user_PE_fit_params,
             user_PE_predict, user_PE_predict_params,
             replace, verbose, return_pe, return_bf,
             early_stop_params,
             missing_data, missing_holdout,
             missing_holdout_imputations,
             impute_with_outcome, impute_with_treatment)

  # Map everything to factor
  cov_inds_data <-
    which(!(colnames(data) %in% c(treated_column_name, outcome_column_name)))
  data[, cov_inds_data] <-
    lapply(data[, cov_inds_data, drop = FALSE], as.factor)

  ord <- c(cov_inds_data,
           which(colnames(data) == outcome_column_name),
           which(colnames(data) == treated_column_name))

  # Number of levels of each covariate
  n_levels <- sapply(data[, cov_inds_data, drop = FALSE], nlevels)
  if (any(n_levels > nrow(data) / 10)) {
    warning(paste("It looks like the variables {",
                  colnames(data)[cov_inds_data][n_levels > nrow(data) / 10],
                  "} might be continuous; you probably won't make many matches on them.",
                  "If this is the case, please either bin variables prior to",
                  "calling `FLAME` or do not include them at all."))
  }

  cov_inds_holdout <-
    which(!(colnames(holdout) %in% c(treated_column_name, outcome_column_name)))
  holdout[, cov_inds_holdout] <-
    lapply(holdout[, cov_inds_holdout, drop = FALSE], as.factor)

  # Map the factor levels of the covariates as supplied by the user to 0:k
  #   so that bit matching works properly.
  #   Store the mapping so you can return the data in the original format.
  remapped_data <- factor_remap(data, treated_column_name, outcome_column_name)
  data <- remapped_data$df

  mapping <- remapped_data$mapping

  # Impute missing data, if requested, else, prepare to deal with missingness
  #   as specified by missing_data
  missing_out <-
    handle_missing_data(data, holdout,
                        treated_column_name, outcome_column_name,
                        missing_data, missing_holdout,
                        missing_holdout_imputations,
                        impute_with_treatment, impute_with_outcome)

  data <- missing_out[[1]]
  holdout <- missing_out[[2]]
  is_missing <- missing_out[[3]]
  orig_missing <- missing_out[[4]]

  orig_missing[, 'col'] <- match(orig_missing[, 'col'], ord)

  # Move treatment and outcome columns to end.
  # Should maybe (probably?) be combined with factor_remap
  sort_cols_out <-
    sort_cols(data, treated_column_name, outcome_column_name,
              type = 'data', is_missing)

  data <- sort_cols_out[[1]]
  covs <- sort_cols_out[[2]]
  cov_names <- sort_cols_out[[3]]

  # Also sort holdout so that data and holdout columns correspond to one another
  holdout <-
    sort_cols(holdout, treated_column_name, outcome_column_name,
              type = 'holdout')[[1]]

  return(list(data = data, holdout = holdout, covs = covs, cov_names = cov_names,
              mapping = mapping, orig_missing = orig_missing))
}
