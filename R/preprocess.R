preprocess <- function(data, holdout, C, algo, weights,
                       treated_column_name, outcome_column_name, n_flame_iters,
                       PE_method, user_PE_fit, user_PE_fit_params,
                       user_PE_predict, user_PE_predict_params,
                       replace, estimate_CATEs, verbose, return_pe, return_bf,
                       early_stop_params,
                       missing_data, missing_holdout,
                       missing_holdout_imputations,
                       impute_with_outcome, impute_with_treatment) {

  # Get matching and holdout data, from input .csv files, if necessary
  read_data_out <-
    read_data(data, holdout, treated_column_name, outcome_column_name, weights)

  data <- read_data_out[[1]]
  holdout <- read_data_out[[2]]

  if (is.null(outcome_column_name) || is.null(data[[outcome_column_name]])) {
    outcome_type <- 'none'
  }
  else if (length(unique(data[[outcome_column_name]])) == 2) {
    outcome_type <- 'binary'
  }
  else if (is.factor(data[[outcome_column_name]])) {
    outcome_type <- 'categorical'
  }
  else {
    outcome_type <- 'continuous'
  }

  info <- list('algo' = algo,
               'treatment' = treated_column_name,
               'outcome' = outcome_column_name,
               'replacement' = replace,
               'estimate_CATEs' = estimate_CATEs,
               'missing_data' = missing_data,
               'missing_holdout' = missing_holdout,
               'outcome_type' = outcome_type)

  # Make sure the user didn't do anything funny
  check_args(data, holdout, C, weights,
             n_flame_iters,
             PE_method, user_PE_fit, user_PE_fit_params,
             user_PE_predict, user_PE_predict_params,
             verbose, return_pe, return_bf,
             early_stop_params,
             missing_holdout_imputations,
             impute_with_outcome, impute_with_treatment, info)

  # Map covariates to factor
  cov_inds_data <-
    which(!(colnames(data) %in% c(treated_column_name, outcome_column_name)))

  data[, cov_inds_data] <-
    lapply(data[, cov_inds_data, drop = FALSE], as.factor)

  # Drop extraneous levels, warning user if so
  data <- warn_level_drop(data, cov_inds_data, type = 'data')

  ord <- c(cov_inds_data,
           which(colnames(data) == outcome_column_name),
           which(colnames(data) == treated_column_name))

  # Number of levels of each covariate
  n_levels <- vapply(data[, cov_inds_data, drop = FALSE], nlevels, numeric(1))
  if (any(n_levels > nrow(data) / 10)) {
    warning(paste0("It looks like the variables { ",
                  colnames(data)[cov_inds_data][n_levels > nrow(data) / 10],
                  " } might be continuous; you probably won't make (m)any ",
                  "matches on them. If this is the case, please either bin the",
                  " variables prior to calling `", algo, "` or do not include ",
                  "them at all.\n"), call. = FALSE)
  }

  # Might be different from data if no outcome
  cov_inds_holdout <-
    which(!(colnames(holdout) %in% c(treated_column_name, outcome_column_name)))
  holdout[, cov_inds_holdout] <-
    lapply(holdout[, cov_inds_holdout, drop = FALSE], as.factor)

  holdout <- warn_level_drop(holdout, cov_inds_holdout, type = 'holdout')

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

  # Sort holdout so that data and holdout columns correspond to one another.
  holdout <-
    sort_cols(holdout, treated_column_name, outcome_column_name,
              type = 'holdout')[[1]]

  return(list(data = data, holdout = holdout, covs = covs,
              cov_names = cov_names, mapping = mapping,
              orig_missing = orig_missing, info = info))
}
