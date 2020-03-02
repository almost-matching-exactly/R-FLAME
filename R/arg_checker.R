check_args <-
  function(data, holdout, C, treated_column_name, outcome_column_name,
           PE_method, user_PE_fun, PE_fun_params,
           replace, verbose, want_pe, want_bf,
           early_stop_iterations, epsilon,
           early_stop_un_c_frac, early_stop_un_t_frac,
           early_stop_pe, early_stop_bf,
           missing_data, missing_holdout,
           missing_data_imputations, missing_holdout_imputations) {

  if (!is.data.frame(data)) {
    stop('data must be a data frame or a character denoting a .csv file in the working directory.')
  }
  if (!is.data.frame(holdout)) {
    stop('holdout must be a data frame, a character denoting a .csv file in the working directory,
         or a numeric proportion of data to use as a holdout set.')
  }

  data_cols <- colnames(data)
  holdout_cols <- colnames(holdout)

  if (!identical(sort(data_cols), sort(holdout_cols))) {
    stop('Data and holdout must contain identical column names.')
  }

  if (!is.numeric(C) | C < 0 | is.infinite(C)) {
    stop('C must be a finite, nonnegative scalar.')
  }

  if (!is.character(treated_column_name)) {
    stop('If you specify treated_column_name, it must be a character.')
  }

  if (!(treated_column_name %in% data_cols)) {
    stop('treated_column_name must be the name of a column in data.')
  }

  if (!(treated_column_name %in% holdout_cols)) {
    stop('treated_column_name must be the name of a column in holdout.')
  }

  if (!is.character(outcome_column_name)) {
    stop('If you specify outcome_column_name, it must be a character.')
  }

  if (!(outcome_column_name %in% data_cols)) {
    stop('outcome_column_name must be the name of a column in data.')
  }

  if (!(outcome_column_name %in% holdout_cols)) {
    stop('outcome_column_name must be the name of a column in holdout.')
  }

  if (!(PE_method %in% c('elasticnet', 'xgb'))) {
    stop("PE_method must be one of 'elasticnet' or 'xgb'.
         To supply your own model to fit, use user_PE_fit.")
  }

  if (!is.logical(replace)) {
    stop('replace must be a logical scalar')
  }

  if (!(verbose %in% c(0, 1, 2, 3))) {
    stop('Verbose must be one of: 0, 1, 2, 3.')
  }

  if (!is.logical(want_pe)) {
    stop('want_pe must be a logical scalar')
  }

  if (!is.logical(want_bf)) {
    stop('want_bf must be a logical scalar')
  }

  ## Early stop parameters
  if (!is.numeric(early_stop_iterations) | early_stop_iterations < 0) {
    stop('early_stop_iterations must be a nonnegative scalar')
  }

  if (!is.numeric(epsilon) | early_stop_iterations <= 0) {
    stop('epsilon must be a positive scalar')
  }

  if (!is.numeric(early_stop_un_c_frac) |
      early_stop_un_c_frac < 0 |
      early_stop_un_c_frac > 1) {
    stop('early_stop_un_c_frac must be a fraction between 0 and 1 (inclusive).')
  }

  if (!is.numeric(early_stop_un_t_frac) |
      early_stop_un_t_frac < 0 |
      early_stop_un_t_frac > 1) {
    stop('early_stop_un_t_frac must be a fraction between 0 and 1 (inclusive).')
  }

  if (!is.numeric(early_stop_pe) | early_stop_pe < 0) {
    stop('early_stop_pe must be a nonnegative scalar')
  }

  if (!is.numeric(early_stop_bf) | early_stop_bf < 0 | early_stop_bf > 2) {
    stop('early_stop_bf must be a scalar between 0 and 2 (inclusive)')
  }

  ## Missing data parameters
  if (!is.numeric(missing_data) | !(missing_data %in% c(0, 1, 2, 3))) {
    stop('missing_data must be one of: 0, 1, 2, 3')
  }
  if (!is.numeric(missing_holdout) | !(missing_holdout %in% c(0, 1, 2))) {
    stop('missing_data must be one of: 0, 1, 2')
  }
  if (!is.numeric(missing_data_imputations) | missing_data_imputations < 1) {
    stop('missing_data_imputations must be an integer greater than 1')
  }
  if (!is.numeric(missing_holdout_imputations) | missing_holdout_imputations < 1) {
    stop('missing_holdout_imputations must be an integer greater than 1')
  }
}
