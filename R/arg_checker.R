check_args <-
  function(data, holdout, outcome_in_data, C,
           treated_column_name, outcome_column_name,
           binning_method,
           PE_method, user_PE_fit, user_PE_fit_params,
           user_PE_predict, user_PE_predict_params,
           replace, verbose, want_pe, want_bf,
           early_stop_iterations, early_stop_epsilon,
           early_stop_un_c_frac, early_stop_un_t_frac,
           early_stop_pe, early_stop_bf,
           missing_data, missing_holdout,
           missing_data_imputations, missing_holdout_imputations,
           impute_with_outcome, impute_with_treatment) {

  if (!is.data.frame(data)) {
    stop('data must be a data frame or a character denoting a .csv file in the working directory.')
  }
  if (!is.data.frame(holdout)) {
    stop('holdout must be a data frame, a character denoting a .csv file in the working directory,
         or a numeric proportion of data to use as a holdout set.')
  }

  data_cols <- colnames(data)
  holdout_cols <- colnames(holdout)

  if (!(outcome_column_name %in% holdout_cols)) { # I think I can kill this
    stop('Holdout must contain outcome column with name outcome_column_name')
  }

  if (!outcome_in_data) {
    if (!identical(data_cols,
                   holdout_cols[-match(outcome_column_name, holdout_cols)])) {
      stop('Non-outcome columns of data and holdout must have identical names.')
    }
  }
  else {
    if (!identical(data_cols, holdout_cols)) {
      stop('If data outcome supplied, data and holdout must contain identical column names.')
    }
  }

  if (!is.numeric(C) | C < 0 | is.infinite(C)) {
    stop('C must be a finite, nonnegative scalar.')
  }

  if (!is.character(treated_column_name)) {
    stop('treated_column_name must be a character.')
  }

  if (!(treated_column_name %in% data_cols)) {
    stop('treated_column_name must be the name of a column in data.')
  }

  if (!(treated_column_name %in% holdout_cols)) {
    stop('treated_column_name must be the name of a column in holdout.')
  }

  if (is.factor(dplyr::pull(data, !!rlang::enquo(treated_column_name)))) {
    stop('Treated variable in data must be numeric binary or logical.')
  }

  if (is.factor(dplyr::pull(holdout, !!rlang::enquo(treated_column_name)))) {
    stop('Treated variable in holdout must be numeric binary or logical.')
  }

  if (!is.character(outcome_column_name)) {
    stop('Outcome_column_name must be a character.')
  }

  if (outcome_in_data & !(outcome_column_name %in% data_cols)) {
    stop('outcome_column_name must be the name of a column in data.')
  }

  if (!(outcome_column_name %in% holdout_cols)) {
    stop('outcome_column_name must be the name of a column in holdout.')
  }

  if (outcome_in_data && is.factor(dplyr::pull(data, !!rlang::enquo(outcome_column_name)))) {
    stop('Outcome variable in data must be numeric binary or continuous.')
  }

  if (is.factor(dplyr::pull(holdout, !!rlang::enquo(outcome_column_name)))) {
    stop('Outcome variable in holdout must be numeric binary or continuous')
  }

  if (!(binning_method %in% c('sturges', 'scott', 'fd'))) {
    stop("binning_method must be one of: 'sturges', 'scott', or 'fd'")
  }

  if (!(PE_method %in% c('ridge', 'xgb'))) {
    stop("PE_method must be one of 'ridge' or 'xgb'.
         To supply your own model to fit, use user_PE_fit.")
  }

  if (!is.logical(replace)) {
    stop('replace must be a logical scalar')
  }

  if (!(verbose %in% c(0, 1, 2, 3))) {
    stop('Verbose must be one of: 0, 1, 2, 3')
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

  if (!is.numeric(early_stop_epsilon) | early_stop_iterations <= 0) {
    stop('early_stop_epsilon must be a positive scalar')
  }

  if (!is.numeric(early_stop_un_c_frac) |
      early_stop_un_c_frac < 0 |
      early_stop_un_c_frac > 1) {
    stop('early_stop_un_c_frac must be a fraction between 0 and 1 (inclusive)')
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

  if (impute_with_outcome & !outcome_in_data) {
    stop('Outcome not present in data; cannot request to use it to impute missingness.')
  }

  if (!is.logical(impute_with_outcome)) {
    stop('impute_with_outcome must be a logical scalar')
  }

  if (!is.logical(impute_with_treatment)) {
    stop('impute_with_outcome must be a logical scalar')
  }
}
