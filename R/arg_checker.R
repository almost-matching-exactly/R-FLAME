check_args <-
  function(data, holdout, C, treatment_column_name, outcome_column_name,
           PE_method, user_PE_fun, PE_fun_params,
           repeats, verbose, want_pe, want_bf,
           early_stop_iterations, epsilon,
           early_stop_un_c_frac, early_stop_un_t_frac,
           early_stop_pe, early_stop_bf,
           missing_data_replace, missing_holdout_replace,
           missing_data_imputations, missing_holdout_imputations) {

  stopifnot(is.data.frame(data))
  stopifnot(is.data.frame(holdout))

  data_cols <- colnames(data)
  holdout_cols <- colnames(holdout)

  if (!identical(sort(data_cols), sort(holdout_cols))) {
    stop('Data and holdout must contain identical column names')
  }

  stopifnot(is.numeric(C) & C >= 0)

  stopifnot(is.character(treatment_column_name) &
              treatment_column_name %in% data_cols &
              treatment_column_name %in% holdout_cols)

  stopifnot(is.character(outcome_column_name) &
              outcome_column_name %in% data_cols &
              outcome_column_name %in% holdout_cols)

  stopifnot(PE_method %in% c('elasticnet'))
  # Checks for user function / params

  stopifnot(is.logical(repeats))
  stopifnot(verbose %in% c(0, 1, 2, 3))

  stopifnot(is.logical(want_pe))
  stopifnot(is.logical(want_bf))

  ## Early stop parameters
  stopifnot(is.numeric(early_stop_iterations) & early_stop_iterations >= 0)

  stopifnot(is.numeric(epsilon) & epsilon > 0)

  stopifnot(is.numeric(early_stop_un_c_frac) &
              (early_stop_un_c_frac >= 0 & early_stop_un_c_frac <= 1))
  stopifnot(is.numeric(early_stop_un_t_frac) &
              (early_stop_un_t_frac >= 0 & early_stop_un_t_frac <= 1))

  stopifnot(is.numeric(early_stop_pe) & early_stop_pe >= 0)
  stopifnot(is.numeric(early_stop_bf) & early_stop_bf >= 0 & early_stop_bf <= 2)

  ## Missing data parameters
  stopifnot(is.numeric(missing_data_replace) & missing_data_replace %in% c(0, 1, 2, 3))
  stopifnot(is.numeric(missing_holdout_replace) & missing_holdout_replace %in% c(0, 1, 2))
  stopifnot(is.numeric(missing_data_imputations) & missing_data_imputations >= 0)
  stopifnot(is.numeric(missing_holdout_imputations) & missing_holdout_imputations >= 0)
}
