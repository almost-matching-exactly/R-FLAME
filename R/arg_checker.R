check_args <- function(data,
                      treatment_column_name,
                      outcome_column_name,
                      C, holdout,
                      repeats, verbose, want_pe, early_stop_iterations,
                      stop_unmatched_c, early_stop_un_c_frac,
                      stop_unmatched_t, early_stop_un_t_frac,
                      early_stop_pe, early_stop_pe_frac,
                      want_bf, early_stop_bf, early_stop_bf_frac,
                      missing_data_replace, missing_holdout_replace,
                      missing_holdout_imputations, missing_data_imputations) {

  stopifnot(is.data.frame(data))
  stopifnot(is.data.frame(holdout))

  data_cols <- colnames(data)
  holdout_cols <- colnames(holdout)

  if (!identical(sort(data_cols), sort(holdout_cols))) {
    stop('Data and holdout must contain identical column names')
  }
  # Ors or ands...?
  stopifnot(is.character(treatment_column_name) &
              treatment_column_name %in% data_cols &
              treatment_column_name %in% holdout_cols)

  stopifnot(is.character(outcome_column_name) &
              outcome_column_name %in% data_cols &
              outcome_column_name %in% holdout_cols)

  stopifnot(is.numeric(C) & C >= 0)
  stopifnot(is.logical(repeats))
  stopifnot(verbose %in% c(0, 1, 2, 3))

  stopifnot(is.logical(want_pe))
  stopifnot(is.logical(early_stop_pe))

  ## Early stop parameters
  stopifnot(is.numeric(early_stop_iterations) & early_stop_iterations >= 0)
  stopifnot(is.logical(stop_unmatched_c))
  stopifnot(is.numeric(early_stop_un_c_frac) &
              (early_stop_un_c_frac >= 0 & early_stop_un_c_frac <= 1))
  stopifnot(is.logical(stop_unmatched_t))
  stopifnot(is.numeric(early_stop_un_t_frac) &
              (early_stop_un_t_frac >= 0 & early_stop_un_t_frac <= 1))
  stopifnot(is.numeric(early_stop_pe_frac) & early_stop_pe_frac >= 0 & early_stop_pe_frac <= 1)
  stopifnot(is.logical(want_bf))
  stopifnot(is.logical(early_stop_bf))
  stopifnot(is.numeric(early_stop_bf_frac) & early_stop_bf_frac >= 0 & early_stop_bf_frac <= 1)

  ## Missing data parameters
  stopifnot(is.numeric(missing_data_replace) & missing_data_replace %in% c(0, 1, 2, 3))
  stopifnot(is.numeric(missing_holdout_replace) & missing_holdout_replace %in% c(0, 1, 2))
  stopifnot(is.numeric(missing_data_imputations) & missing_data_imputations >= 0)
  stopifnot(is.numeric(missing_holdout_imputations) & missing_holdout_imputations >= 0)
}
