check_args <-
  function(data, holdout, C,
           treated_column_name, outcome_column_name, n_flame_iters,
           PE_method, user_PE_fit, user_PE_fit_params,
           user_PE_predict, user_PE_predict_params,
           replace, verbose, want_pe, want_bf,
           early_stop_params,
           missing_holdout_imputations,
           impute_with_outcome, impute_with_treatment) {

  if (!is.data.frame(data)) {
    stop('`data` must be a data frame or a character denoting a .csv file ',
         'in the working directory.')
  }
  if (!is.data.frame(holdout)) {
    stop('`holdout` must be a data frame, a character denoting a .csv ' ,
         'file in the working directory, or a numeric proportion of ',
         '`data` to use as a holdout set.')
  }

  outcome_in_data <- !is.null(data[[outcome_column_name]])

  data_cols <- colnames(data)
  holdout_cols <- colnames(holdout)

  # Seems like I can remove but maybe there's a reason I put it up here...?
  if (!(outcome_column_name %in% holdout_cols)) {
    stop('`holdout` must contain outcome column with name `outcome_column_name`')
  }

  cov_inds_data <- which(!(colnames(data) %in%
                             c(treated_column_name, outcome_column_name)))
  cov_inds_holdout <- which(!(colnames(holdout) %in%
                                c(treated_column_name, outcome_column_name)))

  if (!identical(colnames(data)[cov_inds_data],
                 colnames(holdout[cov_inds_holdout]))) {
    stop('Covariate columns of `data` and `holdout` ',
         'must have have the same names and be in the same order.')
  }

  # if (!outcome_in_data) {
  #   if (!identical(data_cols,
  #                  holdout_cols[-match(outcome_column_name, holdout_cols)])) {
  #     stop(paste('Non-outcome columns of data and holdout',
  #                'must have identical names.'))
  #   }
  #
  #   data_cov_inds <-
  #     which(!(colnames(data) %in% c(treated_column_name)))
  #   holdout_cov_inds <-
  #     which(!(colnames(data) %in% c(treated_column_name, outcome_column_name)))
  #
  #   sapply(seq_along(data_cov_inds), function(x) {
  #     if (!identical(sort(levels(data[[data_cov_inds[x]]])),
  #                    sort(levels(holdout[[holdout_cov_inds[x]]])))) {
  #       stop('Levels of `data` and `holdout` must be identical')
  #     }
  #   })
  # }
  # else {
  #   if (!identical(data_cols, holdout_cols)) {
  #     stop(paste('If data outcome supplied, data and holdout must contain',
  #                'identical column names.'))
  #   }
  #   cov_inds <-
  #     which(!(colnames(data) %in% c(treated_column_name, outcome_column_name)))
  #   sapply(cov_inds, function(x) {
  #     if (!identical(sort(levels(data[[x]])), sort(levels(holdout[[x]])))) {
  #       stop('Levels of `data` and `holdout` must be identical')
  #     }
  #     })
  # }

  if (!is.numeric(n_flame_iters) | n_flame_iters < 0) {
    stop('`n_flame_iters` must be a nonnegative scalar')
  }

  if (!is.numeric(C) | C < 0 | is.infinite(C)) {
    stop('C must be a finite, nonnegative scalar.')
  }

  if (!is.character(treated_column_name)) {
    stop('`treated_column_name` must be a character.')
  }

  if (!(treated_column_name %in% data_cols)) {
    stop('`treated_column_name` must be the name of a column in `data.`')
  }

  if (!(treated_column_name %in% holdout_cols)) {
    stop('`treated_column_name` must be the name of a column in `holdout.`')
  }

  if (is.factor(data[[treated_column_name]])) {
    stop('Treated variable in `data` must be numeric binary or logical.')
  }

  if (is.factor(holdout[[treated_column_name]])) {
    stop('Treated variable in `holdout` must be numeric binary or logical.')
  }

  if (!is.character(outcome_column_name)) {
    stop('`outcome_column_name` must be a character.')
  }

  if (outcome_in_data & !(outcome_column_name %in% data_cols)) {
    stop('`outcome_column_name` must be the name of a column in `data.` ',
         'if outcome is supplied.')
  }

  if (!(outcome_column_name %in% holdout_cols)) {
    stop('`outcome_column_name` must be the name of a column in `holdout.`')
  }

  if (!(PE_method %in% c('ridge', 'xgb'))) {
    stop("`PE_method` must be one of 'ridge' or 'xgb'. ",
         "To supply your own model to fit, use user_PE_fit.")
  }

  if (!is.null(user_PE_fit_params) & is.null(user_PE_fit)) {
    stop("May not have `user_PE_fit` be NULL if ",
         "`user_PE_fit_params` is not NULL.")
  }

  if (!is.null(user_PE_predict_params) & is.null(user_PE_predict)) {
    stop("May not have `user_PE_predict` be NULL if ",
         "`user_PE_predict_params` is not NULL.")
  }

  if (!is.null(user_PE_fit_params)) {
    if (!is.list(user_PE_fit_params)) {
      stop("`user_PE_fit_params` must be a list")
    }
  }

  if (!is.null(user_PE_predict_params)) {
    if (!is.list(user_PE_predict_params)) {
      stop("`user_PE_predict_params` must be a list")
    }
  }

  if (!is.logical(replace)) {
    stop('`replace` must be a logical scalar')
  }

  if (!(verbose %in% c(0, 1, 2, 3))) {
    stop('`verbose` must be one of: 0, 1, 2, 3')
  }

  if (!is.logical(want_pe)) {
    stop('`want_pe` must be a logical scalar')
  }

  if (!is.logical(want_bf)) {
    stop('`want_bf` must be a logical scalar')
  }

  ## Early stop parameters
  if (!is.numeric(early_stop_params$iterations) | early_stop_params$iterations < 0) {
    stop('`early_stop_iterations` must be a nonnegative scalar')
  }

  if (!is.numeric(early_stop_params$epsilon) | early_stop_params$epsilon <= 0) {
    stop('`early_stop_epsilon` must be a positive scalar')
  }

  if (!is.numeric(early_stop_params$control) |
      early_stop_params$control < 0 |
      early_stop_params$control > 1) {
    stop('`early_stop_control` must be a fraction between 0 and 1 (inclusive)')
  }

  if (!is.numeric(early_stop_params$treated) |
      early_stop_params$treated < 0 |
      early_stop_params$treated > 1) {
    stop('`early_stop_treated` must be a fraction between 0 and 1 (inclusive).')
  }

  if (!is.numeric(early_stop_params$PE) | early_stop_params$PE < 0) {
    stop('`early_stop_pe` must be a nonnegative scalar')
  }

  if (!is.numeric(early_stop_params$BF) |
      early_stop_params$BF < 0 |
      early_stop_params$BF > 2) {
    stop('`early_stop_bf` must be a scalar between 0 and 2 (inclusive)')
  }

  ## Missing data parameters
  if (!is.numeric(missing_holdout_imputations) |
      missing_holdout_imputations < 1) {
    stop('`missing_holdout_imputations` must be an integer greater than 1')
  }

  if (impute_with_outcome & !outcome_in_data) {
    stop('Outcome not present in `data`; ',
         'cannot request to use it to impute missingness.')
  }

  if (!is.logical(impute_with_outcome)) {
    stop('`impute_with_outcome` must be a logical scalar')
  }

  if (!is.logical(impute_with_treatment)) {
    stop('`impute_with_outcome` must be a logical scalar')
  }
}
