order_cov_names <- function(subset, cov_names, sorting_order) {
  return(subset[order(match(subset, cov_names[order(sorting_order)]))])
}

sort_cols <- function(df, treated_column_name, outcome_column_name,
                      type, is_missing) {

  n_covs <- ncol(df[[1]]) - 2 # Ignore treatment, outcome
  n_df <- length(df) # Always pass in a list of data frames

  # Treatment and outcome will be constant across imputations
  treatment_col <-
    df[[1]] %>%
    dplyr::select(!!rlang::enquo(treated_column_name))

  outcome_col <-
    df[[1]] %>%
    dplyr::select(!!rlang::enquo(outcome_column_name))

  n <- nrow(df[[1]])

  # For all imputed data sets
  for (i in 1:n_df) {
    tmp_df <- df[[i]]

    tmp_df %<>%
      dplyr::select(-c(!!rlang::enquo(treated_column_name),
                       !!rlang::enquo(outcome_column_name))) %>%
      cbind(outcome_col) %>%
      cbind(treatment_col)

    tmp_df[, 1:n_covs] %<>%
      lapply(as.factor)

    # Number of levels of each covariate
    n_levels <- sapply(tmp_df[, 1:n_covs, drop = FALSE], nlevels)

    # To sort covariates in increasing order of number of levels
    sorting_order <- order(n_levels)
    # To make sure the column names are also reordered
    cov_names <- colnames(tmp_df)[1:n_covs][sorting_order]

    # Sorted number of levels of each covariate
    n_levels <- n_levels[sorting_order]

    # Data sorted by n_levels
    tmp_df[, 1:n_covs] <- tmp_df[, sorting_order]
    # Sorting data column names
    colnames(tmp_df) <- c(cov_names, 'outcome', 'treated')

    if (type == 'data') {
      for (j in 1:n_covs) {
        levels(tmp_df[, j]) %<>% c('*')
      }
    }

    # covs denotes the covariates currently being matched on
    covs <- 1:n_covs

    # Denote whether a unit is matched and to how many others, respectively
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
              n_levels = n_levels,
              cov_names = cov_names,
              sorting_order = sorting_order))
}
