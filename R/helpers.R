order_cov_names <- function(subset, cov_names, sorting_order) {
  return(subset[order(match(subset, cov_names[order(sorting_order)]))])
}

organize_data <- function(data, holdout,
                          treatment_column_name, outcome_column_name) {

  # Rearrange data
  treatment_col <-
    data %>%
    dplyr::select(!!rlang::enquo(treatment_column_name))

  outcome_col <-
    data %>%
    dplyr::select(!!rlang::enquo(outcome_column_name))

  data %<>%
    dplyr::select(-c(!!rlang::enquo(treatment_column_name),
                     !!rlang::enquo(outcome_column_name))) %>%
    cbind(outcome_col) %>%
    cbind(treatment_col)

  ##
  treatment_col <-
    holdout %>%
    dplyr::select(!!rlang::enquo(treatment_column_name))

  outcome_col <-
    holdout %>%
    dplyr::select(!!rlang::enquo(outcome_column_name))

  holdout %<>%
    dplyr::select(-c(!!rlang::enquo(treatment_column_name),
                     !!rlang::enquo(outcome_column_name))) %>%
    cbind(outcome_col) %>%
    cbind(treatment_col)
  ##

  n_covs <- ncol(data) - 2 # ignore treatment and outcome

  data[, 1:n_covs] %<>%
    lapply(as.factor)
  holdout[, 1:n_covs] %<>%
    lapply(as.factor)

  # Number of levels of each covariate
  n_levels <- sapply(data[, 1:n_covs, drop = FALSE], nlevels)

  # Change levels to allow for 'unmatched on this covariate' indicator: '*'
  for (i in 1:n_covs) {
    levels(data[, i]) %<>% c('*')
  }

  # To sort covariates in increasing order of number of levels
  sorting_order <- order(n_levels)
  # To make sure the column names are also reordered
  cov_names <- colnames(data)[1:n_covs][sorting_order]

  # Sorted number of levels of each covariate
  n_levels <- n_levels[sorting_order]
  # Data sorted by n_levels
  data[, 1:n_covs] <- data[, sorting_order]
  # Sorting data column names
  colnames(data) <- c(cov_names, 'outcome', 'treated')

  # Holdout and its column names sorted by n_levels
  holdout[, 1:n_covs] <- holdout[, sorting_order]
  colnames(holdout) <- c(cov_names, 'outcome', 'treated')

  n <- nrow(data)
  # holdout_inds <- sample(1:n, size = round(holdout_data * n))
  # holdout_data <- data[holdout_inds, ]
  # data <- data[-holdout_inds, ]

  # # List of covariates used to match at each level
  # matching_covs <- list()

  # covs denotes the covariates currently being matched on
  covs <- 1:n_covs

  # Denotes the number of covariates each unit is matched on
  data$matched <- rep(FALSE, n)

  return(list(data = data,
              holdout = holdout,
              covs = covs,
              n_covs = n_covs,
              n_levels = n_levels,
              cov_names = cov_names,
              sorting_order = sorting_order))
}
