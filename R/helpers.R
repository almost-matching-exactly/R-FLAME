order_cov_names <- function(subset, cov_names, sorting_order) {
  return(subset[order(match(subset, cov_names[order(sorting_order)]))])
}

show_progress <- function(verbose, iter, data) {
  if (verbose == 1) {
    message(paste('Starting iteration', iter + 1, 'of FLAME'))
  }
  else if (verbose == 2) {
    if ((iter + 1) %% 5 == 0) {
      message(paste('Starting iteration', iter + 1, 'of FLAME'))
      message(paste(sum(!data$matched), 'unmatched units remaining'))
    }
  }
  else if (verbose == 3) {
    message(paste('Starting iteration', iter + 1, 'of FLAME'))
    message(paste(sum(!data$matched), 'unmatched units remaining'))
  }
}

organize_data <- function(data, holdout,
                          treatment_column_name, outcome_column_name) {

  # Rearrange data
  treatment_col <-
    data %>%
    dplyr::select(!!enquo(treatment_column_name))

  outcome_col <-
    data %>%
    dplyr::select(!!enquo(outcome_column_name))

  data %<>%
    dplyr::select(-c(!!enquo(treatment_column_name),
                     !!enquo(outcome_column_name))) %>%
    cbind(outcome_col) %>%
    cbind(treatment_col)

  ##
  treatment_col <-
    holdout %>%
    dplyr::select(!!enquo(treatment_column_name))

  outcome_col <-
    holdout %>%
    dplyr::select(!!enquo(outcome_column_name))

  holdout %<>%
    dplyr::select(-c(!!enquo(treatment_column_name),
                     !!enquo(outcome_column_name))) %>%
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

early_stop <- function(verbose, iter, covs, data, early_stop_iterations,
                       store_pe, store_bf,
                       stop_unmatched_c, early_stop_un_c_frac,
                       stop_unmatched_t, early_stop_un_t_frac,
                       early_stop_bf, early_stop_bf_frac,
                       early_stop_pe, early_stop_pe_frac) {

  type <- early_stop_type(iter, covs, data, early_stop_iterations,
                          store_pe, store_bf,
                          stop_unmatched_c, early_stop_un_c_frac,
                          stop_unmatched_t, early_stop_un_t_frac,
                          early_stop_bf, early_stop_bf_frac,
                          early_stop_pe, early_stop_pe_frac)
  if (verbose == 0) {
    if (type == 0) {
      return(FALSE)
    }
    return(TRUE)
  }

  if (type == 0) {
    show_progress(verbose, iter, data)
    return(FALSE)
  }

  if (type == 1) {
    message('FLAME stopping: all covariates dropped')
  }
  if (type == 2) {
    message('FLAME stopping: all units matched')
  }
  if (type == 3) {
    message('FLAME stopping: completed ', iter, ' iterations')
  }
  if (type == 4) {
    message('FLAME stopping: balancing factor dropped below ', early_stop_bf_frac)
  }
  if (type == 5) {
    message('FLAME stopping: predictive error dropped below ', early_stop_pe_frac)
  }
  if (type == 6) {
    message('FLAME stopping: all control units have been matched')
  }
  if (type == 7) {
    message('FLAME stopping: all treatment units have been matched')
  }
  if (type == 8) {
    message('FLAME stopping: proportion of control units ',
    'that are unmatched dropped below ', early_stop_un_c_frac)
  }
  if (type == 9) {
    message('FLAME stopping: proportion of treatment units ',
    'that are unmatched dropped below ', early_stop_un_t_frac)
  }
  return(TRUE)
}

early_stop_type <-
  function(iter, covs, data, early_stop_iterations,
           store_pe, store_bf,
           stop_unmatched_c, early_stop_un_c_frac,
           stop_unmatched_t, early_stop_un_t_frac,
           early_stop_bf, early_stop_bf_frac,
           early_stop_pe, early_stop_pe_frac) {

  if (length(covs) == 1) {
    return(1)
  }

  if (all(data$matched)) {
    return(2)
  }

  if (iter >= early_stop_iterations) {
    return(3)
  }
  if (early_stop_bf) {
    if (store_bf[iter] < early_stop_bf_frac) {
      return(4)
    }
  }
  if (early_stop_pe) {
    if (store_pe[iter] < early_stop_pe_frac) { ###### PE is bad so should this be > than?
      return(5)
    }
  }
  if (stop_unmatched_c) {
    if (all(data$treated[!data$matched] == 0)) {
      return(6)
    }
  }
  if (stop_unmatched_t) {
    if (all(data$treated[!data$matched] == 1)) {
      return(7)
    }
  }
  if (sum(data$treated == 0 & !data$matched) / sum(data$treated == 0) <
      early_stop_un_c_frac) {
    return(8)
  }
  if (sum(data$treated == 1 & !data$matched) / sum(data$treated == 1) <
      early_stop_un_t_frac) {
    return(9)
  }
  return(0)
}
