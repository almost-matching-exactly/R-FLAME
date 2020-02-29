aggregate_table <- function(vec, opt) {
  tab = table(as.character(vec))
  tab = unclass(tab)
  name = names(tab)
  list_val = as.character(vec)
  # browser()
  if (opt == 1) {
    return(as.vector(tab[match(list_val, name)]))
  }
  message('got here. tf?')
  return(as.vector(tab[sapply(list_val, function(x) which(name ==  x))]))
}

# update_matched_bit takes a dataframe, a set of covariates to match on,
# the treatment indicator column and the matched indicator column.
# it returns the array indicating whether each unit is matched (the first return
# value), and a list of indices for the matched units (the second return value)

update_matched_bit <- function(data, covs, n_levels, opt) {
  data_wo_t <- gmp::as.bigz(as.matrix(data[, covs]))
######## Do a massive dataset to check this

  # Compute b_u
  multiplier <- gmp::pow.bigz(n_levels, seq_along(n_levels) - 1)

  b_u <-
    gmp::`%*%`(data_wo_t, multiplier) %>%
    as.vector()

  stopifnot(gmp::is.bigz(b_u))

  # Compute b_u+
  multiplier <- gmp::pow.bigz(n_levels, seq_along(n_levels))

  b_u_plus <-
    gmp::`%*%`(data_wo_t, multiplier) %>%
    gmp::add.bigz(data$treated) %>%
    as.vector()
  stopifnot(gmp::is.bigz(b_u_plus))

  # Compute c_u
  c_u = aggregate_table(b_u, opt = opt)

  # Compute c_u+
  c_u_plus = aggregate_table(b_u_plus, opt = opt)

  match_index = mapply(function(x,y) (x != y) && (x >= 2) && (y >= 1), c_u, c_u_plus)
  index = b_u[match_index]

  return(list(match_index = match_index,
              index = index))
}

#match_quality function takes holdout dataset, number of total covariates,
#list of current covariates, covariate c to temporily remove from, and trafeoff
#parameter as input. The function then computes Balancing Factor and Predictive Error,
#returning Match Quality.

get_match_quality <- function(cov_to_drop, data, repeats, holdout, covs, n_levels, C,
                              PE_method, alpha, opt) {
  ### Compute Balancing Factor

  # Calculate number of units unmatched (available)
  num_control <- sum(data$treated == 0)
  num_treated <- sum(data$treated == 1)

  # Number of matched units
  if (repeats) {
    match_index <-
      update_matched_bit(data, setdiff(covs, covs_to_drop), n_levels[-which(covs == cov_to_drop)], opt = opt) %>%
      magrittr::extract2('match_index')

    # match_index <-
    #   update_matched_bit(data, covs[-cov_to_drop], n_levels[-cov_to_drop], opt = opt) %>%
    #   magrittr::extract2('match_index')
    units_matched <- which(match_index)
  }
  else {
    match_index <-
      update_matched_bit(dplyr::filter(data, !matched), setdiff(covs, cov_to_drop),
                         n_levels[-which(covs == cov_to_drop)], opt = opt) %>%
      magrittr::extract2('match_index')
    units_matched <- which(!data$matched)[match_index]
  }

  num_control_matched <- sum(data$treated[units_matched] == 0)
  num_treated_matched <- sum(data$treated[units_matched] == 1)

  BF <- dplyr::if_else(num_control == 0 | num_treated == 0, # Is this if_else really necessary? We should catch it earlier
                0,
                num_control_matched / num_control + num_treated_matched / num_treated)

  ### Compute Predictive Error

  # Default PE - elastic net regression with 0.1 regularization parameter
  if (PE_method == 'elasticnet') {
    PE <- predict_elasticnet(holdout, cov_to_drop, alpha)
  }
  else {
    stop("I don't recognize this prediction method")
  }
  return(list(BF = BF,
              PE = PE))
}

make_MGs <- function(data, index, matched_units, covs, cov_names) {
  # Takes all the units that were matched on these p' covariates and separates them
  # into matched groups based off their unique values of those covariates
  # Returns a list with three items:
  ## MGs: a list, each entry of which corresponds to a different MG and contains the indices of the corresponding units
  ## CATEs: a vector the same length as MGs, each entry of which is the CATE for the corresponding MG
  ## matched_on: a list the same length as MGs, each entry of which is a named vector detailing the covariates
  ##  and their values that the units in the corresponding MG matched on
  unique_MGs <- unique(index)
  n_MGs <- length(unique_MGs)

  MGs <- vector('list', length = n_MGs)
  CATEs <- vector('numeric', length = n_MGs)
  matched_on <- vector('list', length = n_MGs)

  for (i in 1:n_MGs) {
    members <- matched_units[which(index == unique_MGs[i])]
    MGs[[i]] <- members
    treated <- intersect(members, which(data$treated == 1))
    control <- intersect(members, which(data$treated == 0))
    CATEs[i] <- mean(data$outcome[treated]) - mean(data$outcome[control])
    matched_on[[i]] <-
      data[members[1], covs, drop = FALSE] %>%
      `rownames<-`(NULL)
    names(matched_on[[i]]) <- cov_names[covs]
  }
  return(list(MGs = MGs,
              CATEs = CATEs,
              matched_on = matched_on))
}

process_matches <- function(data, repeats, covs, n_levels, MGs, matched_on, matching_covs, CATE, cov_names, opt) {
  if (repeats) {
    match_out <- update_matched_bit(data, covs, n_levels, opt = opt)
    match_index <- match_out[[1]]
    index <- match_out[[2]]
    # c(match_index, index) %<-% update_matched_bit(data, covs, n_levels, opt = opt)
    units_matched <- which(match_index)
  }
  else {
    # c(match_index, index) %<-% update_matched_bit(dplyr::filter(data, !matched), covs, n_levels, opt = opt)
    match_out <- update_matched_bit(dplyr::filter(data, !matched), covs, n_levels, opt = opt)
    match_index <- match_out[[1]]
    index <- match_out[[2]]
    # adjusts for fact that match_index returned by update_matched bit is
    # with respect to the unmatched subset of data
    units_matched <- which(!data$matched)[match_index]
  }

  made_matches <- sum(match_index) > 0

  if (made_matches) {
    new_MGs <- make_MGs(data, index, units_matched, covs, cov_names)
    MGs <- c(MGs, new_MGs$MGs)
    CATE <- c(CATE, new_MGs$CATEs)
    matched_on <- c(matched_on, new_MGs$matched_on)
  }
  matching_covs[[length(matching_covs) + 1]] <- cov_names[covs]

  return(list(CATE = CATE,
              MGs = MGs,
              matched_on = matched_on,
              units_matched = units_matched,
              made_matches = made_matches))
}

get_PE <- function(cov_to_drop, covs, holdout, PE_method, user_PE_fun, PE_fun_params) {
  if (!is.null(user_PE_fun)) {
    PE_func <- user_PE_fun
  }
  else {
    if (PE_method == 'elasticnet') {
      PE_func <- glmnet::glmnet
      if (length(unique(holdout$outcome)) == 2) {
        family <- 'binomial'
        lambda <- 1
      }
      else {
        family <- 'gaussian'
        lambda <- 0.1
      }
      PE_fun_params <- list(alpha = 0, lambda = lambda, family = family)
    }
    else if (PE_method == 'xgb') {
      PE_func <- xgboost::xgboost
      PE_fun_params <- list(nrounds = 100, verbose = 0)
    }
    else {
      stop('PE_method not recognized. To supply your own function, use user_PE_fun')
    }
  }

  PE <- predict_master(holdout, covs, cov_to_drop, PE_func, PE_fun_params)
  return(PE)
}

get_BF <- function(cov_to_drop, data, repeats, covs, n_levels, opt) {
  # Calculate number of units unmatched (available)
  num_control <- sum(data$treated == 0)
  num_treated <- sum(data$treated == 1)

  # Number of matched units
  if (repeats) {
    match_index <-
      update_matched_bit(data, setdiff(covs, cov_to_drop), n_levels[-which(covs == cov_to_drop)], opt) %>%
      magrittr::extract2('match_index')
    # match_index <-
    #   update_matched_bit(data, covs[-cov_to_drop], n_levels[-cov_to_drop], opt) %>%
    #   magrittr::extract2('match_index')
    units_matched <- which(match_index)
  }
  else {
    match_index <-
      update_matched_bit(dplyr::filter(data, !matched), setdiff(covs, cov_to_drop),
                         n_levels[-which(covs == cov_to_drop)], opt) %>%
      magrittr::extract2('match_index')
    # match_index <-
    #   update_matched_bit(dplyr::filter(data, !matched), covs[-cov_to_drop], n_levels[-cov_to_drop], opt) %>%
    #   magrittr::extract2('match_index')
    units_matched <- which(!data$matched)[match_index]
  }

  # Newly matched
  num_control_matched <- sum(data$treated[units_matched] == 0)
  num_treated_matched <- sum(data$treated[units_matched] == 1)

  # All matched units; for stopping rule purposes
  all_unmatched <-
    setdiff(1:nrow(data), union(units_matched, which(data$matched)))

  n_control_unmatched <- sum(all_unmatched %in% which(data$treated == 0))
  n_treated_unmatched <- sum(all_unmatched %in% which(data$treated == 1))

  prop_control_unmatched <- n_control_unmatched / num_control
  prop_treated_unmatched <- n_treated_unmatched / num_treated

  # Is this if_else really necessary? We should catch it earlier
  BF <- dplyr::if_else(num_control == 0 | num_treated == 0,
                0,
                num_control_matched / num_control + num_treated_matched / num_treated)

  return(list(BF = BF,
              prop_unmatched =
                list(control = prop_control_unmatched,
                     treated = prop_treated_unmatched)))
}

#' Bit Vectors Implementation
#'
#' \code{FLAME_bit} runs the FLAME matching algorithm implemented using bit
#' vectors. The user is required to pass in \code{data}.
#' This is a longer description of what FLAME does.
#'
#' @param data Data to be matched. Either a dataframe or a path to a .csv file
#'   to be read into a dataframe.
#' @param holdout Data to be used to compute predictive error. If a numeric
#'   scalar between 0 and 1 (default = 0.1), that proportion of \code{data} will
#'   made into a holdout set to compute predictive error and only the remaining
#'   proportion of \code{data} will be matched. Otherwise, a dataframe or a path
#'   to a csv file.
#' @param treatment_column_name A character with the name of the treatment
#'   column in \code{data}.
#' @param outcome_column_name A character with the name of the outcome column
#'   in \code{data}.
#' @param PE_method One of "elasticnet", .... Denotes the method to be used for
#'   computation of PE.
#' @param alpha A positive scalar to be passed to glmnet for computation of PE
#'   via elastic nets.
#' @param C A positive scalar denoting the tradeoff between BF and PE. Higher C
#'   will prioritize more matches and lower C will prioritize not dropping
#'   important covariates.
#' @param repeats A logical scalar. If true, allows the same unit to be matched
#'   multiple times, on different numbers of covariates.
#' @param verbose Controls output while FLAME is running. If 0, no output. If 1,
#'   outputs the iteration every iteration. If 2, outputs the iteration and
#'   number of unmatched units every 5 iterations. If 3, outputs the iteration
#'   and number of unmatched units every 5 iterations.
#' @param want_pe A logical scalar. If TRUE, the predictive error (PE) at each
#'   iteration will be returned.
#' @param want_bf A logical scalar. If TRUE, the balancing factor (BF) at each
#'   iteration will be returned.
#' @section Parameters for Early Stopping:
#' @param early_stop_iterations A nonnegative integer, denoting the number of
#'   iterations of FLAME to be performed. If 0, one round of exact matching is
#'   performed before terminating.
#' @param stop_unmatched_c A logical scalar. If TRUE, FLAME will terminate when
#'   there are no more control units to be matched.
#' @param early_stop_un_c_frac A numeric value between 0 and 1 (inclusive). If
#'   the proportion of control units that are unmatched falls below this value,
#'   FLAME terminates.
#' @param stop_unmatched_t A logical scalar. If TRUE, FLAME will terminate when
#'   there are no more treated units to be matched.
#' @param early_stop_un_t_frac A numeric value between 0 and 1 (inclusive). If
#'   the proportion of treatment units that are unmatched falls below this
#'   value, FLAME terminates.
#' @param early_stop_pe A numeric value between 0 and 1 (inclusive).
#' @param early_stop_bf A logical scalar. If TRUE....
#' @param early_stop_bf A numeric value between 0 and 1 (inclusive).
#' @section Parameters for Missing Data:
#' @param missing_data_replace If 0, assumes no missingness in \code{data}. If
#'   1, eliminates units with missingness from \code{data}. If 2, will not match
#'   a unit on a covariate that it is missing. If 3, performs
#'   \code{missing_data_imputations} of MICE to impute the missing data.
#' @param missing_holdout_replace If 0, assumes no missing data in
#'   \code{holdout}. If 1, eliminates units with missingness from
#'   \code{holdout}. If 2,  performs \code{missing_holdout_imputations} of MICE
#'   to impute the missing data.
#' @param missing_holdout_imputations If \code{missing_holdout_replace} = 2,
#'   performs this many imputations of the missing data in \code{holdout} using
#'   MICE.
#' @param missing_data_imputations If \code{missing_data_replace} = 2, performs
#'   this many imputations of the missing data in \code{data} using MICE.

#' @examples
#' data <- gen_data()
#' holdout <- gen_data()
#' FLAME_out <- FLAME_bit(data = data, holdout = holdout)
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @importFrom zeallot %<-%
#' @importFrom rlang !!
#' @export
FLAME_bit <- function(data,
           treatment_column_name = 'treated',
           outcome_column_name='outcome',
           PE_method = 'elasticnet',
           alpha = 0.1, C = 0.1, holdout = FALSE,
           repeats = FALSE, verbose = 2, want_pe = TRUE, early_stop_iterations = Inf,
           stop_unmatched_c = FALSE, early_stop_un_c_frac = 0.1,
           stop_unmatched_t = FALSE, early_stop_un_t_frac = 0.1,
           early_stop_pe = 0.01,
           want_bf = FALSE, early_stop_bf = 0.01,
           missing_data_replace = 0, missing_holdout_replace = 0,
           missing_holdout_imputations = 10, missing_data_imputations = 0, opt = 0) {

  read_data_out <- read_data(data, holdout)
  data <- read_data_out[[1]]
  holdout <- read_data_out[[2]]
  # c(data, holdout) %<-% read_data(data, holdout)

  check_args(data, treatment_column_name, outcome_column_name,
             C, holdout,
             repeats, verbose, want_pe, early_stop_iterations,
             stop_unmatched_c, early_stop_un_c_frac,
             stop_unmatched_t, early_stop_un_t_frac,
             early_stop_pe, early_stop_pe,
             want_bf, early_stop_bf, early_stop_bf,
             missing_data_replace, missing_holdout_replace,
             missing_holdout_imputations, missing_data_imputations)

  # c(data, holdout, covs, n_covs, n_levels,
  #   cov_names, sorting_order) %<-%
  #   organize_data(data, holdout, treatment_column_name, outcome_column_name)

  organized_data <- organize_data(data, holdout, treatment_column_name, outcome_column_name)
  data <- organized_data[[1]]
  holdout <- organized_data[[2]]
  covs <- organized_data[[3]]
  n_covs <- organized_data[[4]]
  n_levels <- organized_data[[5]]
  cov_names <- organized_data[[6]]
  sorting_order <- organized_data[[7]]

  missing_data_out <- handle_missing_data(data, holdout,
                                          missing_data_replace, missing_holdout_replace,
                                          missing_data_imputations, missing_holdout_imputations)
  data <- missing_data_out[[1]]
  holdout <- missing_data_out[[2]]
  # c(data, holdout) %<-%
  #   handle_missing_data(data, holdout,
  #                       missing_data_replace, missing_holdout_replace,
  #                       missing_data_imputations, missing_holdout_imputations)

  # List of MGs, each entry contains the corresponding MGs entries
  MGs <- list()
  # List of CATEs, each entry contains the corresponding MGs CATE
  CATE <- vector('numeric')
  # List of covariates and their values matched on for each corresponding MG
  matched_on <- list()

  # List of covariates used to match at each level
  matching_covs <- list()

  # Try and make matches on all covariates
  c(CATE, MGs, matched_on, units_matched, made_matches) %<-%
    process_matches(data, repeats, covs, n_levels, MGs,
                    matched_on, matching_covs, CATE, cov_names, opt = opt)

  if (made_matches) {
    data$matched[units_matched] <- TRUE
  }
  store_pe <- NULL
  store_bf <- NULL

  iter <- 0
  while (length(covs) > 1 & !(all(data$matched))) {
    iter <- iter + 1

    # Compute the match quality associated with dropping each covariate
    MQ <- lapply(covs, get_match_quality, data, repeats, holdout, covs, n_levels,
                 C, PE_method, alpha, opt = opt)
    # PE <- sapply(covs, get_PE, holdout, alpha, PE_method)
    #
    # best_lower_bound <- max(-PE)
    # upper_bound <- 2 * C - PE
    # drop_candidates <- which(upper_bound > best_lower_bound)
    # PE <- PE[drop_candidates]
    #
    # BF <- sapply(drop_candidates, get_BF, data, repeats, covs, n_levels)
    # browser()
    # (First, in unlikely case of ties) covariate yielding the highest match quality
    drop <-
      sapply(MQ, function(x) C * x$BF - x$PE) %>%
      which.max()
    # MQ <- C * BF - PE
    # drop <- which.max(MQ)

    store_pe %<>% c(MQ[[drop]]$PE)
    store_bf %<>% c(MQ[[drop]]$BF)
    # store_pe %<>% c(PE[drop])
    # store_bf %<>% c(BF[drop])

    # Update covariates to match on
    covs <- covs[-drop]
    n_levels <- n_levels[-drop]
    # covs <- covs[-drop_candidates[drop]]
    # n_levels <- n_levels[-drop_candidates[drop]]

    # Make new matches having dropped a covariate
    ## Ideally should just return this from MQ so you don't have to redo it
    c(CATE, MGs, matched_on, units_matched, made_matches) %<-%
      process_matches(data, repeats, covs, n_levels, MGs, matched_on, matching_covs, CATE, cov_names, opt = opt)
    if (made_matches) {
      data[units_matched, setdiff(1:n_covs, covs)] <- '*' ## Same as covs_dropped?
      data$matched[units_matched] <- TRUE
    }
    show_progress(verbose, iter, data)
    if (early_stop(iter, data, early_stop_iterations,
                   store_pe, store_bf,
                   stop_unmatched_c, early_stop_un_c_frac,
                   stop_unmatched_t, early_stop_un_t_frac,
                   early_stop_bf, early_stop_bf,
                   early_stop_pe, early_stop_pe)) {
      break
    }
  }

  # Done matching!

  # Reorder the data according to the original column order
  data[, 1:n_covs] %<>% dplyr::select(order(sorting_order))
  # colnames(data) <- c(original_colnames, 'matched')

  ret_list <- list(MGs = MGs,
                   CATE = CATE,
                   matched_on = matched_on,
                   data = data,
                   matching_covs = matching_covs)

  if (want_pe) {
    ret_list %<>% c('PE' = list(store_pe))
  }
  if (want_bf) {
    ret_list %<>% c('BF' = list(store_bf))
  }

  return(ret_list)
}

# Check covariate order and also return_df

#' @export
FLAME_bit_new <-
  function(data, holdout = 0.1, C = 0.1, ## Algorithmic arguments
           treatment_column_name = 'treated', outcome_column_name = 'outcome',
           PE_method = 'elasticnet', user_PE_fun = NULL, PE_fun_params = NULL,
           repeats = FALSE, verbose = 2, want_pe = TRUE, want_bf = FALSE,
           early_stop_iterations = Inf, epsilon = 0.05, ## Early stopping arguments
           early_stop_un_c_frac = 0, early_stop_un_t_frac = 0,
           early_stop_pe = Inf, early_stop_bf = 0,
           missing_data_replace = 0, missing_holdout_replace = 0, ## Missing data arguments
           missing_data_imputations = 10, missing_holdout_imputations = 10, opt = 0) {

  read_data_out <- read_data(data, holdout)
  data <- read_data_out[[1]]
  holdout <- read_data_out[[2]]

  check_args(data, holdout, C,
             treatment_column_name, outcome_column_name,
             PE_method, user_PE_fun, PE_fun_params,
             repeats, verbose, want_pe, want_bf,
             early_stop_iterations, epsilon,
             early_stop_un_c_frac, early_stop_un_t_frac,
             early_stop_pe, early_stop_bf,
             missing_data_replace, missing_holdout_replace,
             missing_data_imputations, missing_holdout_imputations)

  missing_out <-
    handle_missing_data(data, holdout,
                        missing_data_replace, missing_holdout_replace,
                        missing_data_imputations, missing_holdout_imputations)

  data <- missing_out[[1]]
  holdout <- missing_out[[2]]

  c(data, covs, n_covs, n_levels, cov_names, sorting_order) %<-%
    sort_cols(data, treatment_column_name, outcome_column_name, type = 'data')

  holdout <-
    sort_cols(holdout, treatment_column_name, outcome_column_name,
              type = 'holdout') %>%
    magrittr::extract2(1)

  n_iters <- length(data)

  FLAME_out <- vector(mode = 'list', length = n_iters)
  for (i in 1:n_iters) {
    if (missing_data_replace == 3) {
      message('Running FLAME on imputed dataset ', i, ' of ', n_iters)
    }
    FLAME_out[[i]] <-
      FLAME_internal(data[[i]], holdout, covs, n_covs, n_levels, cov_names, sorting_order, C, PE_method, user_PE_fun, PE_fun_params,
                     repeats, verbose, want_pe, want_bf, early_stop_iterations,
                     epsilon = epsilon, early_stop_un_c_frac, early_stop_un_t_frac,
                     early_stop_pe, early_stop_bf, opt = opt)
  }

  if (n_iters == 1) {
    return(FLAME_out[[1]])
  }
  return(FLAME_out)
}

FLAME_internal <- function(data, holdout = 0.1, covs, n_covs, n_levels, cov_names, sorting_order, C = 0.1, ## Algorithmic arguments
                           PE_method = 'elasticnet', user_PE_fun = NULL, PE_fun_params = NULL,
                           repeats = FALSE, verbose = 2, want_pe = TRUE, want_bf = FALSE,
                           early_stop_iterations = Inf, epsilon = 0.05, ## Early stopping arguments
                           early_stop_un_c_frac = 0, early_stop_un_t_frac = 0,
                           early_stop_pe = Inf, early_stop_bf = 0, opt = 0) {
  # List of MGs, each entry contains the corresponding MGs entries
  MGs <- list()
  # List of CATEs, each entry contains the corresponding MGs CATE
  CATE <- vector('numeric')
  # List of covariates and their values matched on for each corresponding MG
  matched_on <- list()

  # List of covariates used to match at each level
  matching_covs <- list()
  covs_dropped <- NULL

  # Try and make matches on all covariates
  # c(CATE, MGs, matched_on, units_matched, made_matches) %<-%
    # process_matches(data, repeats, covs, n_levels, MGs,
    #                 matched_on, matching_covs, CATE, cov_names, opt = 1)
    #
  processed_matches <- process_matches(data, repeats, covs, n_levels, MGs,
                                       matched_on, matching_covs, CATE, cov_names, opt = 1)
  CATE <- processed_matches[[1]]
  MGs <- processed_matches[[2]]
  matched_on <- processed_matches[[3]]
  units_matched <- processed_matches[[4]]
  made_matches <- processed_matches[[5]]

  if (made_matches) {
    data$matched[units_matched] <- TRUE
    data$weights[units_matched] %<>% magrittr::add(1)
    matching_covs %<>% c(list(order_cov_names(cov_names[covs],
                                              cov_names,
                                              sorting_order)))
  }
  store_pe <- NULL
  store_bf <- NULL

  iter <- 0
  baseline_PE <- get_PE(cov_to_drop = NULL, covs, holdout,
                        PE_method, user_PE_fun, PE_fun_params)

  while (!early_stop(iter, data, covs, early_stop_iterations)) {
    iter <- iter + 1
    show_progress(verbose, iter, data)

    # Compute the PE associated with dropping each covariate
    PE <- sapply(covs, get_PE, covs, holdout, PE_method, user_PE_fun, PE_fun_params)
    if (early_stop_PE(min(PE), early_stop_pe, epsilon, baseline_PE)) {
      break
    }

    if (C != 0) {
      best_lower_bound <- max(-PE)
      upper_bound <- 2 * C - PE

      drop_candidates <- which(upper_bound >= best_lower_bound) # Should this be strictly greater than?
      PE <- PE[drop_candidates]

      BF_out <- lapply(covs[drop_candidates], get_BF, data, repeats, covs, n_levels, opt = opt)
      BF <- sapply(BF_out, function(x) x[['BF']])

      MQ <- C * BF - PE
    }
    else {
      MQ <- PE
    }

    # (First, in unlikely case of ties) covariate yielding the highest match quality
    drop <- which.max(MQ)
    prop_unmatched <- BF_out[[drop]][['prop_unmatched']]

    if (early_stop_BF(BF[drop], early_stop_bf,
                      prop_unmatched[['control']], prop_unmatched[['treated']],
                      early_stop_un_c_frac, early_stop_un_t_frac)) {
      break
    }

    store_pe %<>% c(PE[drop])
    store_bf %<>% c(BF[drop])

    covs_dropped <- c(covs_dropped, cov_names[covs[drop_candidates[drop]]])
    covs <- covs[-drop_candidates[drop]]
    n_levels <- n_levels[-drop_candidates[drop]]

    matching_covs %<>% c(list(order_cov_names(cov_names[covs],
                                              cov_names,
                                              sorting_order)))

    # Make new matches having dropped a covariate
    ## Ideally should just return this from MQ so you don't have to redo it
    # c(CATE, MGs, matched_on, units_matched, made_matches) %<-%
    #   process_matches(data, repeats, covs, n_levels,
    #                   MGs, matched_on, matching_covs, CATE, cov_names, opt = opt)
    processed_matches <- process_matches(data, repeats, covs, n_levels, MGs,
                                         matched_on, matching_covs, CATE, cov_names, opt = 1)
    CATE <- processed_matches[[1]]
    MGs <- processed_matches[[2]]
    matched_on <- processed_matches[[3]]
    units_matched <- processed_matches[[4]]
    made_matches <- processed_matches[[5]]

    if (made_matches) {
      data[units_matched, setdiff(1:n_covs, covs)] <- '*'
      data$matched[units_matched] <- TRUE
      data$weights[units_matched] %<>% magrittr::add(1)
    }
  }

  # Done matching!

  # Substitute covariate values of all unmatched units with the unmatched
  # covariate symbol '*'
  data[!data$matched, 1:n_covs] <- '*'

  # Reorder the data according to the original column order
  data[, 1:n_covs] %<>% dplyr::select(order(sorting_order))
  colnames(data) <-
    c(colnames(data)[1:n_covs][order(sorting_order)],
      'outcome', 'treated', 'matched', 'weights')

  ret_list <- list(MGs = MGs,
                   CATE = CATE,
                   matched_on = matched_on,
                   data = data,
                   matching_covs = matching_covs,
                   dropped = covs_dropped)

  if (want_pe) {
    ret_list %<>% c('PE' = list(store_pe))
  }
  if (want_bf) {
    ret_list %<>% c('BF' = list(store_bf))
  }

  return(ret_list)
}
