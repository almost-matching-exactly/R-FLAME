aggregate_table <- function(vec, opt) {
  tab = table(as.character(vec))
  tab = unclass(tab)
  name = names(tab)
  list_val = as.character(vec)
  # browser()
  if (opt == 1) {
    return(as.vector(tab[match(list_val, name)]))
  }
  return(as.vector(tab[sapply(list_val, function(x) which(name ==  x))]))
}

# update_matched_bit takes a dataframe, a set of covariates to match on,
# the treatment indicator column and the matched indicator column.
# it returns the array indicating whether each unit is matched (the first return
# value), and a list of indices for the matched units (the second return value)

update_matched_bit <- function(data, covs, n_levels, opt) {
  data_wo_t <- as.bigz(as.matrix(data[, covs]))
######## Do a massive dataset to check this

  # Compute b_u
  multiplier <- pow.bigz(n_levels, seq_along(n_levels) - 1)

  b_u <- gmp::`%*%`(data_wo_t, multiplier) %>%
    as.vector()
  # browser()
  # Compute b_u+
  multiplier <- pow.bigz(n_levels, seq_along(n_levels))

  b_u_plus <-
    gmp::`%*%`(data_wo_t, multiplier) %>%
    add.bigz(data$treated) %>%
    as.vector()

  # browser()
  # browser()
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
      extract2('match_index')
    # match_index <-
    #   update_matched_bit(data, covs[-cov_to_drop], n_levels[-cov_to_drop], opt = opt) %>%
    #   extract2('match_index')
    units_matched <- which(match_index)
  }
  else {
    match_index <-
      update_matched_bit(dplyr::filter(data, !matched), setdiff(covs, cov_to_drop),
                         n_levels[-which(covs == cov_to_drop)], opt = opt) %>%
      extract2('match_index')
    units_matched <- which(!data$matched)[match_index]
  }

  num_control_matched <- sum(data$treated[units_matched] == 0)
  num_treated_matched <- sum(data$treated[units_matched] == 1)

  BF <- if_else(num_control == 0 | num_treated == 0, # Is this if_else really necessary? We should catch it earlier
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
  column <- colnames(data)
  if (repeats) {
    c(match_index, index) %<-% update_matched_bit(data, covs, n_levels, opt = opt)
    units_matched <- which(match_index)
  }
  else {
    c(match_index, index) %<-% update_matched_bit(dplyr::filter(data, !matched), covs, n_levels, opt = opt)
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

get_PE <- function(cov_to_drop, covs, holdout, alpha, PE_method) {
  if (PE_method == 'elasticnet') {
    PE <- predict_elasticnet(holdout, covs, cov_to_drop, alpha)
  }
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
      extract2('match_index')
    # match_index <-
    #   update_matched_bit(data, covs[-cov_to_drop], n_levels[-cov_to_drop], opt) %>%
    #   extract2('match_index')
    units_matched <- which(match_index)
  }
  else {
    match_index <-
      update_matched_bit(dplyr::filter(data, !matched), setdiff(covs, cov_to_drop),
                         n_levels[-which(covs == cov_to_drop)], opt) %>%
      extract2('match_index')
    # match_index <-
    #   update_matched_bit(dplyr::filter(data, !matched), covs[-cov_to_drop], n_levels[-cov_to_drop], opt) %>%
    #   extract2('match_index')
    units_matched <- which(!data$matched)[match_index]
  }

  num_control_matched <- sum(data$treated[units_matched] == 0)
  num_treated_matched <- sum(data$treated[units_matched] == 1)

  BF <- if_else(num_control == 0 | num_treated == 0, # Is this if_else really necessary? We should catch it earlier
                0,
                num_control_matched / num_control + num_treated_matched / num_treated)
  return(BF)
}

#' Bit Vectors Implementation
#'
#' \code{FLAME_bit} runs the FLAME matching algorithm implemented using bit
#' vectors. The user is required to pass in \code{data}.
#' This is a longer description of what FLAME does.
#'
#' @param data Data to be matched. Either a dataframe or a path to a .csv file
#'   to be read into a dataframe.
#' @param holdout Data to be used to compute predictive error. If FALSE
#'   (default), 10% of \code{data} will be used for this purpose and only the
#'   remaining 90% of \code{data} will be matched. Otherwise, a dataframe or a
#'   path to a csv file.
#' @param treatment_column_name A character with the name of the treatment
#'   column in \code{data}.
#' @param outcome_column_name A character with the name of the outcomee column
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
#' @param early_stop_pe A logical scalar. If TRUE....
#' @param early_stop_pe_frac A numeric value between 0 and 1 (inclusive).
#' @param early_stop_bf A logical scalar. If TRUE....
#' @param early_stop_bf_frac A numeric value between 0 and 1 (inclusive).
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
#' data <- gen_data(100, 5)
#' FLAME_bit(data = data,)

FLAME_bit <- function(data,
           treatment_column_name = 'treated',
           outcome_column_name='outcome',
           PE_method = 'elasticnet',
           alpha = 0.1, C = 0.1, holdout = FALSE,
           repeats = FALSE, verbose = 2, want_pe = TRUE, early_stop_iterations = Inf,
           stop_unmatched_c = FALSE, early_stop_un_c_frac = 0.1,
           stop_unmatched_t = FALSE, early_stop_un_t_frac = 0.1,
           early_stop_pe = FALSE, early_stop_pe_frac = 0.01,
           want_bf = FALSE, early_stop_bf = FALSE, early_stop_bf_frac = 0.01,
           missing_data_replace = 0, missing_holdout_replace = 0,
           missing_holdout_imputations = 10, missing_data_imputations = 0, opt = 0) {
  require(gmp)
  require(glmnet)
  require(dplyr)
  require(magrittr)
  require(zeallot)

  c(data, holdout) %<-% read_data(data, holdout)

  check_args(data, treatment_column_name, outcome_column_name,
             alpha, C, holdout,
             repeats, verbose, want_pe, early_stop_iterations,
             stop_unmatched_c, early_stop_un_c_frac,
             stop_unmatched_t, early_stop_un_t_frac,
             early_stop_pe, early_stop_pe_frac,
             want_bf, early_stop_bf, early_stop_bf_frac,
             missing_data_replace, missing_holdout_replace,
             missing_holdout_imputations, missing_data_imputations)

  c(data, holdout, covs, n_covs, n_levels,
    cov_names, original_colnames, sorting_order) %<-%
    organize_data(data, holdout, treatment_column_name, outcome_column_name)

  c(data, holdout) %<-%
    handle_missing_data(data, holdout,
                        missing_data_replace, missing_holdout_replace,
                        missing_data_imputations, missing_holdout_imputations)

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
      data[units_matched, setdiff(1:n_covs, covs)] <- '*' ## Same as covs?
      data$matched[units_matched] <- TRUE
    }
    show_progress(verbose, iter, data)
    if (early_stop(iter, data, early_stop_iterations,
                   store_pe, store_bf,
                   stop_unmatched_c, early_stop_un_c_frac,
                   stop_unmatched_t, early_stop_un_t_frac,
                   early_stop_bf, early_stop_bf_frac,
                   early_stop_pe, early_stop_pe_frac)) {
      break
    }
  }

  # Done matching!

  # Reorder the data according to the original column order
  data[, 1:n_covs] %<>% dplyr::select(order(sorting_order))
  colnames(data) <- c(original_colnames, 'matched')

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

FLAME_bit_new <- function(data,
                      treatment_column_name = 'treated',
                      outcome_column_name='outcome',
                      PE_method = 'elasticnet',
                      alpha = 0.1, C = 0.1, holdout = FALSE,
                      repeats = FALSE, verbose = 2, want_pe = TRUE, early_stop_iterations = Inf,
                      stop_unmatched_c = FALSE, early_stop_un_c_frac = 0.1,
                      stop_unmatched_t = FALSE, early_stop_un_t_frac = 0.1,
                      early_stop_pe = FALSE, early_stop_pe_frac = 0.01,
                      want_bf = FALSE, early_stop_bf = FALSE, early_stop_bf_frac = 0.01,
                      missing_data_replace = 0, missing_holdout_replace = 0,
                      missing_holdout_imputations = 10, missing_data_imputations = 0, opt = 0) {
  require(gmp)
  require(glmnet)
  require(dplyr)
  require(magrittr)
  require(zeallot)

  c(data, holdout) %<-% read_data(data, holdout)

  check_args(data, treatment_column_name, outcome_column_name,
             alpha, C, holdout,
             repeats, verbose, want_pe, early_stop_iterations,
             stop_unmatched_c, early_stop_un_c_frac,
             stop_unmatched_t, early_stop_un_t_frac,
             early_stop_pe, early_stop_pe_frac,
             want_bf, early_stop_bf, early_stop_bf_frac,
             missing_data_replace, missing_holdout_replace,
             missing_holdout_imputations, missing_data_imputations)

  c(data, holdout, covs, n_covs, n_levels,
    cov_names, original_colnames, sorting_order) %<-%
    organize_data(data, holdout, treatment_column_name, outcome_column_name)

  c(data, holdout) %<-%
    handle_missing_data(data, holdout,
                        missing_data_replace, missing_holdout_replace,
                        missing_data_imputations, missing_holdout_imputations)

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
    # MQ <- lapply(covs, get_match_quality, data, repeats, holdout, covs, n_levels,
    #              C, PE_method, alpha)
    PE <- sapply(covs, get_PE, covs, holdout, alpha, PE_method)

    best_lower_bound <- max(-PE)
    upper_bound <- 2 * C - PE
    # browser()
    drop_candidates <- which(upper_bound >= best_lower_bound) # Should this be strictly greater than?
    PE <- PE[drop_candidates]

    BF <- sapply(covs[drop_candidates], get_BF, data, repeats, covs, n_levels, opt = opt)
    # (First, in unlikely case of ties) covariate yielding the highest match quality
    # drop <-
    #   sapply(MQ, function(x) C * x$BF - x$PE) %>%
    #   which.max()
    MQ <- C * BF - PE
    drop <- which.max(MQ)

    # store_pe %<>% c(MQ[[drop]]$PE)
    # store_bf %<>% c(MQ[[drop]]$BF)
    store_pe %<>% c(PE[drop])
    store_bf %<>% c(BF[drop])
######## drop_candidates[drop] is the **entry** of covs that we drop
    # Update covariates to match on
    # covs <- covs[-drop]
    # n_levels <- n_levels[-drop]
    # browser()
    covs_dropped <- c(covs_dropped, cov_names[drop_candidates[drop]])
    covs <- covs[-drop_candidates[drop]]
    n_levels <- n_levels[-drop_candidates[drop]] # drop covs[dropcands[drop]]

    # Make new matches having dropped a covariate
    ## Ideally should just return this from MQ so you don't have to redo it
    c(CATE, MGs, matched_on, units_matched, made_matches) %<-%
      process_matches(data, repeats, covs, n_levels, MGs, matched_on, matching_covs, CATE, cov_names, opt = opt)
    if (made_matches) {
      data[units_matched, setdiff(1:n_covs, covs)] <- '*' ## Same as covs?
      data$matched[units_matched] <- TRUE
    }
    show_progress(verbose, iter, data)
    if (early_stop(iter, data, early_stop_iterations,
                   store_pe, store_bf,
                   stop_unmatched_c, early_stop_un_c_frac,
                   stop_unmatched_t, early_stop_un_t_frac,
                   early_stop_bf, early_stop_bf_frac,
                   early_stop_pe, early_stop_pe_frac)) {
      break
    }
  }

  # Done matching!

  # Reorder the data according to the original column order
  data[, 1:n_covs] %<>% dplyr::select(order(sorting_order))
  colnames(data) <- c(original_colnames, 'matched')

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

# Check covariate order and also return_df

