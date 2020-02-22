aggregate_table <- function(list) {
  tab = table(as.character(list))
  tab = unclass(tab)
  name = names(tab)
  list_val = as.character(list)
  return(as.vector(tab[sapply(list_val, function(x) which(name ==  x))]))
}

# update_matched_bit takes a dataframe, a set of covariates to match on,
# the treatment indicator column and the matched indicator column.
# it returns the array indicating whether each unit is matched (the first return value),
# and a list of indices for the matched units (the second return value)

update_matched_bit <- function(data, covs, n_levels) {
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
  c_u = aggregate_table(b_u)

  # Compute c_u+
  c_u_plus = aggregate_table(b_u_plus)

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
                              PE_method, alpha) {
  ### Compute Balancing Factor

  # Calculate number of units unmatched (available)
  num_control <- sum(data$treated == 0)
  num_treated <- sum(data$treated == 1)

  # Number of matched units
  if (repeats) {
    match_index <-
      update_matched_bit(data, covs[-cov_to_drop], n_levels[-cov_to_drop]) %>%
      extract2('match_index')
    units_matched <- which(match_index)
  }
  else {
    match_index <-
      update_matched_bit(dplyr::filter(data, !matched), covs[-cov_to_drop], n_levels[-cov_to_drop]) %>%
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

process_matches <- function(data, repeats, covs, n_levels, MGs, matched_on, matching_covs, CATE, cov_names) {
  column <- colnames(data)
  if (repeats) {
    c(match_index, index) %<-% update_matched_bit(data, covs, n_levels)
    units_matched <- which(match_index)
  }
  else {
    c(match_index, index) %<-% update_matched_bit(dplyr::filter(data, !matched), covs, n_levels)
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


#' Bit Vectors Implementation
#'
#' \code{FLAME_bit} applies FLAME matching algorithm based on bit vectors.
#' The required arguments include (1) data and (2) holdout. The default model
#' for Match Quality is set to Ridge regression with 0.1 regularization parameter.
#'
#' @param data input data
#' @param holdout holdout training data
#' @param C Match Quality C parameter (optional, default =
#'   0.1)
#' @param PE_function user defined function to compute predictive error
#'   (optional)
#' @param model Linear, Ridge, or Lasso (optional)

#' Sizes, conditional average treatment effects (CATEs)
#' of matches at each iteration, (3) match quality at each iteration, and (4) the original
#' data with additional column *matched*, indicating the number of covariates each unit is
#' matched on. If a unit is never matched, then *matched* will be 0.
#' @examples
#' data(toy_data)
#' FLAME_bit(data = toy_data, holdout = toy_data)
#' @import dplyr
#' @import gmp
#' @import glmnet
#' @importFrom rlang .data
#' @importFrom graphics boxplot
#' @importFrom stats rbinom rnorm runif setNames
#' @importFrom stats lm var
#' @export

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
           missing_holdout_imputations = 10, missing_data_imputations = 0) {
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
                    matched_on, matching_covs, CATE, cov_names)

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
                 C, PE_method, alpha)
    # (First, in unlikely case of ties) covariate yielding the highest match quality
    drop <-
      sapply(MQ, function(x) C * x$BF - x$PE) %>%
      which.max()

    store_pe %<>% c(MQ[[drop]]$PE)
    store_bf %<>% c(MQ[[drop]]$BF)

    # Update covariates to match on
    covs <- covs[-drop]
    n_levels <- n_levels[-drop]

    # Make new matches having dropped a covariate
    ## Ideally should just return this from MQ so you don't have to redo it
    c(CATE, MGs, matched_on, units_matched, made_matches) %<-%
      process_matches(data, repeats, covs, n_levels, MGs, matched_on, matching_covs, CATE, cov_names)
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

