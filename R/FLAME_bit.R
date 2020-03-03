aggregate_table <- function(vec) {
  tab = table(as.character(vec))
  tab = unclass(tab)
  name = names(tab)
  list_val = as.character(vec)
  return(as.vector(tab[match(list_val, name)]))
}

# update_matched_bit takes a dataframe, a set of covariates to match on,
# the treatment indicator column and the matched indicator column.
# it returns the array indicating whether each unit is matched (the first return
# value), and a list of indices for the matched units (the second return value)

update_matched_bit <- function(data, covs, n_levels) {
  data_wo_t <- gmp::as.bigz(as.matrix(data[, covs]))

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
  c_u <- aggregate_table(b_u)

  # Compute c_u+
  c_u_plus <- aggregate_table(b_u_plus)

  match_index <-
    mapply(function(x,y) (x != y) && (x >= 2) && (y >= 1), c_u, c_u_plus)
  index <- b_u[match_index]

  return(list(match_index = match_index,
              index = index))
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

process_matches <-
  function(data, replace, covs, n_levels, MGs,
           matched_on, matching_covs, CATE, cov_names) {
  if (replace) {
    match_out <- update_matched_bit(data[!data$missing, ], covs, n_levels)
    match_index <- match_out[[1]]
    index <- match_out[[2]]
    units_matched <- which(!data$missing)[match_index]
  }
  else {
    match_out <-
      update_matched_bit(data[!data$matched & !data$missing, ], covs, n_levels)
    match_index <- match_out[[1]]
    index <- match_out[[2]]

    # Adjusts for fact that match_index returned by update_matched bit is
    # with respect to the unmatched subset of data
    units_matched <- which(!data$matched & !data$missing)[match_index]
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

get_PE <- function(cov_to_drop, covs, holdout, PE_method,
                   user_PE_fit, user_PE_fit_params,
                   user_PE_predict, user_PE_predict_params) {

  if (!is.null(user_PE_fit)) {
    PE_fit <- user_PE_fit
    PE_fit_params <- user_PE_fit_params
  }
  else {
    if (PE_method == 'elasticnet') {
      PE_fit <- glmnet::cv.glmnet
      if (length(unique(holdout$outcome)) == 2) {
        family <- 'binomial'
      }
      else {
        family <- 'gaussian'
      }
      PE_fit_params <- list(family = family, nfolds = 5)
    }
    else if (PE_method == 'xgb') {
      PE_fit <- cv_xgboost
      PE_fit_params <- list()
    }
    else {
      stop('PE_method not recognized.
           To supply your own function, use user_PE_fit and user_PE_predict')
    }
  }

  if (!is.null(user_PE_predict)) {
    PE_predict <- user_PE_predict
    PE_predict_params <- user_PE_predict_params
  }
  else {
    PE_predict <- predict
    PE_predict_params <- list()
  }

  PE <- predict_master(holdout, covs, cov_to_drop,
                       PE_fit, PE_predict, PE_fit_params, PE_predict_params)
  return(PE)
}

get_BF <- function(cov_to_drop, data, replace, covs, n_levels) {
  # Calculate number of units unmatched (available)
  num_control <- sum(data$treated == 0)
  num_treated <- sum(data$treated == 1)

  # Number of matched units
  if (replace) {
    match_index <-
      update_matched_bit(data[!data$missing, ], setdiff(covs, cov_to_drop),
                         n_levels[-which(covs == cov_to_drop)])[['match_index']]
    units_matched <- which(!data$missing)[match_index]
  }
  else {
    match_index <-
      update_matched_bit(data[!data$matched & !data$missing, ],
                         setdiff(covs, cov_to_drop),
                         n_levels[-which(covs == cov_to_drop)])[['match_index']]
    units_matched <- which(!data$matched & !data$missing)[match_index]
  }

  if (!replace & any(units_matched %in% which(data$matched))) {
    browser()
    stop('Rematching a matched unit')
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
                num_control_matched / num_control +
                  num_treated_matched / num_treated)

  return(list(BF = BF,
              prop_unmatched =
                list(control = prop_control_unmatched,
                     treated = prop_treated_unmatched)))
}

#' Bit Vectors Implementation of FLAME
#'
#' \code{FLAME} runs the bit-vectors implementation of the FLAME algorithm.
#'
#' @section Introduction:
#' FLAME is a matching algorithm for causal inference that matches units if they
#' match exactly on certain covariates. It starts by making any possible matches
#' on all covariates. It then drops a covariate, makes any possible matches on
#' the remaining covariates, and repeats this process until stopping. The
#' covariate dropped at any given iteration is the one yielding the greatest
#' match quality \eqn{MQ}, defined as \eqn{MQ = C \times BF - PE}. Here,
#' \eqn{BF} denotes the balancing factor, defined as the proportion of control
#' units, plus the proportion of treated units, that can be newly matched by
#' dropping that covariate. And \eqn{PE} denotes the prediction error, defined
#' as the training error incurred when predicting the outcome from covariates on
#' a separate, holdout set. In this way, FLAME encourages making many matches
#' and also matching on covariates important to the outcome. The hyperparameter
#' \eqn{C} controls the balance between these two objectives. For more details,
#' please see the FLAME paper \href{https://arxiv.org/pdf/1707.06315.pdf}{here}.
#'
#' @section Stopping Rules:
#' By default, \code{FLAME} stops when 1. all covariates have been dropped or 2.
#' all treatment or control units have been matched. This behavior can be
#' modified by the arguments whose prefix is "early_stop". With the exception of
#' \code{early_stop_iterations}, all the rules come into play \emph{before} the
#' offending covariate is dropped. That is, if \code{early_stop_BF = 0.2} and at
#' the current iteration, dropping the covariate leading to highest match
#' quality is associated with a balancing factor of 0.1, FLAME stops
#' \emph{without} dropping this covariate.
#'
#' @section Missing Data:
#' \code{FLAME} offers functionality for handling missing data in the XXXXX, for
#' both the \code{data} and \code{holdout} sets. This functionality can be
#' specified via the arguments whose prefix is "missing". The simplest option,
#' for missingness in \code{data} is to set \code{missing_data = 1}, in which
#' case any units with missing data are dropped and not used in the algorithm.
#' The same goes for \code{holdout} and \code{holdout_data}. Missing values can
#' also be imputed via \code{mice::mice} by specifying \code{missing_data = 2}
#' and \code{missing_holdout = 2}, respectively. In this case,
#' \code{missing_data_imputations} many datasets will be imputed for \code{data}
#' (and similarly for \code{holdout}). If more than one imputation of
#' \code{data} is requested, the FLAME algorithm will be run on all imputations.
#' If more than one imputation of \code{holdout} is requested, the predictive
#' error at an iteration will be the average of predictive  errors across all
#' imputed \code{holdout} datasets.
#'
#' @param data Data to be matched. Either a dataframe or a path to a .csv file
#'   to be read (via \code{read.csv}) into a dataframe. Treatment must be
#'   described in a logical or binary column with name
#'   \code{treated_column_name}. Outcome must be either binary (logical or
#'   binary column) or continuous (numeric). All other columns will be assumed
#'   to be covariates to be used for matching and must be categorical as they
#'   will be coerced to factors. No default.
#' @param holdout Holdout data to be used to compute predictive error. If a
#'   numeric scalar between 0 and 1, that proportion of \code{data} will be made
#'   into a holdout set and only the remaining proportion of \code{data} will be
#'   matched. Otherwise, a dataframe or a path to a .csv file. Must have the
#'   same column names as \code{data}. This data will \emph{not} be matched.
#'   Defaults to 0.1
#' @param C A finite, positive scalar denoting the tradeoff between BF and PE in
#'   the FLAME algorithm. Higher C prioritizes more matches and lower C
#'   prioritizes not dropping important covariates. Defaults to 0.1.
#' @param treated_column_name A character with the name of the treatment column
#'   in \code{data} and \code{holdout}. Defaults to 'treated'.
#' @param outcome_column_name A character with the name of the outcome column in
#'   \code{data} and \code{holdout}. Defaults to 'outcome'.
#' @param binning_method The method to be used to bin continuous covariates in the
#' data. One of: "sturges", "scott", or "fd", denoting Sturges' rule, Scott's rule,
#' or the Freedman-Diaconis rule for determining number of bins in a histogram.
#' Each continuous covariate will be binned into the corresponding number of bins.
#' @param PE_method Either "elasticnet" or "xgb". Denotes the method to be used
#'   to compute PE. If "elasticnet", uses \code{glmnet::cv.glmnet} with default
#'   parameters and then the default predict method to estimate the outcome. If
#'   "xgb", uses \code{xgboost::xgb.cv} on a wide range of parameter values to
#'   cross-validate and find the best with respect to RMSE (for continuous
#'   outcomes) or binary misclassification rate (for binary outcomes). Then uses
#'   default predict method to estimate the outcome. Defaults to "elasticnet".
#' @param user_PE_fit An optional function supplied by the user that can be used
#'   instead of those allowed for by \code{PE_method} to fit a model fitting the
#'   outcome from the covariates. Must take in a matrix of covariates as its
#'   first argument and a vector outcome as its second argument. Defaults to
#'   \code{NULL}.
#' @param user_PE_fit_params A named list of optional parameters to be used by
#'   \code{user_PE_fit}. Defaults to \code{NULL}.
#' @param user_PE_predict An optional function supplied by the user that can be
#'   used to generate predictions from the output of \code{user_PE_fit}. If not
#'   supplied, defaults to \code{predict}. Defaults to \code{NULL}.
#' @param user_PE_predict_params A named list of optional parameters to be used
#'   by \code{user_PE_params}. Defaults to \code{NULL}.
#' @param replace A logical scalar. If \code{TRUE}, allows the same unit to be
#'   matched multiple times, on different sets of covariates. Defaults to
#'   \code{FALSE}.
#' @param verbose Controls output while FLAME is running. If 0, no output. If 1,
#'   outputs the iteration every iteration and the stopping condition. If 2,
#'   outputs the iteration and number of unmatched units every 5 iterations, and
#'   the stopping condition. If 3, outputs the iteration and number of unmatched
#'   units every 5 iterations, and the stopping condition. Defaults to 2.
#' @param want_pe A logical scalar. If TRUE, the predictive error (PE) at each
#'   iteration will be returned. Defaults to \code{FALSE}.
#' @param want_bf A logical scalar. If TRUE, the balancing factor (BF) at each
#'   iteration will be returned. Defaults to \code{FALSE}.
#' @param early_stop_iterations A nonnegative integer, denoting the number of
#'   iterations of FLAME to be performed. If 0, one round of exact matching is
#'   performed before stopping. Defaults to \code{Inf}.
#' @param early_stop_epsilon A nonnegative numeric. If FLAME attemts to drop a
#'   covariate that would raise the PE above (1 + early_stop_epsilon) times the
#'   baseline PE (the PE ebefore any covariates have been dropped), FLAME will
#'   stop. Defaults to 0.25.
#' @param early_stop_un_c_frac A numeric value between 0 and 1 (inclusive). If
#'   the proportion of control units that are unmatched falls below this value,
#'   FLAME stops. Defaults to 0.
#' @param early_stop_un_t_frac A numeric value between 0 and 1 (inclusive). If
#'   the proportion of treatment units that are unmatched falls below this
#'   value, FLAME stops. Defaults to 0.
#' @param early_stop_pe A numeric value between 0 and 1 (inclusive). If FLAME
#'   attempts to drop a covariate that would lead to a PE above this value,
#'   FLAME stops. Defaults to \code{Inf}.
#' @param early_stop_bf A numeric value between 0 and 1 (inclusive). If FLAME
#'   attempts to drop a covariate that would lead to a BF below this value,
#'   FLAME stops. Defaults to 0.
#' @param missing_data If 0, assumes no missingness in \code{data}. If 1,
#'   eliminates units with missingness from \code{data}. If 2, performs
#'   \code{missing_data_imputations} of MICE to impute the missing data. In this
#'   case, FLAME will be run on each imputed dataset and a list of the results
#'   for each dataset will be returned. If 3, will not match a unit on a
#'   covariate that it is missing. Defaults to 0.
#' @param missing_holdout If 0, assumes no missing data in \code{holdout}. If 1,
#'   eliminates units with missingness from \code{holdout}. If 2, performs
#'   \code{missing_holdout_imputations} of MICE to impute the missing data. In
#'   this latter case, all imputations will be used to compute PE, and the PE at
#'   an iteration will be the average across all imputations. Defaults to 0.
#' @param missing_holdout_imputations If \code{missing_holdout} = 2, performs
#'   this many imputations of the missing data in \code{holdout} using MICE.
#'   Defaults to 5.
#' @param missing_data_imputations If \code{missing_data} = 2, performs this
#'   many imputations of the missing data in \code{data} using MICE. Defaults to
#'   5.

#' @examples
#' data <- gen_data()
#' holdout <- gen_data()
#' FLAME_out <- FLAME(data = data, holdout = holdout)
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @importFrom zeallot %<-%
#' @importFrom rlang !!
#' @export
FLAME <-
  function(data, holdout = 0.1, C = 0.1,
           treated_column_name = 'treated', outcome_column_name = 'outcome',
           binning_method = 'sturges',
           PE_method = 'elasticnet',
           user_PE_fit = NULL, user_PE_fit_params = NULL,
           user_PE_predict = NULL, user_PE_predict_params = NULL,
           replace = FALSE, verbose = 2, want_pe = FALSE, want_bf = FALSE,
           early_stop_iterations = Inf, early_stop_epsilon = 0.25,
           early_stop_un_c_frac = 0, early_stop_un_t_frac = 0,
           early_stop_pe = Inf, early_stop_bf = 0,
           missing_data = 0, missing_holdout = 0,
           missing_data_imputations = 5, missing_holdout_imputations = 5) {

  read_data_out <- read_data(data, holdout)
  data <- read_data_out[[1]]
  holdout <- read_data_out[[2]]

  check_args(data, holdout, C,
             treated_column_name, outcome_column_name,
             binning_method,
             PE_method, user_PE_fit, user_PE_fit_params,
             user_PE_predict, user_PE_predict_params,
             replace, verbose, want_pe, want_bf,
             early_stop_iterations, early_stop_epsilon,
             early_stop_un_c_frac, early_stop_un_t_frac,
             early_stop_pe, early_stop_bf,
             missing_data, missing_holdout,
             missing_data_imputations, missing_holdout_imputations)

  missing_out <-
    handle_missing_data(data, holdout,
                        missing_data, missing_holdout,
                        missing_data_imputations, missing_holdout_imputations)

  data <- missing_out[[1]]
  holdout <- missing_out[[2]]
  is_missing <- missing_out[[3]]

  c(data, covs, n_covs, n_levels, cov_names, sorting_order) %<-%
    sort_cols(data, treated_column_name, outcome_column_name,
              binning_method, type = 'data', is_missing)

  holdout <-
    sort_cols(holdout, treated_column_name, outcome_column_name,
              binning_method, type = 'holdout')[[1]]

  n_iters <- length(data)

  browser()

  FLAME_out <- vector(mode = 'list', length = n_iters)
  for (i in 1:n_iters) {
    if (missing_data == 2) {
      message('Running FLAME on imputed dataset ', i, ' of ', n_iters)
    }
    FLAME_out[[i]] <-
      FLAME_internal(data[[i]], holdout, covs, n_covs, n_levels,
                     cov_names, sorting_order, C,
                     PE_method, user_PE_fit, user_PE_fit_params,
                     user_PE_predict, user_PE_predict_params,
                     replace, verbose, want_pe, want_bf,
                     early_stop_iterations, early_stop_epsilon,
                     early_stop_un_c_frac, early_stop_un_t_frac,
                     early_stop_pe, early_stop_bf)
  }

  if (n_iters == 1) {
    return(FLAME_out[[1]])
  }
  return(FLAME_out)
}

FLAME_internal <-
  function(data, holdout, covs, n_covs, n_levels, cov_names, sorting_order, C,
           PE_method, user_PE_fit, user_PE_fit_params ,
           user_PE_predict, user_PE_predict_params,
           replace, verbose, want_pe, want_bf,
           early_stop_iterations, early_stop_epsilon,
           early_stop_un_c_frac, early_stop_un_t_frac,
           early_stop_pe, early_stop_bf) {

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
  processed_matches <-
    process_matches(data, replace, covs, n_levels, MGs,
                    matched_on, matching_covs, CATE, cov_names)

  CATE <- processed_matches[[1]]
  MGs <- processed_matches[[2]]
  matched_on <- processed_matches[[3]]
  units_matched <- processed_matches[[4]]
  made_matches <- processed_matches[[5]]

  if (made_matches) {
    data$matched[units_matched] <- TRUE
    data$weight[units_matched] <- data$weight[units_matched] + 1
    matching_covs %<>% c(list(order_cov_names(cov_names[covs],
                                              cov_names,
                                              sorting_order)))
  }
  store_pe <- NULL
  store_bf <- NULL

  iter <- 0
  baseline_PE <- get_PE(cov_to_drop = NULL, covs, holdout,
                        PE_method, user_PE_fit, user_PE_fit_params,
                        user_PE_predict, user_PE_predict_params)

  while (!early_stop(iter, data, covs, early_stop_iterations, verbose)) {
    iter <- iter + 1
    show_progress(verbose, iter, data)

    # Compute the PE associated with dropping each covariate
    PE <- sapply(covs, get_PE, covs, holdout,
                 PE_method, user_PE_fit, user_PE_fit_params,
                 user_PE_predict, user_PE_predict_params)

    if (early_stop_PE(min(PE), early_stop_pe,
                      early_stop_epsilon, baseline_PE, verbose)) {
      break
    }

    if (C != 0) {
      best_lower_bound <- max(-PE)
      upper_bound <- 2 * C - PE

      drop_candidates <- which(upper_bound >= best_lower_bound)
      PE <- PE[drop_candidates]

      BF_out <-
        lapply(covs[drop_candidates], get_BF, data, replace, covs, n_levels)
      BF <- sapply(BF_out, function(x) x[['BF']])

      MQ <- C * BF - PE
    }
    else {
      MQ <- PE
    }

    # (First, in unlikely case of ties) covariate yielding highest match quality
    drop <- which.max(MQ)
    prop_unmatched <- BF_out[[drop]][['prop_unmatched']]

    if (early_stop_BF(BF[drop], early_stop_bf,
                      prop_unmatched[['control']], prop_unmatched[['treated']],
                      early_stop_un_c_frac, early_stop_un_t_frac,
                      verbose)) {
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
    processed_matches <-
      process_matches(data, replace, covs, n_levels, MGs,
                      matched_on, matching_covs, CATE, cov_names)
    CATE <- processed_matches[[1]]
    MGs <- processed_matches[[2]]
    matched_on <- processed_matches[[3]]
    units_matched <- processed_matches[[4]]
    made_matches <- processed_matches[[5]]

    if (made_matches) {
      if (replace) {
        # Only use * to refer to main matched group
        # weight == 0 implies never matched before so this match is their MMG
        data[intersect(units_matched, which(data$weight == 0)),
                       setdiff(1:n_covs, covs)] <- '*'
      }
      else {
        data[units_matched, setdiff(1:n_covs, covs)] <- '*'
      }
      data$matched[units_matched] <- TRUE
      data$weight[units_matched] <- data$weight[units_matched] + 1
    }
  }

  # Done matching!

  # Substitute covariate values of all unmatched units with the unmatched
  # covariate symbol '*'
  data[!data$matched, 1:n_covs] <- '*'

  # Reorder the data according to the original column order
  data[, 1:n_covs] %<>% dplyr::select(order(sorting_order))
  data[, ncol(data)] <- NULL
  colnames(data) <-
    c(colnames(data)[1:n_covs][order(sorting_order)],
      'outcome', 'treated', 'matched', 'weight')

  if (!replace) {
    matched <- which(data$matched)
    stopifnot(all(data$weight[matched] == 1))
    stopifnot(all(data$weight[-matched] == 0))
  }

  ret_list <-
    list(data = data,
         MGs = MGs,
         CATE = CATE,
         matched_on = matched_on,
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
