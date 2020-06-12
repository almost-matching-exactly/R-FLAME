aggregate_table <- function(vals) {
  vals <- as.character(vals)
  tab <- table(vals)
  name <- names(tab)
  return(as.vector(tab[match(vals, name)]))
}

# bit_match takes a data frame, a set of covariates to match on, the
# treatment indicator column and the matched indicator column. it returns the
# array indicating whether each unit is matched (the first return value), and a
# list of indices for the matched units (the second return value)

bit_match <- function(data, covs) {
  # gmp::as.bigz is for handling lots and lots of covariates so we don't
  # have trouble with overflow issues
  n_levels <- sapply(data[, covs, drop = FALSE], nlevels) - 1
  data_wo_t <- gmp::as.bigz(as.matrix(data[, covs[order(n_levels)]]))
  n_levels <- sort(n_levels) # isn't this already sorted?

  # Compute b_u
  multiplier <- gmp::pow.bigz(n_levels, seq_along(n_levels) - 1)

  b_u <- as.vector(gmp::`%*%`(data_wo_t, multiplier))

  # Compute b_u+
  multiplier <- gmp::pow.bigz(n_levels, seq_along(n_levels))

  b_u_plus <-
    as.vector(gmp::add.bigz(gmp::`%*%`(data_wo_t, multiplier), data$treated))

  # Compute c_u
  c_u <- aggregate_table(b_u)

  # Compute c_u+
  c_u_plus <- aggregate_table(b_u_plus)

  match_index <- (c_u != c_u_plus) & (c_u >= 2)
  index <- b_u[match_index]

  return(list(match_index = match_index,
              index = index))
}

make_MGs <-
  function(data, outcome_in_data, index, matched_units, covs, cov_names) {
  # Takes all the units that were matched on these p covariates and separates
  # them into matched groups based off their unique values of those covariates
  # Returns a list with three items
    # MGs: a list, each entry of which corresponds to a different MG and
    #   contains the indices of the corresponding units
    # CATEs: a vector the same length as MGs, each entry of which is the CATE
    #   for the corresponding MG
    # matched_on: a list the same length as MGs, each entry of which is a named
    #   vector detailing the covariates and their values that the units in the
    #   corresponding MG matched on

  unique_MGs <- unique(index)
  n_MGs <- length(unique_MGs)

  # MGs <- vector('list', length = n_MGs)
  MGs <- list()
  # CATEs <- vector('numeric', length = n_MGs)
  # matched_on <- vector('list', length = n_MGs)

  for (i in 1:n_MGs) {
    members <- matched_units[which(index == unique_MGs[i])]
    if (any(data$weight[members] == 0)) { # Only store new MGs if replace = TRUE
      MGs <- c(MGs, list(members))
    }
    # treated <- intersect(members, which(data$treated == 1))
    # control <- intersect(members, which(data$treated == 0))
    # if (outcome_in_data) {
    #   CATEs[i] <- mean(data$outcome[treated]) - mean(data$outcome[control])
    # }
    # else {
    #   CATEs[i] <- NA
    # }

    # matched_on[[i]] <- data[members[1], covs, drop = FALSE]
    # rownames(matched_on[[i]]) <- NULL
    #
    # names(matched_on[[i]]) <- cov_names[covs]
  }
  return(list(MGs = MGs))
  # return(list(MGs = MGs,
  #             CATEs = CATEs,
  #             matched_on = matched_on))
}

# process_matches <-
#   function(data, outcome_in_data, replace, covs, MGs,
#            matched_on, matching_covs, CATE, cov_names) {
process_matches <-
  function(data, outcome_in_data, replace, covs, MGs,
           matching_covs, cov_names) {
  # Make matches implied by covs in data
    # Store any relevant information

  # Probably cleaner to relegate this if-statement (and ones like it elsewhere in the code)
    # to the bit_match function
  if (replace) {
    match_out <- bit_match(data[!data$missing, ], covs)
    match_index <- match_out[[1]]
    index <- match_out[[2]]
    units_matched <- which(!data$missing)[match_index]
  }
  else {
    match_out <-
      bit_match(data[!data$matched & !data$missing, ], covs)
    match_index <- match_out[[1]]
    index <- match_out[[2]]

    # Adjusts for fact that match_index returned by update_matched bit is
    # with respect to the unmatched subset of data
    units_matched <- which(!data$matched & !data$missing)[match_index]
  }

  made_matches <- sum(match_index) > 0

  if (made_matches) {
    new_MGs <-
      make_MGs(data, outcome_in_data, index, units_matched, covs, cov_names)
    MGs <- c(MGs, new_MGs$MGs)
    # CATE <- c(CATE, new_MGs$CATEs)
    # matched_on <- c(matched_on, new_MGs$matched_on)
  }
  # matching_covs[[length(matching_covs) + 1]] <- cov_names[covs]
  return(list(MGs = MGs,
              units_matched = units_matched,
              made_matches = made_matches))
  # return(list(CATE = CATE,
  #             MGs = MGs,
  #             matched_on = matched_on,
  #             units_matched = units_matched,
  #             made_matches = made_matches))
}

get_PE <- function(cov_to_drop, covs, holdout, PE_method,
                   user_PE_fit, user_PE_fit_params,
                   user_PE_predict, user_PE_predict_params) {

  PE_predict_params <- list()

  if (!is.null(user_PE_fit)) {
    PE_fit <- user_PE_fit
    PE_fit_params <- user_PE_fit_params
  }
  else {
    if (PE_method == 'ridge') {
      PE_fit <- glmnet::cv.glmnet
      if (length(unique(holdout[[1]]$outcome)) == 2) {
        family <- 'binomial'
        PE_predict_params <- list(type = 'class')
      }
      else if (is.factor(holdout[[1]]$outcome)) {
        family <- 'multinomial'
        PE_predict_params <- list(type = 'class')
      }
      else {
        family <- 'gaussian'
      }
      PE_fit_params <- list(family = family, nfolds = 5)
    }
    else if (PE_method == 'xgb') {
      PE_fit <- cv_xgboost
      if (length(unique(holdout[[1]]$outcome)) == 2) {
        obj <- 'binary:logistic'
      }
      else if (is.factor(holdout[[1]]$outcome)) {
        obj <- 'multi:softmax'
      }
      else {
        obj <- 'reg:squarederror'
      }
      PE_fit_params <- list(obj = obj)
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
  }

  PE <- predict_master(holdout, covs, cov_to_drop,
                       PE_fit, PE_predict, PE_fit_params, PE_predict_params)
  return(PE)
}

get_BF <- function(cov_to_drop, data, replace, covs) {

  # Calculate number of units eligible to be matched
  if (replace) {
    n_control <- sum(data$treated[!data$missing] == 0)
    n_treated <- sum(data$treated[!data$missing] == 1)
  }
  else {
    n_control <- sum(data$treated[!data$matched & !data$missing] == 0)
    n_treated <- sum(data$treated[!data$matched & !data$missing] == 1)
  }

  if (replace) {
    match_index <-
      bit_match(data[!data$missing, ], setdiff(covs, cov_to_drop))$match_index
    units_matched <- which(!data$missing)[match_index]
  }
  else {
    match_index <-
      bit_match(data[!data$matched & !data$missing, ],
                         setdiff(covs, cov_to_drop))$match_index
    units_matched <- which(!data$matched & !data$missing)[match_index]
  }

  # Newly matched
  n_control_matched <- sum(data$treated[units_matched] == 0)
  n_treated_matched <- sum(data$treated[units_matched] == 1)

  # All matched units; for stopping rule purposes
  all_unmatched <-
    setdiff(1:nrow(data), union(units_matched, which(data$matched)))

  n_control_unmatched <- sum(all_unmatched %in% which(data$treated == 0))
  n_treated_unmatched <- sum(all_unmatched %in% which(data$treated == 1))

  prop_control_unmatched <- n_control_unmatched / n_control
  prop_treated_unmatched <- n_treated_unmatched / n_treated

  # Is this if_else really necessary? We should catch it earlier
  if (n_control == 0 | n_treated == 0) {
    BF <- 0
  }
  else {
    BF <- n_control_matched / n_control + n_treated_matched / n_treated
  }

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
#' \eqn{BF} denotes the balancing factor, defined as the proportion of unmatched
#' control units, plus the proportion of unmatched treated units, that can now
#' be matched by dropping that covariate. And \eqn{PE} denotes the prediction
#' error, defined as the training error incurred when predicting the outcome
#' from covariates on a separate, holdout set. In this way, FLAME encourages
#' making many matches and also matching on covariates important to the outcome.
#' The hyperparameter \eqn{C} controls the balance between these two objectives.
#' For more details, please see the FLAME paper
#' \href{https://arxiv.org/pdf/1707.06315.pdf}{here}.
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
#' \code{FLAME} offers functionality for handling missing data in the
#' covariates, for both the \code{data} and \code{holdout} sets. This
#' functionality can be specified via the arguments whose prefix is "missing" or
#' "impute". It allows for ignoring missing data, imputing it, or (for
#' \code{data}) not matching on missing values. If \code{data} is imputed, the
#' FLAME algorithm will be run on all imputations. If \code{holdout} is imputed,
#' the predictive error at an iteration will be the average of predictive errors
#' across all imputed \code{holdout} datasets.
#'
#'
#' @param data Data to be matched. Either a data frame or a path to a .csv file
#'   to be read into a data frame. Treatment must be described by a logical or
#'   binary numeric column with name \code{treated_column_name}. Outcome must be
#'   described by a column with name \code{outcome_column_name}. The outcome
#'   will be treated as continuous if numeric, as binary if a two-level factor
#'   or numeric with two unique values, and as multi-class if a factor with more
#'   than two levels. The outcome column may be omitted, in which case matching
#'   will be performed but treatment effect estimation will not be possible. All
#'   columns not containing outcome or treatment will be treated as covariates
#'   for matching. Covariates are assumed to be categorical and will be coerced
#'   to factors, though they may be passed as either factors or numeric. If you wish
#'   to use continuous covariates for matching, they should be binned prior to being
#'   passed to \code{FLAME}. There is no default for \code{data}.
#' @param holdout Holdout data to be used to compute predictive error. If a
#'   numeric scalar between 0 and 1, that proportion of \code{data} will be made
#'   into a holdout set and only the remaining proportion of \code{data} will be
#'   matched. Otherwise, a data frame or a path to a .csv file. The holdout data must contain
#'   an outcome column with name \code{outcome_column_name}; other restrictions on
#'   column types are as for \code{data}. Covariate columns must have the same column
#'   names and order as \code{data}. This data will \emph{not} be matched.
#'   Defaults to 0.1.
#' @param C A finite, positive scalar denoting the tradeoff between BF and PE in
#'   the FLAME algorithm. Higher C prioritizes more matches and lower C
#'   prioritizes not dropping important covariates. Defaults to 0.1.
#' @param treated_column_name A character with the name of the treatment column
#'   in \code{data} and \code{holdout}. Defaults to 'treated'.
#' @param outcome_column_name A character with the name of the outcome column in
#'   \code{holdout} and also in \code{data}, if supplied in the latter.
#'   Defaults to 'outcome'.
#' @param PE_method Either "ridge" or "xgb". Denotes the method to be used to
#'   compute PE. If "ridge", uses \code{glmnet::cv.glmnet} with default
#'   parameters and then the default predict method to estimate the outcome. If
#'   "xgb", uses \code{xgboost::xgb.cv} on a wide range of parameter values to
#'   cross-validate and find the best with respect to RMSE (for continuous
#'   outcomes) or misclassification rate (for binary/multi-class outcomes). Then uses
#'   the default \code{predict} method to estimate the outcome. Defaults to "ridge".
#' @param user_PE_fit An optional function supplied by the user that can be used
#'   instead of those allowed for by \code{PE_method} to fit a model for the
#'   outcome from the covariates. Must take in a matrix of covariates as its
#'   first argument and a vector outcome as its second argument. Defaults to
#'   \code{NULL}.
#' @param user_PE_fit_params A named list of optional parameters to be used by
#'   \code{user_PE_fit}. Defaults to \code{NULL}.
#' @param user_PE_predict An optional function supplied by the user that can be
#'   used to generate predictions from the output of \code{user_PE_fit}. As its
#'   first argument, must take an object of the type returned by
#'  \code{user_PE_fit} and as its second, a matrix of values for which to
#'  generate predictions. When the outcome is binary or multi-class, must
#'  return the maximum probability class label. If not supplied,
#'  defaults to \code{predict}.
#' @param user_PE_predict_params A named list of optional parameters to be used
#'   by \code{user_PE_predict}. Defaults to \code{NULL}.
#' @param replace A logical scalar. If \code{TRUE}, allows the same unit to be
#'   matched multiple times, on different sets of covariates. In this case,
#'   balancing factor is computing by dividing by the total number of treatment
#'   (control) units, instead of the number of unmatched treatment (control)
#'   units. Defaults to \code{FALSE}.
#' @param verbose Controls how FLAME displays progress while running. If 0, no
#'   output. If 1, only outputs the stopping condition. If 2, outputs the
#'   iteration and number of unmatched units every 5 iterations, and the
#'   stopping condition. If 3, outputs the iteration and number of unmatched
#'   units every iteration, and the stopping condition. Defaults to 2.
#' @param return_pe A logical scalar. If \code{TRUE}, the predictive error (PE)
#'   at each iteration will be returned. Defaults to \code{FALSE}.
#' @param return_bf A logical scalar. If \code{TRUE}, the balancing factor (BF)
#'   at each iteration will be returned. Defaults to \code{FALSE}.
#' @param early_stop_iterations A nonnegative integer, denoting an upper bound
#'   on the number of iterations of FLAME to be performed. If 0, one round of
#'   exact matching is performed before stopping. Defaults to \code{Inf}.
#' @param early_stop_epsilon A nonnegative numeric. If FLAME attemts to drop a
#'   covariate that would raise the PE above (1 + early_stop_epsilon) times the
#'   baseline PE (the PE before any covariates have been dropped), FLAME will
#'   stop. Defaults to 0.25.
#' @param early_stop_control A numeric value between 0 and 1. If
#'   the proportion of control units that are unmatched falls below this value,
#'   FLAME stops. Defaults to 0.
#' @param early_stop_treated A numeric value between 0 and 1. If
#'   the proportion of treatment units that are unmatched falls below this
#'   value, FLAME stops. Defaults to 0.
#' @param early_stop_pe A numeric value between 0 and 1. If FLAME
#'   attempts to drop a covariate that would lead to a PE above this value,
#'   FLAME stops. Defaults to \code{Inf}.
#' @param early_stop_bf A numeric value between 0 and 1. If FLAME
#'   attempts to drop a covariate that would lead to a BF below this value,
#'   FLAME stops. Defaults to 0.
#' @param missing_data If 0, assumes no missingness in \code{data}. If 1, does
#'   not match units with missingness in \code{data}. In this case, the
#'   balancing factor is computed ignoring units with missingness. If 2,
#'   generates \code{missing_data_imputations} imputed datasets via
#'   \code{mice::mice}. In this case, the results of running \code{FLAME} on
#'   each imputed dataset will be returned in a list. Within each of these list
#'   entries, the \code{data} entry will contain the imputed, not missing,
#'   values. If 3, will not match a unit on a covariate that it is missing.
#'   Defaults to 0.
#' @param missing_holdout If 0, assumes no missing data in \code{holdout}. If 1,
#'   eliminates units with missingness from \code{holdout}. If 2, generates
#'   \code{missing_holdout_imputations} imputed datasets via \code{mice::mice}.
#'   In this latter case, all imputations will be used to compute PE, and the PE
#'   at an iteration will be the average across all imputations. Defaults to 0.
#' @param missing_holdout_imputations If \code{missing_holdout} = 2, performs
#'   this many imputations of the missing data in \code{holdout} via
#'   \code{mice::mice}. Defaults to 5.
#' @param missing_data_imputations If \code{missing_data} = 2, performs this
#'   many imputations of the missing data in \code{data} via \code{mice::mice}.
#'   Defaults to 5.
#' @param impute_with_treatment A logical scalar. If \code{TRUE}, uses treatment
#'   assignment to impute covariates when \code{missing_data = 2} or
#'   \code{missing_holdout = 2}. Defaults to \code{TRUE}.
#' @param impute_with_outcome A logical scalar. If \code{TRUE}, uses outcome
#'   information to impute covariates when \code{missing_data = 2} or
#'   \code{missing_holdout = 2}. Defaults to \code{FALSE}.
#'
#' @return The basic object returned by \code{FLAME} is a list of 6 entries:
#' \describe{
#' \item{data}{The original data frame with several modifications:
#'   \enumerate{
#'     \item An extra logical column, \code{data$matched},
#'     that indicates whether or not a unit was matched.
#'     \item An extra numeric column, \code{data$weight},
#'     that denotes on how many different sets of covariates a unit was matched.
#'     This will only be greater than 1 when \code{replace = TRUE}.
#'     \item Regardless of their original names, the columns denoting treatment
#'     and outcome in the data will be renamed 'treated' and 'outcome' and they
#'     are moved to be located after all the covariate data.
#'     \item Units that were not matched on all covariates will have a *
#'     in place of their covariate value for all covariates on which they
#'     were not matched.
#'     }
#'  }
#'  \item{MGs}{A list of all the matched groups formed by FLAME. Each entry
#'  contains the units in a single matched group}
#'  \item{CATE}{A numeric vector with the conditional average treatment effect
#'    of every matched group in \code{MGs}. Returned only if the outcome is
#'    numeric.}
#'  \item{matched_on}{A list corresponding to \code{MGs} that gives the
#'  covariates, and their values, on which units in each matched group were
#'  matched.}
#'  \item{matching_covs}{A list with the covariates used for matching on every
#'  iteration of FLAME}
#'  \item{dropped}{A vector with the covariate dropped at each iteration of
#'  FLAME}
#' }
#'
#' @examples
#' data <- gen_data()
#' holdout <- gen_data()
#' FLAME_out <- FLAME(data = data, holdout = holdout)
#' @importFrom stats model.matrix predict rbinom rnorm var
#' @importFrom utils flush.console read.csv write.csv
#' @importFrom devtools load_all
#' @export
FLAME <-
  function(data, holdout = 0.1, C = 0.1,
           treated_column_name = 'treated', outcome_column_name = 'outcome',
           PE_method = 'ridge',
           user_PE_fit = NULL, user_PE_fit_params = NULL,
           user_PE_predict = NULL, user_PE_predict_params = NULL,
           replace = FALSE, verbose = 2, return_pe = FALSE, return_bf = FALSE,
           early_stop_iterations = Inf, early_stop_epsilon = 0.25,
           early_stop_control = 0, early_stop_treated = 0,
           early_stop_pe = Inf, early_stop_bf = 0,
           missing_data = 0, missing_holdout = 0,
           missing_data_imputations = 5, missing_holdout_imputations = 5,
           impute_with_treatment = TRUE, impute_with_outcome = FALSE) {

  # Get matching and holdout data, from input .csv files, if necessary
  read_data_out <-
    read_data(data, holdout, treated_column_name, outcome_column_name)

  data <- read_data_out[[1]]
  holdout <- read_data_out[[2]]

  # Was outcome supplied by user (in matching data)?
  # If so, CATEs will be computed as FLAME runs
  if (!(outcome_column_name %in% colnames(data))) {
    outcome_in_data <- FALSE
  }
  else {
    outcome_in_data <- TRUE
  }

  # Make sure the user didn't do anything funny
  check_args(data, holdout, outcome_in_data, C,
             treated_column_name, outcome_column_name,
             PE_method, user_PE_fit, user_PE_fit_params,
             user_PE_predict, user_PE_predict_params,
             replace, verbose, return_pe, return_bf,
             early_stop_iterations, early_stop_epsilon,
             early_stop_control, early_stop_treated,
             early_stop_pe, early_stop_bf,
             missing_data, missing_holdout,
             missing_data_imputations, missing_holdout_imputations,
             impute_with_outcome, impute_with_treatment)

  # Map everything to factor
  cov_inds_data <- which(!(colnames(data) %in% c(treated_column_name, outcome_column_name)))
  data[, cov_inds_data] <- lapply(data[, cov_inds_data, drop = FALSE], as.factor)

  # Number of levels of each covariate
  n_levels <- sapply(data[, cov_inds_data, drop = FALSE], nlevels)
  if (any(n_levels > nrow(data) / 10)) {
    warning(paste("It looks like the variables {",
                  colnames(data)[cov_inds_data][n_levels > nrow(data) / 10],
                  "} might be continuous; you probably won't make many matches on them.",
                  "If this is the case, please either bin variables prior to",
                  "calling FLAME or do not include them at all."))
  }

  cov_inds_holdout <- which(!(colnames(holdout) %in% c(treated_column_name, outcome_column_name)))
  holdout[, cov_inds_holdout] <- lapply(holdout[, cov_inds_holdout, drop = FALSE], as.factor)

  # browser()
  # Map the factor levels of the covariates as supplied by the user to 0:k
  #   so that bit matching works properly.
  #   Store the mapping so you can return the data in the original format.
  remapped_data <- factor_remap(data, treated_column_name, outcome_column_name)
  data <- remapped_data$df

  mapping <- remapped_data$mapping

  # remapped_holdout <-
  #   factor_remap(holdout, treated_column_name, outcome_column_name, mapping)
  # holdout <- remapped_holdout$df

  # browser()
  # Impute missing data, if requested, else, prepare to deal with missingness
  #   as specified by missing_data
  #   Returns a LIST of data frames (see definition of n_iters)
  missing_out <-
    handle_missing_data(data, holdout, outcome_in_data,
                        treated_column_name, outcome_column_name,
                        missing_data, missing_holdout,
                        missing_data_imputations, missing_holdout_imputations,
                        impute_with_treatment, impute_with_outcome)

  data <- missing_out[[1]]
  holdout <- missing_out[[2]]
  is_missing <- missing_out[[3]]

  # Put the columns in increasing order of number of levels so that bit-matching
  #  works. Should maybe (probably?) be combined with factor_remap
  sort_cols_out <-
    sort_cols(data, outcome_in_data, treated_column_name, outcome_column_name,
              type = 'data', is_missing)

  data <- sort_cols_out[[1]]
  covs <- sort_cols_out[[2]]
  n_covs <- sort_cols_out[[3]]
  cov_names <- sort_cols_out[[4]]

  # Also sort holdout so that data and holdout columns correspond to one another
  holdout <-
    sort_cols(holdout, outcome_in_data = TRUE,
              treated_column_name, outcome_column_name,
              type = 'holdout')[[1]]

  # data is now a list of data frames so as to accomodate multiple imputations
  n_iters <- length(data)

  # For each imputed data set (just 1 if no missingness), run FLAME
  FLAME_out <- vector(mode = 'list', length = n_iters)
  for (i in 1:n_iters) {
    if (missing_data == 2) {
      message('Running FLAME on imputed dataset ', i, ' of ', n_iters, '\r',
              appendLF = FALSE)
      flush.console()
    }
    FLAME_out[[i]] <-
      FLAME_internal(data[[i]], outcome_in_data,
                     holdout, covs, n_covs,
                     cov_names, C,
                     PE_method, user_PE_fit, user_PE_fit_params,
                     user_PE_predict, user_PE_predict_params,
                     replace, verbose, return_pe, return_bf,
                     early_stop_iterations, early_stop_epsilon,
                     early_stop_control, early_stop_treated,
                     early_stop_pe, early_stop_bf, mapping)
  }

  # If 0 or 1 imputations, don't return a list
  if (n_iters == 1) {
    return(FLAME_out[[1]])
  }
  return(FLAME_out)
}

FLAME_internal <-
  function(data, outcome_in_data, holdout, covs, n_covs,
           cov_names, C,
           PE_method, user_PE_fit, user_PE_fit_params ,
           user_PE_predict, user_PE_predict_params,
           replace, verbose, return_pe, return_bf,
           early_stop_iterations, early_stop_epsilon,
           early_stop_control, early_stop_treated,
           early_stop_pe, early_stop_bf, mapping) {

  # The lines after processed_matches suggest these several lines should be removed.
  #   To do.
  # List of MGs, each entry contains the corresponding MG's entries
  MGs <- list()
  # List of CATEs, each entry contains the corresponding MG's CATE
  # CATE <- vector('numeric')
  # List of covariates and their values matched on for each corresponding MG
  # matched_on <- list()

  # List of covariates used to match at each level
  matching_covs <- list()
  covs_dropped <- NULL

  # Try and make matches on all covariates
  # processed_matches <-
  #   process_matches(data, outcome_in_data, replace, covs, MGs,
  #                   matched_on, matching_covs, CATE, cov_names)
  processed_matches <-
    process_matches(data, outcome_in_data, replace, covs, MGs,
                    matching_covs, cov_names)

  # CATE <- processed_matches[[1]]
  MGs <- processed_matches[[1]]
  # matched_on <- processed_matches[[3]]
  # units_matched <- processed_matches[[4]]
  units_matched <- processed_matches[[2]]
  # made_matches <- processed_matches[[5]]
  made_matches <- processed_matches[[3]]

  if (made_matches) {
    data$matched[units_matched] <- TRUE
    data$weight[units_matched] <- data$weight[units_matched] + 1
    matching_covs <- c(matching_covs, list(cov_names[covs]))
  }
  store_pe <- NULL
  store_bf <- NULL

  # Predictive error using all covariates. Used for stopping condition.
  baseline_PE <- get_PE(cov_to_drop = NULL, covs, holdout,
                        PE_method, user_PE_fit, user_PE_fit_params,
                        user_PE_predict, user_PE_predict_params)
  iter <- 0

  while (!early_stop(iter, data, covs, early_stop_iterations, verbose)) {
    iter <- iter + 1
    # Progress messages
    show_progress(verbose, iter, data)

    # Compute the PE associated with dropping each covariate
    PE <- sapply(covs, get_PE, covs, holdout,
                 PE_method, user_PE_fit, user_PE_fit_params,
                 user_PE_predict, user_PE_predict_params)

    ## Note for Vittorio: min(PE)?
    if (early_stop_PE(min(PE), early_stop_pe,
                      early_stop_epsilon, baseline_PE, verbose)) {
      break
    }

    # Because 0 < BF < 2, we have that -PE < MQ < 2 * C - PE
    #   Thus if 2 * C - PE associated with a covariate X is not higher than
    #   the highest -PE across all covariates, we'll never end up dropping X,
    #   no matter the associated BF. So we can simply not compute it.
    #   This scenario comes (e.g.) up when you have "pretty irrelevant" covariates because the
    #   maximum -PE will be quite large so we'll never even consider dropping "pretty relevant" covariates.
    best_lower_bound <- max(-PE)
    upper_bound <- 2 * C - PE

    # Based off the above, the covariates we'll consider dropping.
    drop_candidates <- which(upper_bound >= best_lower_bound)

    PE <- PE[drop_candidates]
    BF_out <-
      lapply(covs[drop_candidates], get_BF, data, replace, covs)
    BF <- sapply(BF_out, function(x) x[['BF']])

    MQ <- C * BF - PE

    # (First, in unlikely case of ties) covariate yielding highest MQ
    drop <- which.max(MQ)
    prop_unmatched <- BF_out[[drop]][['prop_unmatched']]

    if (early_stop_BF(BF[drop], early_stop_bf,
                      prop_unmatched[['control']], prop_unmatched[['treated']],
                      early_stop_control, early_stop_treated,
                      verbose)) {
      break
    }

    # Store the PE and BF of the associated dropped covariate, to return if desired
    store_pe <- c(store_pe, PE[drop])
    store_bf <- c(store_bf, BF[drop])

    covs_dropped <- c(covs_dropped, cov_names[covs[drop_candidates[drop]]])
    covs <- covs[-drop_candidates[drop]]

    matching_covs <- c(matching_covs, list(cov_names[covs]))

    # Make new matches having dropped a covariate
    ## Ideally should just return this from MQ so you don't have to redo it
    # processed_matches <-
    #   process_matches(data, outcome_in_data, replace, covs, MGs,
    #                   matched_on, matching_covs, CATE, cov_names)
    processed_matches <-
      process_matches(data, outcome_in_data, replace, covs, MGs,
                      matching_covs, cov_names)
    # CATE <- processed_matches[[1]]
    MGs <- processed_matches[[1]]
    # matched_on <- processed_matches[[3]]
    units_matched <- processed_matches[[2]]
    made_matches <- processed_matches[[3]]

    if (made_matches) {
      if (replace) { # Note for VittoriO: shouldn't need the if-else bc data$weight equivalent with unmatched for !replace
        # Only use * to refer to main matched group
        # weight == 0 implies never matched before so this match is their MMG
        data[intersect(units_matched, which(data$weight == 0)),
                       setdiff(1:n_covs, covs)] <- '*' # can we just use dropped again here?
      }
      else {
        data[units_matched, setdiff(1:n_covs, covs)] <- '*'
      }
      data$matched[units_matched] <- TRUE
      data$weight[units_matched] <- data$weight[units_matched] + 1
    }
  }

  # Done matching!
  # is_missing <- data$missing

  # Substitute covariate values of all unmatched units with the unmatched
  # covariate symbol '*'
  data[!data$matched, 1:n_covs] <- '*'

  data[, ncol(data)] <- NULL

  # Swap back the original factor levels
  rev_mapping <- lapply(mapping, function(x) {
    tmp <- c(names(x), '*')
    names(tmp) <- c(x, '*')
    return(tmp)
  })

  data <-
    factor_remap(data, mapping = rev_mapping)$df

  # Undoes rownames conflicting with functions in post_matching.R
  # when holdout is taken from data
  rownames(data) <- 1:nrow(data)

  matched_on <- lapply(MGs, function(x) {
    cols_matched_on <- sapply(data[x, 1:n_covs, drop = FALSE], function(y) all(y != '*'))
    tmp <- data[x[1], 1:n_covs, drop = FALSE][1, cols_matched_on, drop = FALSE]
    # if (n_covs == 1) { # weird behavior where column is renamed if drop == F and n_covs == 1
    #   colnames(tmp) <- colnames(data)[1]
    # }
    rownames(tmp) <- NULL
    return(tmp)
  })

  ret_list <-
    list(data = data,
         MGs = MGs,
         matched_on = matched_on,
         matching_covs = matching_covs,
         dropped = covs_dropped)

  # implicitly checks that outcome exists in data; otherwise data$outcome returns NULL
  if (is.numeric(data$outcome)) {
    CATE <-
      sapply(MGs, function(x) {
        treated <- data$treated[x]
        outcomes <- data$outcome[x]
        mean(outcomes[treated == 1]) - mean(outcomes[treated == 0])
      })
    ret_list$CATE <- CATE
  }

  if (return_pe) {
    ret_list <- c(ret_list, 'PE' = list(store_pe))
  }
  if (return_bf) {
    ret_list <- c(ret_list, 'BF' = list(store_bf))
  }
  return(ret_list)
}
