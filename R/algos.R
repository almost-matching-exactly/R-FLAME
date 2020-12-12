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
#' @section Missing Data: \code{FLAME} offers functionality for handling missing
#'   data in the covariates, for both the \code{data} and \code{holdout} sets.
#'   This functionality can be specified via the arguments whose prefix is
#'   "missing" or "impute". It allows for ignoring missing data, imputing it, or
#'   (for \code{data}) not matching on missing values. If \code{data} is
#'   imputed, imputation will be done once and the FLAME algorithm will be run
#'   on all imputations. If \code{holdout} is imputed, the predictive error at
#'   an iteration will be the average of predictive errors across all imputed
#'   \code{holdout} datasets.
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
#'   to factors, though they may be passed as either factors or numeric. If you
#'   wish to use continuous covariates for matching, they should be binned prior
#'   to being passed to \code{FLAME}. There is no default for \code{data}.
#' @param holdout Holdout data to be used to compute predictive error. If a
#'   numeric scalar between 0 and 1, that proportion of \code{data} will be made
#'   into a holdout set and only the remaining proportion of \code{data} will be
#'   matched. Otherwise, a data frame or a path to a .csv file. The holdout data
#'   must contain an outcome column with name \code{outcome_column_name}; other
#'   restrictions on column types are as for \code{data}. Covariate columns must
#'   have the same column names and order as \code{data}. This data will
#'   \emph{not} be matched. Defaults to 0.1.
#' @param C A finite, positive scalar denoting the tradeoff between BF and PE in
#'   the FLAME algorithm. Higher C prioritizes more matches and lower C
#'   prioritizes not dropping important covariates. Defaults to 0.1.
#' @param treated_column_name A character with the name of the treatment column
#'   in \code{data} and \code{holdout}. Defaults to 'treated'.
#' @param outcome_column_name A character with the name of the outcome column in
#'   \code{holdout} and also in \code{data}, if supplied in the latter.
#'   Defaults to 'outcome'.
#' @param PE_method Denotes how predictive error (PE) is to be computed. Either
#'   a string -- one of "ridge" (default) or "xgb" -- or a function. If "ridge",
#'   ridge regression is used to fit a an outcome regression model via
#'   \code{glmnet::cv.glmnet} with default parameters. If "xgb", gradient
#'   boosting with a wide range of parameter values to cross-validate is used
#'   via \code{xgboost::xgb.cv} and the best parameters with respect to RMSE
#'   (for continuous outcomes) or misclassification rate (for binary/multi-class
#'   outcomes) are chosen. In both cases, the default \code{predict} method is
#'   used to generate in-sample predictions. If a function, denotes a
#'   user-supplied function that should be used for computing PE. Must take in a
#'   matrix of covariates as its first argument (model matrix?) and a vector
#'   outcome as its second argument. Must return a vector of in-sample
#'   predictions, which, if the outcome is binary or multi-class, must be
#'   maximum probability class labels.
#' @param user_PE_fit Deprecated; use argument `PE_method` instead. An optional
#'   function supplied by the user that can be used instead of those allowed for
#'   by \code{PE_method} to fit a model for the outcome from the covariates.
#'   Must take in a matrix of covariates as its first argument and a vector
#'   outcome as its second argument. Defaults to \code{NULL}.
#' @param user_PE_fit_params Deprecated; use argument `PE_method` instead. A
#'   named list of optional parameters to be used by \code{user_PE_fit}.
#'   Defaults to \code{NULL}.
#' @param user_PE_predict Deprecated; use argument `PE_method` instead. An
#'   optional function supplied by the user that can be used to generate
#'   predictions from the output of \code{user_PE_fit}. As its first argument,
#'   must take an object of the type returned by \code{user_PE_fit} and as its
#'   second, a matrix of values for which to generate predictions. When the
#'   outcome is binary or multi-class, must return the maximum probability class
#'   label. If not supplied, defaults to \code{predict}.
#' @param user_PE_predict_params Deprecated; use argument `PE_method` instead. A
#'   named list of optional parameters to be used by \code{user_PE_predict}.
#'   Defaults to \code{NULL}.
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
#' @param early_stop_pe Deprecated. A positive numeric. If FLAME
#'   attempts to drop a covariate that would lead to a PE above this value,
#'   FLAME stops. Defaults to \code{Inf}.
#' @param early_stop_bf Deprecated. A numeric value between 0 and 2. If FLAME
#'   attempts to drop a covariate that would lead to a BF below this value,
#'   FLAME stops. Defaults to 0.
#' @param missing_data Specifies how to handle missingness in \code{data}. If
#'   'none' (default), assumes no missing data; if 'drop', effectively drops
#'   units with missingness from the data and does not match them; if 'ignore',
#'   ignores matches on missing covariates; and if 'impute', imputes the missing
#'   data via \code{mice::mice}.
#' @param missing_holdout Specifies how to handle missingness in \code{holdout}.
#'   If 'none' (default), assumes no missing data; if 'drop', effectively drops
#'   units with missingness and does not use them to compute PE; and if 'impute',
#'   imputes the missing data via \code{mice::mice}. In this last case, the PE
#'   at an iteration will be given by the average PE across all imputations.
#' @param missing_holdout_imputations If \code{missing_holdout} = 'impute', performs
#'   this many imputations of the missing data in \code{holdout} via
#'   \code{mice::mice}. Defaults to 5.
#' @param missing_data_imputations Defunct. If \code{missing_data} = 'impute', one
#'   round of imputation will be performed on \code{data} via \code{mice::mice}.
#'   To view results for multiple imputations, please wrap calls to \code{FLAME}
#'   in a loop. This argument will be removed in a future release.
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
           weights = NULL,
           PE_method = 'ridge',
           user_PE_fit = NULL, user_PE_fit_params = NULL,
           user_PE_predict = NULL, user_PE_predict_params = NULL,
           replace = FALSE, verbose = 2, return_pe = FALSE, return_bf = FALSE,
           early_stop_iterations = Inf, early_stop_epsilon = 0.25,
           early_stop_control = 0, early_stop_treated = 0,
           early_stop_pe = Inf, early_stop_bf = 0,
           missing_data = c('none', 'drop', 'ignore', 'impute'),
           missing_holdout = c('none', 'drop', 'impute'),
           missing_data_imputations = 1, missing_holdout_imputations = 5,
           impute_with_treatment = TRUE, impute_with_outcome = FALSE) {

    missing_data <- match.arg(missing_data)
    missing_holdout <- match.arg(missing_holdout)

    # Check the !identical clause
    if (!is.null(user_PE_fit) | !is.null(user_PE_fit_params) |
        !identical(user_PE_predict, predict)|!is.null(user_PE_predict_params)) {
      warning('Arguments `user_PE_fit`, `user_PE_fit_params` ',
              '`user_PE_predict`, and `user_PE_predict_params` are ',
              'deprecated and will be removed in a future release. Please use ',
              '`PE_method` instead.', call. = FALSE)
    }

    if (missing_data_imputations != 1) {
      missing_data_imputations <- 1
      warning('Argument `missing_data_imputations` is defunct and will be ',
              'removed in a future release. One round of imputation will be ',
              'performed for `data`. To perform multiple rounds of imputation',
              ', please wrap calls to `FLAME` in a loop.', call. = FALSE)
    }

    if (early_stop_pe < Inf) {
      warning('Argument `early_stop_pe` is deprecated and will ',
              'be removed in a future release. Please use `early_stop_epsilon`',
              ' to early stop based off predictive error.', call. = FALSE)
    }

    if (early_stop_bf > 0) {
      warning('Argument `early_stop_bf` is deprecated and will ',
              'be removed in a future release. Please use `early_stop_epsilon`',
              ' to early stop based off predictive error.', call. = FALSE)
    }
    input_args <- as.list(environment())

    return(do.call(AME, c(list(algo = 'FLAME'), input_args)))
}

DAME <-
  function(data, holdout = 0.1, C = 0.1,
           treated_column_name = 'treated', outcome_column_name = 'outcome',
           weights = NULL,
           PE_method = 'ridge', n_flame_iters = 0,
           user_PE_fit = NULL, user_PE_fit_params = NULL,
           user_PE_predict = NULL, user_PE_predict_params = NULL,
           replace = FALSE, verbose = 2, return_pe = FALSE, return_bf = FALSE,
           early_stop_iterations = Inf, early_stop_epsilon = 0.25,
           early_stop_control = 0, early_stop_treated = 0,
           early_stop_pe = Inf, early_stop_bf = 0,
           missing_data = c('none', 'drop', 'ignore', 'impute'),
           missing_holdout = c('none', 'drop', 'impute'),
           missing_data_imputations = 1, missing_holdout_imputations = 5,
           impute_with_treatment = TRUE, impute_with_outcome = FALSE) {

    missing_data <- match.arg(missing_data)
    missing_holdout <- match.arg(missing_holdout)

    if (missing_data_imputations != 1) {
      missing_data_imputations <- 1
      warning('Argument `missing_data_imputations` is defunct and will be ',
              'removed in a future release. One round of imputation will be ',
              'performed for `data`. To perform multiple rounds of imputation',
              ', please wrap calls to `DAME` in a loop.', call. = FALSE)
    }

    if (early_stop_pe < Inf) {
      warning('Argument `early_stop_pe` is deprecated and will ',
              'be removed in a future release. Please use `early_stop_epsilon`',
              ' to early stop based off predictive error.', call. = FALSE)
    }

    if (early_stop_bf > 0) {
      warning('Argument `early_stop_bf` is deprecated and will ',
              'be removed in a future release. Please use `early_stop_epsilon`',
              ' to early stop based off predictive error.', call. = FALSE)
    }

    input_args <- as.list(environment())
    return(do.call(AME, c(list(algo = 'DAME'), input_args)))
  }

# Take out defaults?
AME <- function(algo, data, holdout = 0.1, C = 0.1,
            treated_column_name = 'treated', outcome_column_name = 'outcome',
            weights = NULL,
            PE_method = 'ridge', n_flame_iters = 0,
            user_PE_fit = NULL, user_PE_fit_params = NULL,
            user_PE_predict = NULL, user_PE_predict_params = NULL,
            replace = FALSE, verbose = 2, return_pe = FALSE, return_bf = FALSE,
            early_stop_iterations = Inf, early_stop_epsilon = 0.25,
            early_stop_control = 0, early_stop_treated = 0,
            early_stop_pe = Inf, early_stop_bf = 0,
            missing_data = 0, missing_holdout = 0,
            missing_data_imputations = 1, missing_holdout_imputations = 5,
            impute_with_treatment = TRUE, impute_with_outcome = FALSE) {

  early_stop_params <-
    list(iterations = early_stop_iterations,
         epsilon = early_stop_epsilon,
         control = early_stop_control,
         treated = early_stop_treated,
         PE = early_stop_pe,
         BF = early_stop_bf)

  out <-
    preprocess(data, holdout, C,
               treated_column_name, outcome_column_name, n_flame_iters,
               PE_method, user_PE_fit, user_PE_fit_params,
               user_PE_predict, user_PE_predict_params,
               replace, verbose, return_pe, return_bf,
               early_stop_params,
               missing_data, missing_holdout,
               missing_holdout_imputations,
               impute_with_outcome, impute_with_treatment)

  data <- out$data[[1]]
  holdout <- out$holdout
  covs <- out$covs
  mapping <- out$mapping
  orig_missing <- out$orig_missing
  cov_names <- out$cov_names

  n_covs <- ncol(data) - 4 - !is.null(data$outcome)

  # List of MGs, each entry contains the corresponding MG's entries
  MGs <- vector('list', nrow(data))

  # Try and make matches on all covariates
  matches_out <- update_matches(data, replace, c(), n_covs, MGs)
  data <- matches_out$data
  MGs <- matches_out$MGs
  ###### DO SMM FOR PE STOPPING IF WEIGHTS NOT NULL

  active_cov_sets <- as.list(1:n_covs)
  processed_cov_sets <- list()

  # Predictive error using all covariates. Used for stopping condition.
  if (is.null(weights)) {
    baseline_PE <- get_PE(c(), 1:n_covs, holdout,
                          PE_method, user_PE_fit, user_PE_fit_params,
                          user_PE_predict, user_PE_predict_params)
    early_stop_params$baseline_PE <- baseline_PE
  }

  AME_out <- run_AME(data, active_cov_sets, processed_cov_sets, early_stop_params,
                     verbose, C, algo, weights, MGs, replace, n_flame_iters,
                     return_pe, return_bf, n_covs,
                     holdout,
                     PE_method, user_PE_fit, user_PE_fit_params,
                     user_PE_predict, user_PE_predict_params)

  data <- AME_out$data
  MGs <- AME_out$MGs
  store_pe <- AME_out$store_pe
  store_bf <- AME_out$store_bf

  AME_out <-
    postprocess(data, MGs, n_covs, mapping, orig_missing, return_pe, return_bf,
                store_pe, store_bf)

  return(AME_out)

}
