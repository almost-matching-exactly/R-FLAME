#' Almost Matching Exactly (AME) Algorithms for Discrete, Observational Data
#' @param data Data to be matched. Either a data frame or a path to a .csv file
#'   to be read into a data frame. Treatment must be described by a logical or
#'   binary numeric column with name \code{treated_column_name}. If supplied,
#'   outcome must be described by a column with name \code{outcome_column_name}.
#'   The outcome will be treated as continuous if numeric with more than two
#'   values, as binary if a two-level factor or numeric with values 0 and 1
#'   exclusively, and as multi-class if a factor with more than two levels. If
#'   the outcome column is omitted, matching will be performed but treatment
#'   effect estimation will not be possible. All columns not containing outcome
#'   or treatment will be treated as covariates for matching. Covariates are
#'   assumed to be categorical and will be coerced to factors, though they may
#'   be passed as either factors or numeric; if the former, unused levels will
#'   automatically be dropped. If you wish to use continuous covariates for
#'   matching, they should be binned prior to matching.
#' @param holdout Holdout data to be used to compute predictive error, if
#'   \code{weights} is not supplied. If a numeric scalar between 0 and 1, that
#'   proportion of \code{data} will be made into a holdout set and only the
#'   \emph{remaining proportion} of \code{data} will be matched. Otherwise, a
#'   data frame or a path to a .csv file. The holdout data must contain an
#'   outcome column with name \code{outcome_column_name}; other restrictions on
#'   column types are as for \code{data}. Covariate columns must have the same
#'   column names and order as \code{data}. This data will \emph{not} be
#'   matched. Defaults to 0.1.
#' @param treated_column_name Name of the treatment column in \code{data} and
#'   \code{holdout}. Defaults to 'treated'.
#' @param outcome_column_name Name of the outcome column in \code{holdout} and
#'   also in \code{data}, if supplied in the latter. Defaults to 'outcome'.
#' @param weights A positive numeric vector representing covariate importances.
#'   Supplying this argument prevents PE from being computed as it determines
#'   dropping order by forcing covariate subsets with lower weights to be
#'   dropped first. The weight of a covariate subset is defined to be the sum of
#'   the weights of the constituent covariates. Ties are broken at random.
#' @param PE_method Denotes how predictive error (PE) is to be computed. Either
#'   a string -- one of "ridge" (default) or "xgb" -- or a function. If "ridge",
#'   ridge regression is used to fit a an outcome regression model via
#'   \code{glmnet::cv.glmnet} with default parameters. If "xgb", gradient
#'   boosting with a wide range of parameter values to cross-validate is used
#'   via \code{xgboost::xgb.cv} and the best parameters with respect to RMSE
#'   (for continuous outcomes) or misclassification rate (for binary/multi-class
#'   outcomes) are chosen. In both cases, the default \code{predict} method is
#'   used to generate in-sample predictions. If a function, denotes a
#'   user-supplied function that should be used for computing PE. This function
#'   must be passed a data frame of covariates as its first argument and a
#'   vector of outcome values as its second argument. It must return a vector of
#'   in-sample predictions, which, if the outcome is binary or multi-class, must
#'   be maximum probability class labels. See below for examples.
#' @param user_PE_fit Deprecated; use argument `PE_method` instead. An optional
#'   function supplied by the user that can be used instead of those allowed for
#'   by \code{PE_method} to fit a model for the outcome from the covariates.
#'   This function will be passed a data frame of covariates
#'   as its first argument and a vector of outcome values as its
#'   second argument. See below for examples. Defaults to \code{NULL}.
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
#'   matched multiple times, on different sets of covariates. In this case, the
#'   balancing factor for \code{FLAME} is computing by dividing by the total
#'   number of treatment (control) units, instead of the number of unmatched
#'   treatment (control) units. Defaults to \code{FALSE}.
#' @param estimate_CATEs A logical scalar. If \code{TRUE}, CATEs for each unit
#'   are estimated throughout the matching procedure, which will be much faster
#'   than computing them after a call to \code{FLAME} or \code{DAME} for very
#'   large inputs. Defaults to \code{FALSE}.
#' @param verbose Controls how FLAME displays progress while running. If 0, no
#'   output. If 1, only outputs the stopping condition. If 2, outputs the
#'   iteration and number of unmatched units every 5 iterations, and the
#'   stopping condition. If 3, outputs the iteration and number of unmatched
#'   units every iteration, and the stopping condition. Defaults to 2.
#' @param return_pe A logical scalar. If \code{TRUE}, the predictive error (PE)
#'   at each iteration will be returned. Defaults to \code{FALSE}.
#' @param return_bf A logical scalar. If \code{TRUE}, the balancing factor (BF)
#'   at each iteration will be returned. Defaults to \code{FALSE}.
#' @param early_stop_iterations A positive integer, denoting an upper bound
#'   on the number of matching rounds to be performed. If 1, one round of
#'   exact matching is performed before stopping. Defaults to \code{Inf}.
#' @param early_stop_epsilon A nonnegative numeric. If fixed covariate weights
#'   are passed via \code{weights}, then the algorithm will stop before matching
#'   on a covariate set whose error is above \code{early_stop_epsilon}, where in
#'   this case the error is defined as: \eqn{1 - weight(covariate set matched
#'   on) / weight(all covariates)}. Otherwise, if \code{weights} is \code{NULL},
#'   if FLAME or DAME attempts to drop a covariate set that would raise the PE
#'   above (1 + \code{early_stop_epsilon}) times the baseline PE (the PE before
#'   any covariates have been dropped), the algorithm will stop. Defaults to
#'   0.25.
#' @param early_stop_control,early_stop_treated If the proportion of control,
#'   treated units, respectively, that are unmatched falls below this value, the
#'   matching algorithm will stop. Default to 0.
#' @param early_stop_pe Deprecated. A positive numeric. If FLAME attempts to
#'   drop a covariate that would lead to a PE above this value, FLAME stops.
#'   Defaults to \code{Inf}.
#' @param early_stop_bf Deprecated. A numeric value between 0 and 2. If FLAME
#'   attempts to drop a covariate that would lead to a BF below this value,
#'   FLAME stops. Defaults to 0.
#' @param missing_data Specifies how to handle missingness in \code{data}. If
#' 'none' (default), assumes no missing data. If 'drop', effectively drops units
#' with missingness from the data and does not match them (they will still
#' appear in the matched dataset that is returned, however). If 'keep', keeps
#' the missing values in the data; in this case, a unit can only match on sets
#' containing covariates it is not missing. If 'impute', imputes the missing
#' data via \code{mice::mice}.
#' @param missing_holdout Specifies how to handle missingness in \code{holdout}.
#'   If 'none' (default), assumes no missing data; if 'drop', drops units with
#'   missingness and does not use them to compute PE; and if 'impute', imputes
#'   the missing data via \code{mice::mice}. In this last case, the PE at an
#'   iteration will be given by the average PE across all imputations.
#' @param missing_holdout_imputations If \code{missing_holdout} = 'impute',
#'   performs this many imputations of the missing data in \code{holdout} via
#'   \code{mice::mice}. Defaults to 5.
#' @param missing_data_imputations Defunct. If \code{missing_data} = 'impute',
#'   one round of imputation will be performed on \code{data} via
#'   \code{mice::mice}. To view results for multiple imputations, please wrap
#'   calls to \code{FLAME} or \code{DAME} in a loop. This argument will be
#'   removed in a future release.
#' @param impute_with_treatment,impute_with_outcome If \code{TRUE}, use
#'   treatment, outcome, respectively, to impute covariates when either
#'   \code{missing_data} or \code{missing_holdout} is equal to \code{'impute'}.
#'   Default to \code{TRUE}, \code{FALSE}, respectively.
#'
#'
#' @section Introduction: FLAME and DAME are matching algorithms for
#'   observational causal inference on data with discrete (categorical)
#'   covariates. They match units that share identical values of certain
#'   covariates, as follows. The algorithms first make any possible \emph{exact}
#'   matches; that is, they match units that share identical values of all
#'   covariates (this is possible because covariates are discrete). They then
#'   iteratively drop a set of covariates and make any possible matches on the
#'   remaining covariates, until stopping. For each unit, DAME solves an
#'   optimization problem that finds the highest quality set of covariates the
#'   unit can be matched to others on, where quality is determined by how well
#'   that set of covariates predicts the outcome. FLAME approximates the
#'   solution to the problem solved by DAME; at each step, it drops the
#'   covariate leading to the smallest drop in match quality \eqn{MQ}, defined
#'   as \eqn{MQ = C · BF - PE}. Here, \eqn{PE} denotes the predictive error,
#'   which measures how important the dropped covariate is for predicting the
#'   outcome. The balancing factor \eqn{BF} measures the number of matches
#'   formed by dropping that covariate. In this way, FLAME encourages matching
#'   on covariates more important to the outcome and also making many matches.
#'   The hyperparameter \eqn{C} controls the balance between these two
#'   objectives. In both cases, a machine learning algorithm trained on a
#'   holdout dataset is responsible for learning the quality / importance of
#'   covariates. For more details on the algorithms, please see the vignette,
#'   the FLAME paper \href{https://arxiv.org/pdf/1707.06315.pdf}{here} and/or
#'   the DAME paper \href{https://arxiv.org/pdf/1806.06802.pdf}{here}.
#'
#' @section Stopping Rules: By default, both \code{FLAME} and \code{DAME} stop
#'   when 1. all covariates have been dropped or 2. all treatment or control
#'   units have been matched. This behavior can be modified by the arguments
#'   whose prefix is "early_stop". With the exception of
#'   \code{early_stop_iterations}, all the rules come into play \emph{before}
#'   the offending covariate set is dropped. For example, if
#'   \code{early_stop_control = 0.2} and at the current iteration, dropping the
#'   covariate leading to highest match quality is associated with a unmatched
#'   control proportion of 0.1, FLAME will stop \emph{without} dropping this
#'   covariate.
#'
#' @section Missing Data: \code{FLAME} and \code{DAME} offer functionality for
#'   handling missing data in the covariates, for both the \code{data} and
#'   \code{holdout} sets. This functionality can be specified via the arguments
#'   whose prefix is "missing" or "impute". It allows for ignoring missing data,
#'   imputing it, or (for \code{data}) not matching on missing values. If
#'   \code{data} is imputed, imputation will be done once and the matching
#'   algorithm will be run on the imputed dataset. If \code{holdout} is imputed,
#'   the predictive error at an iteration will be the average of predictive
#'   errors across all imputed \code{holdout} datasets. Units with missingness
#'   in the treatment or outcome will be dropped.
#' @name AME
#' @return An object of type \code{ame}, which by default is a list of 4
#'   entries:
#' \describe{
#' \item{data}{The original data frame with several modifications:
#'   \enumerate{
#'     \item An extra logical column, \code{data$matched}, that indicates
#'     whether or not a unit was matched.
#'     \item An extra numeric column, \code{data$weight}, that denotes on how
#'     many different sets of covariates a unit was matched. This will only be
#'     greater than 1 when \code{replace = TRUE}.
#'     \item The columns denoting treatment and outcome will be moved after all
#'     covariate columns.
#'     \item If \code{replace} is \code{FALSE}, a column containing a matched
#'     group identifier for each unit.
#'     \item If, \code{estimate_CATEs = TRUE}, a column containing the CATE
#'     estimate for each unit.
#'    }
#'  }
#'  \item{MGs}{A list whose \eqn{i}'th entry contains the indices of units in
#'  the main matched group of the \eqn{i}'th unit.}
#'  \item{cov_sets}{A list whose \eqn{i}'th entry contains the covariates set
#'  \strong{not} matched on in the \eqn{i}'th iteration.}
#'  \item{info}{A list containing miscellaneous information about the data and
#'  matching specifications. Primarily for use by \code{*.ame} methods.}
#' }
#'
#' @examples
#' \dontrun{
#' data <- gen_data()
#' holdout <- gen_data()
#' # FLAME with replacement, stopping after dropping a single covariate
#' FLAME_out <- FLAME(data = data, holdout = holdout,
#'                    replace = TRUE, early_stop_iterations = 2)
#'
#' # Use a linear model to compute predictive error. Call DAME without
#' # replacement, returning predictive error at each iteration.
#' my_PE <- function(X, Y) {
#'   return(lm(Y ~ ., as.data.frame(cbind(X, Y = Y)))$fitted.values)
#' }
#' DAME_out <- DAME(data = data, holdout = holdout,
#'                  PE_method = my_PE, return_PE = TRUE)
#' }
#' @importFrom stats model.matrix predict rbinom rnorm var complete.cases median
#'   density
#' @importFrom utils flush.console read.csv write.csv combn
#' @importFrom graphics abline axis barplot image legend par title
NULL
#> NULL
#' @param C A finite, positive scalar denoting the tradeoff between BF and PE in
#'   the FLAME algorithm. Higher C prioritizes more matches and lower C
#'   prioritizes not dropping important covariates. Defaults to 0.1.
#' @rdname AME
#' @export
FLAME <-
  function(data, holdout = 0.1, C = 0.1,
           treated_column_name = 'treated', outcome_column_name = 'outcome',
           weights = NULL,
           PE_method = 'ridge',
           user_PE_fit = NULL, user_PE_fit_params = NULL,
           user_PE_predict = NULL, user_PE_predict_params = NULL,
           replace = FALSE, estimate_CATEs = FALSE,
           verbose = 2, return_pe = FALSE, return_bf = FALSE,
           early_stop_iterations = Inf, early_stop_epsilon = 0.25,
           early_stop_control = 0, early_stop_treated = 0,
           early_stop_pe = Inf, early_stop_bf = 0,
           missing_data = c('none', 'drop', 'keep', 'impute'),
           missing_holdout = c('none', 'drop', 'impute'),
           missing_data_imputations = 1, missing_holdout_imputations = 5,
           impute_with_treatment = TRUE, impute_with_outcome = FALSE) {

    missing_data <- match.arg(missing_data)
    missing_holdout <- match.arg(missing_holdout)

    # Check the !identical clause
    if (!is.null(user_PE_fit) | !is.null(user_PE_fit_params) |
        !is.null(user_PE_predict) | !is.null(user_PE_predict_params)) {
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

    return(do.call(AME, c(list(algo = 'FLAME',
                               n_flame_iters = Inf), input_args)))
}

#' @param n_flame_iters Specifies that this many iterations of FLAME should be
#'   run before switching to DAME. This can be used to speed up the matching
#'   procedure as FLAME rapidly eliminates irrelevant covariates, after which
#'   DAME will make higher quality matches on the remaining variables.
#' @rdname AME
#' @export
DAME <-
  function(data, holdout = 0.1,
           treated_column_name = 'treated', outcome_column_name = 'outcome',
           weights = NULL,
           PE_method = 'ridge', n_flame_iters = 0,
           user_PE_fit = NULL, user_PE_fit_params = NULL,
           user_PE_predict = NULL, user_PE_predict_params = NULL,
           replace = FALSE,  estimate_CATEs = FALSE, verbose = 2,
           return_pe = FALSE, return_bf = FALSE,
           early_stop_iterations = Inf, early_stop_epsilon = 0.25,
           early_stop_control = 0, early_stop_treated = 0,
           early_stop_pe = Inf, early_stop_bf = 0,
           missing_data = c('none', 'drop', 'keep', 'impute'),
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
    return(do.call(AME, c(list(algo = 'DAME', C = 0.1), input_args)))
  }

AME <- function(algo, data, holdout, C,
            treated_column_name, outcome_column_name,
            weights,
            PE_method, n_flame_iters,
            user_PE_fit, user_PE_fit_params,
            user_PE_predict, user_PE_predict_params,
            replace, estimate_CATEs, verbose, return_pe, return_bf,
            early_stop_iterations, early_stop_epsilon,
            early_stop_control, early_stop_treated,
            early_stop_pe, early_stop_bf,
            missing_data, missing_holdout,
            missing_data_imputations, missing_holdout_imputations,
            impute_with_treatment, impute_with_outcome) {

  early_stop_params <-
    list(iterations = early_stop_iterations,
         epsilon = early_stop_epsilon,
         control = early_stop_control,
         treated = early_stop_treated,
         PE = early_stop_pe,
         BF = early_stop_bf)

  out <-
    preprocess(data, holdout, C, algo, weights,
               treated_column_name, outcome_column_name, n_flame_iters,
               PE_method, user_PE_fit, user_PE_fit_params,
               user_PE_predict, user_PE_predict_params,
               replace, estimate_CATEs, verbose, return_pe, return_bf,
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
  info <- out$info

  n_covs <- sum(!(colnames(data) %in%
                    c('MG', 'missing', 'CATE', 'matched',
                      'weight', 'treated', 'outcome')))

  # List of MGs, each entry contains the corresponding MG's entries
  MGs <- vector('list', nrow(data))

  # Try and make matches on all covariates
  matches_out <- update_matches(data, replace, c(), n_covs, MGs, list(), info)
  data <- matches_out$data
  MGs <- matches_out$MGs

  active_cov_sets <- as.list(1:n_covs)
  processed_cov_sets <- list()

  # Predictive error using all covariates. Used for stopping condition.
  if (is.null(weights)) {
    baseline_PE <- get_PE(c(), 1:n_covs, holdout,
                          PE_method, user_PE_fit, user_PE_fit_params,
                          user_PE_predict, user_PE_predict_params)
    early_stop_params$baseline_PE <- baseline_PE
  }

  AME_out <- run_AME(data, active_cov_sets, processed_cov_sets,
                     early_stop_params, verbose, C, algo, weights, MGs, replace,
                     n_flame_iters, return_pe, return_bf, n_covs, holdout,
                     PE_method, user_PE_fit, user_PE_fit_params,
                     user_PE_predict, user_PE_predict_params, info)

  AME_out <- postprocess(AME_out, n_covs, mapping, orig_missing,
                         return_pe, return_bf, info)

  return(AME_out)

}
