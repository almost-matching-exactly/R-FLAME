get_new_cov_set <- function(active_cov_sets, covs, weights, C, algo,
                            data, holdout,
                            PE_method, user_PE_fit, user_PE_fit_params,
                            user_PE_predict, user_PE_predict_params, replace) {

  if (!is.null(weights)) {
    max_weight <- 0
    p <- length(weights)
    J <- c(1:p)
    for (cov_set in active_cov_sets) {
      if (!is.null(weights)) {
        cov_set_weight <- sum(weights[-cov_set])
      }
      if (cov_set_weight > max_weight) {
        max_weight <- cov_set_weight
        best_cov_set <- cov_set
      }
    }
    return(list(cov_set = best_cov_set,
                PE = NULL,
                BF = NULL))
  }

  PE <- sapply(active_cov_sets, get_PE, covs, holdout,
               PE_method, user_PE_fit, user_PE_fit_params,
               user_PE_predict, user_PE_predict_params)

  if (algo == 'DAME') {
    return(list(cov_set = active_cov_sets[[which.min(PE)]],
                PE = min(PE),
                BF = NULL))
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
    lapply(active_cov_sets[drop_candidates], get_BF, data, replace, covs)
  BF <- sapply(BF_out, function(x) x[['BF']])

  MQ <- C * BF - PE

  # (First, in unlikely case of ties) covariate yielding highest MQ
  drop <- which.max(MQ)

  return(list(cov_set = active_cov_sets[[drop_candidates[drop]]],
              PE = PE[drop],
              BF = BF_out[[drop]]))
}

update_cov_sets <- function(active_cov_sets, processed_cov_sets, covs, weights, C, algo,
                            data, holdout,
                            PE_method, user_PE_fit, user_PE_fit_params,
                            user_PE_predict, user_PE_predict_params, replace) {

  tmp <- get_new_cov_set(active_cov_sets, covs, weights, C, algo,
                         data, holdout,
                         PE_method, user_PE_fit, user_PE_fit_params,
                         user_PE_predict, user_PE_predict_params, replace)
  curr_cov_set <- tmp$cov_set ## actually this is the thing that is just dropped
  PE <- tmp$PE
  BF <- tmp$BF
  if (algo == 'FLAME') {
    active_cov_sets <-
      lapply(setdiff(covs, curr_cov_set), function(x) c(x, curr_cov_set))
  }
  else if (algo == 'DAME') {
    Z_h <- GenerateNewActiveSets(curr_cov_set, processed_cov_sets)
    active_cov_sets <- remove_from_list(active_cov_sets, curr_cov_set)
    active_cov_sets <- append(active_cov_sets, Z_h)
  }
  processed_cov_sets <- append(processed_cov_sets, list(curr_cov_set))

  return(list(current = curr_cov_set,
              active = active_cov_sets,
              processed = processed_cov_sets,
              PE = PE,
              BF = BF))
}
