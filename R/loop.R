run_AME <- function(data, active_cov_sets, processed_cov_sets, early_stop_params,
                    verbose, C, algo, weights, MGs, replace, n_flame_iters,
                    return_pe, return_bf, n_covs,
                    holdout,
                    PE_method, user_PE_fit, user_PE_fit_params,
                    user_PE_predict, user_PE_predict_params) {
  iter <- 0
  store_pe <- NULL
  store_bf <- NULL

  matches_out <- list(data = data, MGs = MGs)

  covs <- 1:n_covs
  if (algo == 'FLAME') {
    n_flame_iters <- Inf
  }
  while (!early_stop(iter, data, n_covs, active_cov_sets, early_stop_params, verbose, algo)) {
    if (iter < n_flame_iters) {
      algo <- 'FLAME'
    }
    else {
      algo <- 'DAME'
    }
    iter <- iter + 1
    show_progress(verbose, iter, data, algo)
    cov_sets <- update_cov_sets(active_cov_sets, processed_cov_sets, covs, weights, C, algo,
                                data, holdout,
                                PE_method, user_PE_fit, user_PE_fit_params,
                                user_PE_predict, user_PE_predict_params, replace)
    curr_cov_set <- cov_sets$current
    active_cov_sets <- cov_sets$active
    processed_cov_sets <- cov_sets$processed

    PE <- cov_sets$PE
    BF <- cov_sets$BF

    ## Note for Vittorio: min(PE)?
    if (early_stop_PE(PE, early_stop_params, verbose)) {
      if (return_pe) {
        matches_out$PE <- store_pe
      }
      if (return_bf) {
        matches_out$BF <- store_bf
      }
      return(matches_out)
    }

    if (!is.null(BF)) {
      prop_unmatched <- BF[['prop_unmatched']]
    }

    if (early_stop_BF(BF$BF,
                      prop_unmatched[['control']], prop_unmatched[['treated']],
                      early_stop_params, verbose)) {

      if (return_pe) {
        matches_out$PE <- store_pe
      }
      if (return_bf) {
        matches_out$BF <- store_bf
      }
      return(matches_out)
    }

    if (!is.null(weights)) {
      store_pe <- c(store_pe, PE)
      store_bf <- c(store_bf, BF$BF) # or something
    }

    # Make new matches having dropped a covariate
    ## Ideally should just return this from MQ so you don't have to redo it
    matches_out <- update_matches(data, replace, curr_cov_set, n_covs, MGs)
    data <- matches_out$data
    MGs <- matches_out$MGs
  }
  if (return_pe) {
    matches_out$PE <- store_pe
  }
  if (return_bf) {
    matches_out$BF <- store_bf
  }
  return(matches_out)
}
