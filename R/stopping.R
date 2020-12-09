show_progress <- function(verbose, iter, data, algo) {
  n <- nrow(data)
  n_digits <- nchar(n)
  padding <- rep(' ', n_digits)
  if (verbose == 2) {
    if (iter %% 5 == 0) {
      message(paste0('Starting iteration ', iter, ' of ', algo, ' (',
                     sum(!data$matched), ' unmatched units remaining)',
                     padding, '\r'),
              appendLF = FALSE)
      flush.console()
    }
  }
  else if (verbose == 3) {
    message(paste0('Starting iteration ', iter, ' of ', algo, ' (',
                  sum(!data$matched), ' unmatched units remaining)',
                  padding, '\r'),
            appendLF = FALSE)
    flush.console()
  }
}

early_stop_BF <-
  function(BF, prop_c_unmatched, prop_t_unmatched, early_stop_params, verbose) {
  if (is.null(BF)) { # DAME or fixed weights; no BF exists
    return(FALSE)
  }
  if (BF < early_stop_params$BF) {
    if (verbose != 0) {
      message('FLAME stopping: balancing factor would have dropped below ',
              early_stop_params$BF)
    }
    return(TRUE)
  }

  if (prop_c_unmatched < early_stop_params$control) {
    if (verbose != 0) {
      message('FLAME stopping: proportion of control units ',
              'that are unmatched would have dropped below ',
              early_stop_params$control)
    }
    return(TRUE)
  }
  if (prop_t_unmatched < early_stop_params$treated) {
    if (verbose != 0) {
      message('FLAME stopping: proportion of treatment units ',
              'that are unmatched would have dropped below ',
              early_stop_params$treated)
    }
    return(TRUE)
  }
  return(FALSE)
}

early_stop_PE <- function(PE, early_stop_params, verbose) {
  if (is.null(PE)) { # Using fixed weights
    return(FALSE)
  }
  if (PE > early_stop_params$PE) {
    if (verbose != 0) {
      message('FLAME stopping: predictive error would have risen above ',
              early_stop_params$PE)
    }
    return(TRUE)
  }
  if (PE > (1 + early_stop_params$epsilon) * early_stop_params$baseline_PE) {
    if (verbose != 0) {
      message('FLAME stopping: predictive error would have risen ',
              100 * early_stop_params$epsilon, '% above the baseline.')
    }
    return(TRUE)
  }
  return(FALSE)
}

pretty_print <- function(to_print, iter, verbose) {
  if (verbose == 0) {
    return()
  }
  if ((verbose == 2 & (iter %% 5 == 0)) | verbose == 3) {
    message('\n', appendLF = FALSE)
  }
  message(to_print)
}

early_stop <- function(iter, data, n_covs, active_cov_sets, early_stop_params, verbose, algo) {
  # Check this
  if (all(sapply(active_cov_sets, length) == n_covs)) {
    pretty_print(paste(algo, 'stopping: only one covariate remaining'),
                 iter, verbose)
    return(TRUE)
  }

  if (all(data$matched)) {
    pretty_print(paste(algo, 'stopping: all units matched'),
                 iter, verbose)
    return(TRUE)
  }

  if (iter >= early_stop_params$iterations) {
    pretty_print(paste0(algo, ' stopping: completed ', iter, ' iterations'),
                 iter, verbose)
    return(TRUE)
  }

  if (sum(!data$matched & data$treated == 0) == 0) {
    pretty_print(paste(algo, 'stopping: all control units matched'),
                 iter, verbose)
    return(TRUE)
  }

  if (sum(!data$matched & data$treated == 1) == 0) {
    pretty_print(paste(algo, 'stopping: all treatment units matched'),
                 iter, verbose)
    return(TRUE)
  }

  return(FALSE)
}
