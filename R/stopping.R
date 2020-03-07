show_progress <- function(verbose, iter, data) {
  if (verbose == 1) {
    message(paste('Starting iteration', iter, 'of FLAME\r'), appendLF = FALSE)
    flush.console()
  }
  else if (verbose == 2) {
    if (iter %% 5 == 0) {
      message(paste0('Starting iteration ', iter, ' of FLAME (',
                     sum(!data$matched), ' unmatched units remaining)\r'),
              appendLF = FALSE)
      flush.console()
    }
  }
  else if (verbose == 3) {
    message(paste0('Starting iteration ', iter, ' of FLAME (',
                  sum(!data$matched), ' unmatched units remaining)\r'),
            appendLF = FALSE)
    flush.console()
  }
}

early_stop_BF <- function(BF, early_stop_bf, prop_c_unmatched, prop_t_unmatched,
                          early_stop_un_c_frac, early_stop_un_t_frac,
                          verbose) {

  if (BF < early_stop_bf) {
    if (verbose != 0) {
      message('FLAME stopping: balancing factor would have dropped below ', early_stop_bf)
    }
    return(TRUE)
  }

  if (prop_c_unmatched < early_stop_un_c_frac) {
    if (verbose != 0) {
      message('FLAME stopping: proportion of control units ',
              'that are unmatched would have dropped below ', early_stop_un_c_frac)
    }
    return(TRUE)
  }
  if (prop_t_unmatched < early_stop_un_t_frac) {
    if (verbose != 0) {
      message('FLAME stopping: proportion of treatment units ',
              'that are unmatched would have dropped below ', early_stop_un_t_frac)
    }
    return(TRUE)
  }
  return(FALSE)
}

early_stop_PE <- function(PE, early_stop_pe, early_stop_epsilon, baseline_PE, verbose) {
  if (PE > early_stop_pe) { # should be >
    if (verbose != 0) {
      message('FLAME stopping: predictive error would have risen above ', early_stop_pe)
    }
    return(TRUE)
  }
  if (PE > (1 + early_stop_epsilon) * baseline_PE) {
    if (verbose != 0) {
      message('FLAME stopping: predictive error would have risen ',
              100 * early_stop_epsilon, '% above the baseline.')
    }
    return(TRUE)
  }
  return(FALSE)
}

early_stop <- function(iter, data, covs, early_stop_iterations, verbose) {
  if (length(covs) == 1) {
    if (verbose != 0) {
      message('FLAME stopping: all covariates dropped')
    }
    return(TRUE)
  }

  if (all(data$matched)) {
    if (verbose != 0) {
      message('FLAME stopping: all units matched')
    }
    return(TRUE)
  }

  if (iter >= early_stop_iterations) {
    if (verbose != 0) {
      message('FLAME stopping: completed ', iter, ' iterations')
    }
    return(TRUE)
  }

  if (sum(!data$matched & data$treated == 0) == 0) {
    if (verbose != 0) {
      message('FLAME stopping: all control units matched')
    }
    return(TRUE)
  }

  if (sum(!data$matched & data$treated == 1) == 0) {
    if (verbose != 0) {
      message('FLAME stopping: all treatment units matched')
    }
    return(TRUE)
  }

  return(FALSE)
}
