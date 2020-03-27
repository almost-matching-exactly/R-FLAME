show_progress <- function(verbose, iter, data) {
  n <- nrow(data)
  n_digits <- nchar(n)
  padding <- rep(' ', n_digits)
  if (verbose == 2) {
    if (iter %% 5 == 0) {
      message(paste0('Starting iteration ', iter, ' of FLAME (',
                     sum(!data$matched), ' unmatched units remaining)',
                     padding, '\r'),
              appendLF = FALSE)
      flush.console()
    }
  }
  else if (verbose == 3) {
    message(paste0('Starting iteration ', iter, ' of FLAME (',
                  sum(!data$matched), ' unmatched units remaining)',
                  padding, '\r'),
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

pretty_print <- function(to_print, iter, verbose) {
  if (verbose == 0) {
    return()
  }
  if ((verbose == 2 & (iter %% 5 == 0)) | verbose == 3) {
    message('\n', appendLF = FALSE)
  }
  message(to_print)
}

early_stop <- function(iter, data, covs, early_stop_iterations, verbose) {
  if (length(covs) == 1) {
    pretty_print('FLAME stopping: only one covariate remaining',
                 iter, verbose)
    return(TRUE)
  }

  if (all(data$matched)) {
    pretty_print('FLAME stopping: all units matched',
                 iter, verbose)
    return(TRUE)
  }

  if (iter >= early_stop_iterations) {
    pretty_print(paste0('FLAME stopping: completed ', iter, ' iterations'),
                 iter, verbose)
    return(TRUE)
  }

  if (sum(!data$matched & data$treated == 0) == 0) {
    pretty_print('FLAME stopping: all control units matched',
                 iter, verbose)
    return(TRUE)
  }

  if (sum(!data$matched & data$treated == 1) == 0) {
    pretty_print('FLAME stopping: all treatment units matched',
                 iter, verbose)
    return(TRUE)
  }

  return(FALSE)
}
