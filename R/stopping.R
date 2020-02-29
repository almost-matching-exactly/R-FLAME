show_progress <- function(verbose, iter, data) {
  if (verbose == 1) {
    message(paste('Starting iteration', iter, 'of FLAME'))
  }
  else if (verbose == 2) {
    if (iter %% 5 == 0) {
      message(paste('Starting iteration', iter, 'of FLAME'))
      message(paste(sum(!data$matched), 'unmatched units remaining'))
    }
  }
  else if (verbose == 3) {
    message(paste('Starting iteration', iter, 'of FLAME'))
    message(paste(sum(!data$matched), 'unmatched units remaining'))
  }
}

early_stop_BF <- function(BF, early_stop_bf, prop_c_unmatched, prop_t_unmatched,
                          early_stop_un_c_frac, early_stop_un_t_frac) {

  if (BF < early_stop_bf) {
    message('FLAME stopping: balancing factor would have dropped below ', early_stop_bf)
    return(TRUE)
  }

  if (prop_c_unmatched < early_stop_un_c_frac) {
    message('FLAME stopping: proportion of control units ',
            'that are unmatched would have dropped below ', early_stop_un_c_frac)
    return(TRUE)
  }
  if (prop_t_unmatched < early_stop_un_t_frac) {
    message('FLAME stopping: proportion of treatment units ',
            'that are unmatched would have dropped below ', early_stop_un_t_frac)
    return(TRUE)
  }
  return(FALSE)
}

early_stop_PE <- function(PE, early_stop_pe, epsilon, baseline_PE) {
  if (PE > early_stop_pe) { # should be >
    message('FLAME stopping: predictive error would have risen above ', early_stop_pe)
    return(TRUE)
  }
  if (PE > (1 + epsilon) * baseline_PE) {
    message('FLAME stopping: predictive error would have risen ',
            100 * epsilon, '% above the baseline.')
    return(TRUE)
  }
  return(FALSE)
}

early_stop <- function(iter, data, covs, early_stop_iterations) {
  if (length(covs) == 1) {
    message('FLAME stopping: all covariates dropped')
    return(TRUE)
  }

  if (all(data$matched)) {
    message('FLAME stopping: all units matched')
    return(TRUE)
  }

  if (iter >= early_stop_iterations) {
    message('FLAME stopping: completed ', iter, ' iterations')
    return(TRUE)
  }

  if (sum(!data$matched & data$treated == 0) == 0) {
    message('FLAME stopping: all control units matched')
    return(TRUE)
  }

  if (sum(!data$matched & data$treated == 1) == 0) {
    message('FLAME stopping: all treatment units matched')
    return(TRUE)
  }

  return(FALSE)
}

# early_stop <- function(verbose, iter, covs, data, early_stop_iterations,
#                        epsilon, baseline_PE, store_pe, store_bf,
#                        early_stop_un_c_frac, early_stop_un_t_frac,
#                        early_stop_bf, early_stop_pe) {
#
#   type <- early_stop_type(iter, covs, data, early_stop_iterations,
#                           epsilon, baseline_PE, store_pe, store_bf,
#                           early_stop_un_c_frac, early_stop_un_t_frac,
#                           early_stop_bf, early_stop_pe)
#   if (verbose == 0) {
#     if (type == 0) {
#       return(FALSE)
#     }
#     return(TRUE)
#   }
#
#   if (type == 0) {
#     show_progress(verbose, iter, data)
#     return(FALSE)
#   }
#
#   if (type == 1) {
#     message('FLAME stopping: all covariates dropped')
#   }
#   if (type == 2) {
#     message('FLAME stopping: all units matched')
#   }
#   if (type == 3) {
#     message('FLAME stopping: completed ', iter, ' iterations')
#   }
#   if (type == 4) {
#     message('FLAME stopping: balancing factor dropped below ', early_stop_bf)
#   }
#   if (type == 5) {
#     message('FLAME stopping: predictive error rose above ', early_stop_pe)
#   }
#   if (type == 6) {
#     message('FLAME stopping: proportion of control units ',
#             'that are unmatched dropped below ', early_stop_un_c_frac)
#   }
#   if (type == 7) {
#     message('FLAME stopping: proportion of treatment units ',
#             'that are unmatched dropped below ', early_stop_un_t_frac)
#   }
#   if (type == 8) {
#     message('FLAME stopping: predictive error rose above ', 100 * epsilon, '% of the baseline.')
#   }
#   return(TRUE)
# }
#
# early_stop_type <-
#   function(iter, covs, data, early_stop_iterations,
#            epsilon, baseline_PE, store_pe, store_bf,
#            early_stop_un_c_frac, early_stop_un_t_frac,
#            early_stop_bf, early_stop_pe) {
#
#     if (length(covs) == 1) {
#       return(1)
#     }
#
#     if (all(data$matched)) {
#       return(2)
#     }
#
#     if (iter >= early_stop_iterations) {
#       return(3)
#     }
#
#     if (iter > 0) {
#       if (store_bf[iter] < early_stop_bf) {
#         return(4)
#       }
#
#       if (store_pe[iter] > early_stop_pe) { ###### PE is bad so should this be > than?
#         return(5)
#       }
#
#       if (store_pe[iter] > (1 + epsilon) * baseline_PE) {
#         return(8)
#       }
#     }
#
#     if (sum(data$treated == 0 & !data$matched) / sum(data$treated == 0) <
#         early_stop_un_c_frac) {
#       return(6)
#     }
#     if (sum(data$treated == 1 & !data$matched) / sum(data$treated == 1) <
#         early_stop_un_t_frac) {
#       return(7)
#     }
#     return(0)
#   }
