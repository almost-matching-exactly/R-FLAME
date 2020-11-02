postprocess <- function(data, MGs, n_covs, mapping, orig_missing, return_pe, return_bf,
                        store_pe, store_bf) {

  # Substitute covariate values of all unmatched units with the unmatched
  # covariate symbol '*'
  data[!data$matched, 1:n_covs] <- '*'

  data[, ncol(data)] <- NULL

  ## Swap back the original factor levels
  rev_mapping <- lapply(mapping, function(x) {
    tmp <- c(names(x), '*')
    names(tmp) <- c(x, '*')
    return(tmp)
  })

  data <-
    factor_remap(data, mapping = rev_mapping)$df

  # Undoes rownames conflicting with functions in post_matching.R
  # when holdout is taken from data
  if (!all(rownames(data) == 1:nrow(data))) {
    data[['original_ind']] <- rownames(data)
    rownames(data) <- 1:nrow(data)
  }

  missing_levels <-
    vapply(unique(data[orig_missing]),
           function(x) paste(x, '(m)'),
           character(1),
           USE.NAMES = FALSE)

  cov_inds <- which(!(colnames(data) %in%
                        c('treated', 'outcome', 'matched', 'weight', 'original_ind')))

  data[, cov_inds] <-
    lapply(data[, cov_inds, drop = FALSE], function(x) {
      levels(x) <- c(levels(x), missing_levels)
      return(x)})

  if (nrow(orig_missing) > 0) {
    for (i in 1:nrow(orig_missing)) {
      data[orig_missing[i, , drop = FALSE]] <-
        paste(data[orig_missing[i, , drop = FALSE]], '(m)')
    }
  }

  data[, cov_inds] <- lapply(data[, cov_inds, drop = FALSE], droplevels)

  ret_list <-
    list(data = data,
         MGs = MGs)

  # Implicitly checks that outcome exists in data;
  #  otherwise data$outcome returns NULL
  outcome_in_data <- is.numeric(data$outcome)
  if (outcome_in_data) {
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
