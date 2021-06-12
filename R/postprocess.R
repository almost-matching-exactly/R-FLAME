postprocess <- function(data, MGs, n_covs, mapping, orig_missing, return_pe, return_bf,
                        store_pe, store_bf, cov_sets) {

  data[['missing']] <- NULL
  data[['MG']][data[['MG']] == 0] <- NA

  ## Swap back the original factor levels
  rev_mapping <- lapply(mapping, function(x) {
    tmp <- names(x)
    names(tmp) <- x
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

  cov_inds <-
    which(!(colnames(data) %in%
              c('treated', 'outcome', 'matched', 'weight', 'original_ind', 'MG', 'missing')))

  data[, cov_inds] <- lapply(data[, cov_inds, drop = FALSE], droplevels)

  ret_list <-
    list(data = data,
         MGs = MGs,
         cov_sets = cov_sets)

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
