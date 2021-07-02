postprocess <- function(AME_out, n_covs, mapping, orig_missing,
                        return_pe, return_bf, info) {

  data <- AME_out$data
  MGs <- AME_out$MGs
  PE <- AME_out$PE
  BF <- AME_out$BF
  cov_sets <- AME_out$cov_sets

  algo <- info$algo
  replace <- info$replacement
  treated_column_name <- info$treatment
  outcome_column_name <- info$outcome
  missing_data <- info$missing_data
  missing_holdout <- info$missing_holdout

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

  cov_inds <- which(!(colnames(data) %in%
                        c('treated', 'outcome', 'matched',
                          'weight', 'MG', 'CATE')))

  # Already done in factor_remap?
  data[, cov_inds] <- lapply(data[, cov_inds, drop = FALSE], droplevels)

  # Replace original variable names
  colnames(data)[colnames(data) == 'outcome'] <- outcome_column_name
  colnames(data)[colnames(data) == 'treated'] <- treated_column_name

  ret_list <- list(data = data, MGs = MGs, cov_sets = cov_sets, info = info)

  if (return_pe) {
    ret_list <- c(ret_list, 'PE' = list(PE))
  }
  if (return_bf) {
    ret_list <- c(ret_list, 'BF' = list(BF))
  }

  if (!info$estimate_CATEs || info$outcome_type != 'continuous') {
    ret_list$data$CATE <- NULL
  }

  class(ret_list) <- 'ame'

  return(ret_list)
}

