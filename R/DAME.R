make_MGs2 <- function(MGs, match_vals, newly_matched_units) {

  # Takes all the units that were matched on these p covariates and separates
  # them into matched groups based off their unique values of those covariates
  # Returns a list with two items
  # MGs: a vector, each entry of which is a list corresponding to a different MG
  #   and contains the indices of the corresponding units
  # matched_on: a vector the same length as MGs, each entry of which is a list
  #  detailing the covariates and their values that the units in the
  #  corresponding MG matched on
  n <- length(MGs)
  for (i in seq_along(newly_matched_units)) {
    MGs[[newly_matched_units[i]]] <-
      (1:n)[match_vals == match_vals[newly_matched_units[i]]]
  }
  return(MGs)

  unique_MGs <- unique(index)
  n_MGs <- length(unique_MGs)

  MGs <- vector('list', length = n_MGs)
  matched_on <- vector('list', length = n_MGs)

  for (i in 1:n_MGs) {
    members <- matched_units[which(index == unique_MGs[i])]
    MGs[[i]] <- members
    matched_on[[i]] <- data[members[1], covs, drop = FALSE]
    rownames(matched_on[[i]]) <- NULL
    names(matched_on[[i]]) <- cov_names[covs]
  }

  return(list(MGs = MGs,
              matched_on = matched_on))
}

DAME <- function(data, w, repeats = FALSE, num_iter) {
  h <- 1
  data$matched <- FALSE
  data_all <- data
  data$weight <- 0
  n <- nrow(data)
  MGs <- vector('list', length = n)
  # matched_on <- vector('list', length = num_iter)
  # data_matched <- vector('list', length = num_iter)
  p <- length(w)

  active_cov_sets <- as.list(1:p)

  processed_cov_sets <- list()
  J <- c(1:p)

  # First try grouping on all covs
  match <- GroupedMR(data_all, J, repeats)
  match_vals <- match[[1]]
  newly_matched <- match[[2]]

  data_all$matched[newly_matched] <- TRUE

  made_matches <- length(newly_matched) > 0

  if (made_matches) {
    # browser()
    MGs <- make_MGs2(MGs, match_vals, newly_matched)
    # new_MGs <- make_MGs(data_all, index, new_matched_indices, J, colnames(data[, J]))
    # MGs[[h]] <- new_MGs$MGs
    # matched_on[[h]] <- new_MGs$matched_on
    # data_matched[[h]] <- new_matched_indices
  }

  unmatched_indices <- which(!data_all$matched)
  h <- h + 1

  # Main loop
  while (length(unmatched_indices) > 0 & h <= num_iter) {

    # Find best covariate-set to drop from active covariate sets
    curr_cov_set <- decide_drop(active_cov_sets, w)

    # Check early stopping condition theta_s * w < 5% sum(w)
    ### EDIT
    if (sum(w) * 0.05 > as.integer(!(J %in% curr_cov_set)) %*% w) {
      break
    }
    cov <- setdiff(J, curr_cov_set)

    match <- GroupedMR(data_all, cov, repeats)
    match_vals <- match[[1]]
    newly_matched <- match[[2]]
    matched <- match[[3]]

    made_matches <- length(newly_matched) > 0

    if (made_matches) {
      browser()
      data_all$matched[newly_matched] <- TRUE
      MGs <- make_MGs2(MGs, match_vals, newly_matched)
      # new_MGs <- make_MGs(data_all, index, new_matched_indices, cov, colnames(data[,cov]))
      # MGs[[h]] <- c(MGs[[h - 1]], new_MGs$MGs)
      # matched_on[[h]] <- c(matched_on[[h - 1]],  new_MGs$matched_on)
      # data_matched[[h]] <- c(data_matched[[h - 1]], new_matched_indices)

      data_all[newly_matched, cov] <- '*'
      data_all$weight[matched] <- data_all$weight[matched] + 1

    }

    Z_h <- GenerateNewActiveSets(curr_cov_set, processed_cov_sets)
    active_cov_sets <- setdiff(active_cov_sets, curr_cov_set)
    active_cov_sets <- append(active_cov_sets, Z_h)
    processed_cov_sets <- append(processed_cov_sets, curr_cov_set)
    # unmatched_indices <- setdiff(unmatched_indices, new_matched_indices)
    h <- h + 1
  }
  return(list(MGs = MGs, data = data_all))
  return(list(data_matched = data_matched,
              MGs = MGs,
              matched_on = matched_on))
}
