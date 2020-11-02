remove_from_list <- function(lst, elem) {
  # Assumes vectors, no duplicates
  for (i in seq_along(lst)) {
    if (length(lst[[i]]) != length(elem)) {
      next
    }
    if (all(lst[[i]] == elem)) {
      lst[[i]] <- NULL
      return(lst)
    }
  }
}

aggregate_table <- function(vals) {
  vals <- as.character(vals)
  tab <- table(vals)
  name <- names(tab)
  return(as.vector(tab[match(vals, name)]))
}

# exact_match_bit takes a data frame, a set of covariates to match on, the
# treatment indicator column and the matched indicator column. it returns the
# array indicating whether each unit is matched (the first return value), and a
# list of indices for the matched units (the second return value)

exact_match_bit <- function(data, covs, replace) {
  # gmp::as.bigz is for handling lots and lots of covariates so we don't
  # have trouble with overflow issues

  if (!replace) {
    valid_matches <- which(!data$matched & !data$missing)
    data <- data[!data$matched, ]
  }
  else {
    valid_matches <- which(!data$missing)
  }

  data <- data[!data$missing, ]

  n_levels <- sapply(data[, covs, drop = FALSE], nlevels)
  data_wo_t <- gmp::as.bigz(as.matrix(data[, covs[order(n_levels)]]))
  n_levels <- sort(n_levels)

  # Compute b_u
  multiplier <- gmp::pow.bigz(n_levels, seq_along(n_levels) - 1)

  b_u <-
    as.vector(gmp::`%*%`(data_wo_t, multiplier))

  # Compute b_u+
  b_u_plus <-
    as.vector(gmp::add.bigz(data$treated, gmp::mul.bigz(b_u, gmp::as.bigz(n_levels))))

  # Compute c_u
  c_u <- aggregate_table(b_u)

  # Compute c_u+
  c_u_plus <- aggregate_table(b_u_plus)

  matched_on_covs <- (c_u != c_u_plus) & (c_u >= 2)

  # Those units matched for the FIRST time on this cov_set

  # If replace:
  #  (1:n)[which units matched on this cov_set AND were not previously matched]
  # If not replace:
  #  we threw out data$matched and so !data$matched is who matched on this cov_set
  #  out of those not previously matched
  newly_matched <- valid_matches[matched_on_covs & !data$matched]

  # Those units matched on this cov_set (same as newly_matched if !replace)
  matched <- valid_matches[matched_on_covs]
  return(list(match_vals = b_u[matched_on_covs],
              newly_matched = newly_matched,
              matched = matched))
}

make_MGs <- function(MGs, match_vals, matched, newly_matched) {
  n <- length(MGs)
  for (i in seq_along(newly_matched)) {
    MGs[[newly_matched[i]]] <- matched[match_vals == match_vals[i]]
  }
  return(MGs)
}

process_matches <- function(data, replace, covs, MGs) {
  match_out <- exact_match_bit(data[!data$missing, ], covs, replace)
  match_vals <- match_out$match_vals
  newly_matched <- match_out$newly_matched
  matched <- match_out$matched

  made_new_matches <- length(newly_matched) > 0

  if (made_new_matches) {
    MGs <- make_MGs(MGs, match_vals, matched, newly_matched)
  }
  return(list(MGs = MGs,
              newly_matched = newly_matched,
              matched = matched))
}

update_matches <- function(data, replace, dropped_cov_set, n_covs, MGs) {
  processed_matches <- process_matches(data, replace, setdiff(1:n_covs, dropped_cov_set), MGs)

  MGs <- processed_matches[[1]]
  newly_matched <- processed_matches[[2]]
  matched <- processed_matches[[3]]

  if (length(newly_matched) > 0) {
    data[newly_matched, dropped_cov_set] <- '*'
    data$matched[newly_matched] <- TRUE
    data$weight[matched] <- data$weight[matched] + 1
  }
  return(list(data = data, MGs = MGs))
}

# NULL default for holdout set
sort_cols <-
  function(df, treated_column_name, outcome_column_name,
           type, is_missing = NULL) {

  outcome_in_data <- !is.null(df[[1]][[outcome_column_name]])

  n_covs <- ncol(df[[1]]) - 1 - outcome_in_data # Ignore treatment, outcome
  n_df <- length(df) # Always pass in a list of data frames

  # Treatment and outcome will be constant across imputations
  treatment_col <- df[[1]][treated_column_name]

  if (type == 'holdout' | (type == 'data' & outcome_in_data)) {
    outcome_col <- df[[1]][outcome_column_name]
  }

  n <- nrow(df[[1]])

  # easier way to do this
  treatment_col_ind <- which(colnames(df[[1]]) == treated_column_name)
  outcome_col_ind <- which(colnames(df[[1]]) == outcome_column_name)
  covariates <-
    which(!(1:ncol(df[[1]]) %in% c(treatment_col_ind, outcome_col_ind)))

  # For all imputed data sets
  for (i in 1:n_df) {
    tmp_df <- df[[i]]

    if (type == 'holdout' | (type == 'data' & outcome_in_data)) {
      tmp_df <- cbind(tmp_df[, covariates, drop = FALSE], outcome_col, treatment_col)
    }
    else {
      tmp_df <- cbind(tmp_df[, covariates, drop = FALSE], treatment_col)
    }

    cov_names <- colnames(tmp_df)[1:n_covs]

    # Sorting data column names
    if (type == 'holdout' | (type == 'data' & outcome_in_data)) {
      colnames(tmp_df) <- c(cov_names, 'outcome', 'treated')
    }
    else {
      colnames(tmp_df) <- c(cov_names, 'treated')
    }

    # Don't I do this twice?
    if (type == 'data') {
      for (j in 1:n_covs) {
        levels(tmp_df[, j]) <- c(levels(tmp_df[, j]), '*')
      }
    }

    # covs denotes the covariates currently being matched on
    covs <- 1:n_covs

    if (type == 'data') {
      tmp_df$matched <- rep(FALSE, n)
      tmp_df$weight <- rep(0, n)
      tmp_df$missing <- is_missing
    }

    df[[i]] <- tmp_df
  }
  return(list(df = df,
              covs = covs,
              cov_names = cov_names))
}
