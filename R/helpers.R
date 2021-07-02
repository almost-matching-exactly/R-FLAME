gen_new_active_sets <- function(s, delta) {
  k <- length(s)
  Z <- list()

  # 3
  if (length(delta) == 0) {
    return (Z)
  }

  delta_k <- list()
  counter <- 1
  for (i in seq_along(delta)) {
    if (length(delta[[i]]) == k) {
      delta_k[[counter]] <- delta[[i]]
      counter <- counter + 1
    }
  }
  delta_k[[counter]] <- s

  # 4, 5
  supp <- table(unlist(delta_k))

  # 6
  omega <- setdiff(strtoi(names(supp[supp >= k])), s)

  # 7
  counter <- 1
  if (all(supp[match(s, names(supp))] >= k)) {
    # 8
    for (a in omega) {
      # 9
      # Necessary to sort?
      r <- sort(c(s, a))
      # 10
      if (all(combn(r, k, simplify = FALSE) %in% delta_k)) {
        # 11
        Z[[counter]] <- r
        counter <- counter + 1
      }
    }
  }
  return(Z)
}


warn_level_drop <- function(data, inds, type) {
  # Used in `preprocess.R` to drop extraneous levels and warn the user
  dropped_levels <- list()
  for (i in inds) {
    level_count <- table(data[, i])
    if (any(level_count == 0)) { # this may happen bc we split
      cov_name <- colnames(data)[i]
      dropped_levels <-
        c(dropped_levels, list(names(level_count)[which(level_count == 0)]))

      names(dropped_levels)[length(dropped_levels)] <- cov_name

      data[, i] <- droplevels(data[, i])
    }
  }

  if (length(dropped_levels) > 0) {
    warning_str <-
      paste0('In `', type, '`: dropping levels ',
            paste(vapply(seq_along(dropped_levels), function(i) {
              paste0('{', paste(dropped_levels[[i]], collapse = ', '),
                     '} from covariate `', names(dropped_levels)[i], '`')},
              character(1)), collapse = ' and '))

    warning(warning_str, call. = FALSE)
  }

  return(data)
}

remove_from_list <- function(lst, elem) {
  # Used in `update_cov_sets`
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
  # For computing b_u, c_u, b_u+, c_u+ values for bit-matching
  vals <- as.character(vals)
  tab <- table(vals)
  name <- names(tab)
  return(as.vector(tab[match(vals, name)]))
}

exact_match_bit <- function(data, covs, replace) {
  # Performs exact matching via bit-vectors of the units in `data` on the
  #   covariates specificed by `covs`.

  # gmp::as.bigz is for handling lots and lots of covariates so we don't
  #   have trouble with overflow issues

  if (!replace) {
    valid_matches <- which(!data$matched & !data$missing)
  }
  else {
    valid_matches <- which(!data$missing)
  }
  data <- data[valid_matches, ]

  n_levels <- vapply(data[, covs, drop = FALSE], nlevels, numeric(1))
  data_wo_t <- gmp::as.bigz(as.matrix(data[, covs[order(n_levels)]]))
  n_levels <- sort(n_levels)

  multiplier <- gmp::pow.bigz(n_levels, seq_along(n_levels) - 1)

  b_u <-
    as.vector(gmp::`%*%`(data_wo_t, multiplier))

  multiplier <- gmp::pow.bigz(n_levels, seq_along(n_levels))

  b_u_plus <-
    as.vector(gmp::add.bigz(gmp::`%*%`(data_wo_t, multiplier), data$treated))

  c_u <- aggregate_table(b_u)

  c_u_plus <- aggregate_table(b_u_plus)

  matched_on_covs <- (c_u != c_u_plus) & (c_u >= 2)

  # Those units matched for the FIRST time on this cov_set

  # If replace:
  #  (1:n)[which units matched on this cov_set AND were not previously matched]
  # If not replace:
  #  we threw out data$matched so !data$matched is who matched on this cov_set
  #  out of those not previously matched
  newly_matched <- valid_matches[matched_on_covs & !data$matched]
  # Those units matched on this cov_set (same as newly_matched if !replace)
  matched <- valid_matches[matched_on_covs]

  return(list(valid_matches = valid_matches,
              match_vals = b_u,
              newly_matched = newly_matched,
              matched = matched))
}

make_MGs <- function(MGs, valid_matches, match_vals,
                     matched, newly_matched, data, info) {

  if (info$estimate_CATEs && info$outcome_type == 'continuous') {
    Tr <- data$treated
    Y <- data$outcome
  }

  # b_u values for those first matched on this cov set
  newly_matched_vals <- match_vals[match(newly_matched, valid_matches)]

  MG_ids <- match(as.character(newly_matched_vals),
                  as.character(unique(newly_matched_vals)))

  MG_counter <- max(data$MG)

  for (i in seq_along(newly_matched)) {
    new_MG <- valid_matches[match_vals == newly_matched_vals[i]]
    MGs[[newly_matched[i]]] <- new_MG

    if (info$estimate_CATEs && info$outcome_type == 'continuous') {
      if (Tr[newly_matched[i]] == 1) {
        data$CATE[newly_matched[i]] <-
          Y[newly_matched[i]] - mean(Y[new_MG[Tr[new_MG] == 0]])
      }
      else {
        data$CATE[newly_matched[i]] <-
          mean(Y[new_MG[Tr[new_MG] == 1]]) - Y[newly_matched[i]]
      }
    }

    data$MG[new_MG] <- MG_counter + MG_ids[i]
  }
  return(list(MGs, data))
}

process_matches <- function(data, replace, covs, MGs, info) {
  match_out <- exact_match_bit(data, covs, replace)
  valid_matches <- match_out$valid_matches
  match_vals <- match_out$match_vals
  newly_matched <- match_out$newly_matched
  matched <- match_out$matched

  made_new_matches <- length(newly_matched) > 0

  if (made_new_matches) {
    MG_out <- make_MGs(MGs, valid_matches, match_vals,
                       matched, newly_matched, data, info)
    MGs <- MG_out[[1]]
    data <- MG_out[[2]]
  }
  return(list(MGs = MGs,
              newly_matched = newly_matched,
              matched = matched,
              data = data))
}

update_matches <- function(data, replace, dropped_cov_set,
                           n_covs, MGs, cov_sets, info) {
  processed_matches <-
    process_matches(data, replace, setdiff(1:n_covs, dropped_cov_set),
                    MGs, info)

  MGs <- processed_matches[[1]]
  newly_matched <- processed_matches[[2]]
  matched <- processed_matches[[3]]
  data <- processed_matches[[4]]

  if (length(newly_matched) > 0) {
    data$matched[newly_matched] <- TRUE
    data$weight[matched] <- data$weight[matched] + 1
  }

  cov_sets <- c(cov_sets, list(colnames(data)[dropped_cov_set]))
  return(list(data = data, MGs = MGs, cov_sets = cov_sets))
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
    which(!(seq_len(ncol(df[[1]])) %in% c(treatment_col_ind, outcome_col_ind)))

  # For all imputed data sets
  for (i in 1:n_df) {
    tmp_df <- df[[i]]

    if (type == 'holdout' | (type == 'data' & outcome_in_data)) {
      tmp_df <-
        cbind(tmp_df[, covariates, drop = FALSE], outcome_col, treatment_col)
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

    # covs denotes the covariates currently being matched on
    covs <- 1:n_covs

    if (type == 'data') {
      tmp_df$matched <- rep(FALSE, n)
      tmp_df$weight <- rep(0, n)
      tmp_df$missing <- is_missing
      tmp_df$MG <- rep(0, n)
      tmp_df$CATE <- rep(NA, n)
    }

    df[[i]] <- tmp_df
  }
  return(list(df = df,
              covs = covs,
              cov_names = cov_names))
}
