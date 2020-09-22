# A maybe not great way of factoring out a lot of common code involving
# iterating through data frames and units for computing MGs and CATEs
# Don't care about for loops because should be super fast unless maybe
# you're asking for all the MGs on a million-unit dataset.

# In any case, I wouldn't worry about going into this code.
iterate <- function(fun, units, FLAME_out, multiple) {
  if (is.null(names(FLAME_out))) { # Is a list of data frames
    n_df <- length(FLAME_out)
  }
  else { # Is a single data frame
    FLAME_out <- list(FLAME_out)
    n_df <- 1
  }

  if (!all(units %in% as.numeric(rownames(FLAME_out[[1]]$data)))) {
    stop('Supplied a unit not originally passed to `FLAME.`')
  }

  out <- vector('list', length = n_df)
  for (k in 1:n_df) {
    MGs <- FLAME_out[[k]]$MGs # MGs made by FLAME

    out[[k]] <- vector('list', length = length(units))

    for (i in seq_along(units)) {
      unit <- units[i]

      # Number of MGs to return for unit
      # Shouldn't need the ifelse
      # Only > 1 if multiple = TRUE and matched multiple times
      n_unit_MGs <- ifelse(multiple, FLAME_out[[k]]$data$weight[unit], 1)

      # To store MGs of unit
      out[[k]][[i]] <- vector('list', length = n_unit_MGs)
      counter <- 1
      for (j in 1:length(MGs)) {
        if (unit %in% MGs[[j]]) {
          out[[k]][[i]][[counter]] <- fun(FLAME_out[[k]], MGs[[j]], j)
          counter <- counter + 1
          if (counter == n_unit_MGs + 1) {
            break
          }
        }
      }
    }
    if (!multiple) {
      out[[k]] <- lapply(out[[k]], `[[`, 1)
    }
  }
  if (length(units) == 1) {
    out <- lapply(out, `[[`, 1)
  }
  if (n_df == 1) {
    out <- out[[1]]
  }
  return(out)
}

CATE_internal <- function(FLAME_out, MG, which_MG = NULL) {
  if (!('`outcome`' %in% colnames(FLAME_out$data))) {
    stop('Outcome not supplied in original call to `FLAME`; ',
         'cannot compute CATE')
  }
  if (is.factor(FLAME_out$data)) {
    stop("Cannot estimate treatment effects with a categorical outcome.")
  }
  outcomes <- FLAME_out$data$outcome[MG]
  treated <- FLAME_out$data$treated[MG] == 1
  CATE <- mean(outcomes[treated]) - mean(outcomes[!treated])
  return(CATE)
}

MG_internal <- function(FLAME_out, MG, which_MG) {
  n_cols <- ncol(FLAME_out$data)
  col_names <- colnames(FLAME_out$data)

  cov_names <-
    col_names[which(!(col_names %in%
                        c('treated', 'outcome', 'weight', 'matched', 'original_ind')))]

  keep_inds <-
    which(!(colnames(FLAME_out$data[MG, ]) %in% c('matched', 'weight', 'original_ind')))
  tmp <- FLAME_out$data[MG, keep_inds]

  # Keep outcome and treatment
  keep <- c(match(names(FLAME_out$matched_on[[which_MG]]), cov_names),
            n_cols - 3, n_cols - 2)

  # Check this
  if (length(keep) == ncol(tmp)) {
    return(tmp)
  }

  # browser()
  # Necessary only in case when replace = TRUE
  # Don't want to overwrite the * (m) with * and also want 5 (m) to go to * (m), not *
  was_na <- apply(tmp[, -keep, drop = FALSE], 2, was_missing)
  any_was_na <- any(was_na)
  was_not_na <- !was_na
  any_was_not_na <- any(was_not_na)

  was_na <- which(was_na, arr.ind = TRUE)
  was_not_na <- which(was_not_na, arr.ind = TRUE)

  if (any_was_na) {
    tmp[, -keep][was_na] <- '* (m)' # Check that this still works even with one not_star_missing
  }
  if (any_was_not_na) {
    tmp[, -keep][was_not_na] <- '*' # Check that this still works even with one not_star_missing
  }
  return(tmp)
}

was_missing <- function(chars) {
  tmp <- strsplit(chars, ' ')
  out <- sapply(tmp, function(char) {
    if (char[length(char)] == '(m)') {
      TRUE
    }
    else {
      FALSE
    }
  })
  return(out)
}

#' Matched Groups
#'
#' \code{MG} returns the matched groups of the supplied units.
#'
#' By default, \code{MG} returns the covariate, treatment, and outcome
#' information for all the units in the relevant matched groups. If only the
#' indices of units in the matched groups are desired, \code{index_only} can be
#' set to \code{TRUE}.
#'
#' Setting \code{multiple = TRUE} will request that all matched groups be
#' returned for each unit -- if \code{\link{FLAME}} was run with \code{replace =
#' TRUE} to generate \code{FLAME_out} in the first place. Otherwise, if
#' \code{\link{FLAME}} was run with \code{replace = TRUE}, but \code{multiple =
#' FALSE}, only main matched groups will be returned. The main matched group of
#' a unit contains the first units it matches with (and therefore those with
#' which it matches on the largest number of covariates). If \code{\link{FLAME}}
#' was run with \code{replace = FALSE}, then the user should only supply
#' \code{multiple = FALSE}.
#'
#' Additionally, if \code{\link{FLAME}} was run with \code{missing_data = 2} to
#' generate \code{FLAME_out}, then \code{MG} will return matched group
#' information for all \code{missing_data_imputations} imputations.
#'
#' @seealso \code{\link{FLAME}}
#'
#' @param units A vector of indices for the units whose matched groups
#' are desired.
#' @param FLAME_out The output of a call to \code{\link{FLAME}}.
#' @param multiple A logical scalar. If \code{FALSE} (default), then \code{MG}
#'   will only return a main matched group for each unit (the first matched
#'   group that unit was a part of). See below for details.
#' @param index_only A logical scalar. If \code{TRUE} then only the indices of
#' units in each matched group are returned.
#'
#' @return \strong{If passing a single set of matched data}
#'
#'   A list of length \code{length(units)}. Each entry is a data frame (if
#'   \code{multiple = FALSE}) or a list of data frames (if \code{multiple =
#'   TRUE}). For a given entry, these data frames are subsets of \code{data}
#'   passed to \code{\link{FLAME}} to generate \code{FLAME_out}, whose rows
#'   correspond to the units in the matched group(s) of that entry. If a unit
#'   is not matched, the corresponding CATE will be \code{NULL}.
#'
#'   The starred entries (*) in the returned data frames have the same meaning
#'   as in \code{FLAME_out$data}, except for if both \code{multiple = TRUE} and
#'   \code{replace = TRUE}. In this case, if \emph{all} units do not match on a
#'   given covariate, all entries of that covariate will be starred, even though
#'   a subset of the units may have matched on them. This is done so that it is
#'   clear on which covariates these units match.
#'
#'   Note that this is the return format also if passing a single set of
#'   imputed data.
#'
#'   \strong{If passing multiple sets of matched, imputed data}
#'
#'   A list of length \code{length(FLAME_out)}, where each entry has the
#'   structure described above, corresponding to that imputed data set.
#'
#' @export
MG <- function(units, FLAME_out, multiple = FALSE, index_only = FALSE) {
  MGs <- iterate(MG_internal, units, FLAME_out, multiple)

  # If multiple imputations,  multiple dfs
  # If multiple units multiple dfs therein
  # If replace = T, multiple dfs therein

  if (index_only) {
    return(MG_rownames(MGs))
  }
  return(MGs)
}

MG_rownames <- function(MG_out) {
  # Recurses over the (possibly complex) structure of MGs generated by
  # iterate(MG_internal, ...) and isolates the rownames
  if (is.data.frame(MG_out)) {
    return(as.numeric(rownames(MG_out)))
  }
  if (length(MG_out) == 0) {
    return(numeric(0))
  }
  out <- list()
  for (i in 1:length(MG_out)) {
    out <- c(out, list(MG_rownames(MG_out[[i]])))
  }
  return(out)
}

#' Conditional Average Treatment Effects
#'
#' \code{CATE} returns the conditional
#' average treatment effects (CATEs) of \code{units}.
#'
#' The CATE of a matched group is defined to be the difference between average
#' treated and control outcomes within that matched group. When we refer to
#' the CATE(s) of a unit, we mean the CATE(s) of its matched group(s).
#'
#' Setting \code{multiple = TRUE} will request that CATEs corresponding to all
#' matched groups be returned for each unit -- if \code{\link{FLAME}} was run
#' with \code{replace = TRUE} to generate \code{FLAME_out} in the first place.
#' Otherwise, if \code{\link{FLAME}} was run with \code{replace = TRUE}, but
#' \code{multiple = FALSE}, only the CATE of the main matched groups will be
#' returned. The main matched group of a unit contains the first units it
#' matches with (and therefore those with which it matches on the largest number
#' of covariates). If \code{\link{FLAME}} was run with \code{replace = FALSE},
#' then the user should only supply \code{multiple = FALSE}.
#'
#' Additionally, if \code{\link{FLAME}} was run with \code{missing_data = 2} to
#' generate \code{FLAME_out}, then \code{CATE} will return CATE
#' information for all \code{missing_data_imputations} imputations.
#'
#' @seealso \code{\link{FLAME}}
#'
#' @param units A vector of indices for the units whose CATEs
#'   are desired.
#' @param FLAME_out The output of a call to \code{\link{FLAME}}.
#' @param multiple A logical scalar. If \code{FALSE} (default), then \code{CATE}
#'   will return CATEs of main matched groups (those with matches on the
#'   greatest number of covariates). See below for details.
#'
#' @return \strong{If passing a single set of matched data}
#'
#'   A list of length \code{length(units)}. Each entry is a CATE (a numeric
#'   scalar) (if \code{multiple = FALSE}) or a list of CATEs (if \code{multiple
#'   = TRUE}). If a unit is not matched, the corresponding CATE will be
#'   \code{NULL}.
#'
#'   Note that this is the return format also if passing a single set of
#'   imputed data.
#'
#'   \strong{If passing multiple sets of matched, imputed data}
#'
#'   A list of length \code{length(FLAME_out)}, where each entry has the
#'   structure described above, corresponding to that imputed data set.
#'
#' @export
CATE <- function(units, FLAME_out, multiple = FALSE) {
  return(iterate(CATE_internal, units, FLAME_out, multiple))
}

#' ATE of a matched dataset
#'
#' \code{ATE} computes the average treatment effect (ATE) of a matched dataset.
#'
#' The ATE is computed as the difference between the weighted treated and the
#' weighted control outcomes in the dataset. A unit's weight is the number of
#' times it was matched.
#'
#' @param FLAME_out An object returned by running \code{\link{FLAME}}
#' @export
ATE <- function(FLAME_out) {

  if (is.null(names(FLAME_out))) { # Is a list of data frames
    n_df <- length(FLAME_out)
  }
  else { # Is a single data frame
    FLAME_out <- list(FLAME_out)
    n_df <- 1
  }

  if (!('outcome' %in% colnames(FLAME_out[[1]]$data))) {
    stop('Outcome was not supplied with `data`; cannot estimate ATE.')
  }

  if (is.factor(FLAME_out[[1]]$data)) {
    stop("Cannot estimate treatment effects with a categorical outcome.")
  }

  out <- vector('numeric', length = n_df)

  for (i in 1:n_df) {
    weight <- FLAME_out[[i]]$data$weight
    CATEs <- FLAME_out[[i]]$CATE
    MGs <- FLAME_out[[i]]$MGs

    weight_sum <- 0
    weighted_CATE_sum <- 0

    for (j in 1:length(MGs)) {
      MG_weight <- sum(weight[MGs[[j]]])
      weight_sum <- weight_sum + MG_weight
      weighted_CATE_sum <- weighted_CATE_sum + MG_weight * CATEs[[j]]
    }

    ATE <- weighted_CATE_sum / weight_sum
    out[i] <- ATE
  }
  return(out)
}

#' ATT of a matched dataset
#'
#' \code{ATT} computes the average treatment effect on the treated (ATT) of a
#' matched dataset.
#'
#' The counterfactual outcome of each treated unit is estimated via the mean
#' outcome of control units in its matched group. The difference between
#' observed treated outcome and estimated counterfactual outcome is then
#' averaged -- with weights given by the number of control units in that treated
#' unit's matched group -- across all treated units to compute the ATT.
#' @param FLAME_out An object returned by running \code{\link{FLAME}}
#' @export
ATT <- function(FLAME_out) {
  if (is.null(names(FLAME_out))) { # Is a list of data frames
    n_df <- length(FLAME_out)
  }
  else { # Is a single data frame
    FLAME_out <- list(FLAME_out)
    n_df <- 1
  }

  if (!('outcome' %in% colnames(FLAME_out[[1]]$data))) {
    stop('Outcome was not supplied with `data`; cannot estimate ATT.')
  }

  if (is.factor(FLAME_out[[1]]$data)) {
    stop("Cannot estimate treatment effects with a categorical outcome.")
  }

  out <- vector('numeric', length = n_df)

  controls <- which(FLAME_out[[1]]$data$treated == 0)
  treated <- which(FLAME_out[[1]]$data$treated == 1)

  outcomes <- FLAME_out[[1]]$data$outcome

  for (i in 1:n_df) {
    weight <- FLAME_out[[i]]$data$weight

    MGs <- FLAME_out[[i]]$MGs

    weight_sum <- 0
    weighted_TT_sum <- 0

    for (j in 1:length(MGs)) {

      MG_controls <- MGs[[j]][MGs[[j]] %in% controls]
      MG_treated <- MGs[[j]][MGs[[j]] %in% treated]

      MG_weight <- sum(weight[MG_controls])
      weight_sum <- weight_sum + MG_weight * length(MG_treated)

      mean_control_outcome <- mean(outcomes[MG_controls])

      for (k in seq_along(MG_treated)) {
        weighted_TT_sum <-
          weighted_TT_sum +
          MG_weight * (outcomes[MG_treated[k]] - mean_control_outcome)
      }
    }

    ATT <- weighted_TT_sum / weight_sum
    out[i] <- ATT
  }
  return(out)
}
