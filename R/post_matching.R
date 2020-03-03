#' Matched groups of units.
#'
#' Returns the matched groups of units, as well as the covariates,
#' treatment, and outcome of the units in each matched group.
#'
#' \code{MG} returns the covariates, treatments, and outcomes of the units in
#' the matched groups of \code{units}. Multiple matched groups per unit may be
#' returned if their creation was allowed when originally running
#' \code{\link{FLAME}}. Additionally, if \code{\link{FLAME}} was run with
#' \code{missing_data = 2} to generate \code{FLAME_out}, then \code{MG} will
#' return matched group information for all \code{missing_data_imputations}
#' imputations.
#'
#' @seealso \code{\link{FLAME}}
#'
#' @param units A vector of integers between 1 and \code{nrow(FLAME_out$data)}
#'   whose matched groups are desired.
#' @param FLAME_out The output of a call to \code{\link{FLAME}}.
#' @param multiple A logical scalar. If \code{FALSE} (default), then will only
#'   return the matched group of each unit with matches on the greatest number
#'   of covariates; otherwise, will return all of them. See below for details.
#'
#' @return \strong{If passing a single set of matched data}
#'
#'   A list of length \code{length(units)}. Each entry is a data frame (if
#'   \code{multiple = FALSE}) or a list of data frames (if \code{multiple =
#'   TRUE}). For a given entry, these data frames are subsets of \code{data}
#'   passed to \code{\link{FLAME}} to generate \code{FLAME_out}, whose rows
#'   correspond to the units in the matched group(s) of that entry.
#'
#'   Each entry thus always corresponds to the matched group(s) of the
#'   corresponding unit in \code{units}, but the structure of the entries
#'   depends on the value of \code{multiple} and on whether or not the call to
#'   \code{\link{FLAME}} that generated \code{FLAME_out} specified \code{repeats
#'   = TRUE} or not. If it did, a unit may have several matched groups,
#'   corresponding to several sets of covariates it was matched on. In this
#'   case, \code{multiple = FALSE} requests that only the matched group with
#'   matches on the greatest number of covariates be returned and \code{multiple
#'   = TRUE} requests that all matched groups involving that unit be returned.
#'   If the call to \code{\link{FLAME}} specified \code{repeats = FALSE}, a unit
#'   only has one matched group that will be returned. Regardless of the value
#'   of \code{repeats}, if \code{multiple = TRUE}, each entry of the returned
#'   list will be a list of one or more data frames (matched groups). And if
#'   \code{multiple = FALSE}, each entry will be a single data frame (matched
#'   group).
#'
#'   The starred entries (*) in the returned data frames have the same meaning
#'   as in \code{FLAME_out$data}, except for if both \code{multiple = TRUE} and
#'   \code{repeats = TRUE}. In this case, if \emph{all} units do not match on a
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
#'
MG <- function(units, FLAME_out, multiple = FALSE) {
  if (is.null(names(FLAME_out))) { # Is a list of data frames
    n_df <- length(FLAME_out)
  }
  else { # Is a single data frame
    FLAME_out <- list(FLAME_out)
    n_df <- 1
  }

  n_cols <- length(colnames(FLAME_out[[1]]$data))
  cov_names <-
    colnames(FLAME_out[[1]]$data)[1:(n_cols - 4)]

  out <- vector('list', length = n_df)

  for (k in 1:n_df) {
    MGs <- FLAME_out[[k]]$MGs # MGs made by FLAME

    out[[k]] <- vector('list', length = length(units))

    for (i in seq_along(units)) {
      unit <- units[i]

      # Number of MGs to return for unit
      # Only > 1 if multiple = TRUE and matched multiple times
      n_unit_MGs <- ifelse(multiple, FLAME_out[[k]]$data$weight[unit], 1)

      # To store MGs of unit
      out[[k]][[i]] <- vector('list', length = n_unit_MGs)
      counter <- 1
      for (j in 1:length(MGs)) {
        if (unit %in% MGs[[j]]) {
          tmp <-
            FLAME_out[[k]]$data[MGs[[j]], ] %>%
            dplyr::select(-c(matched, weight))
          keep <-
            match(names(FLAME_out[[k]]$matched_on[[j]]), cov_names) %>%
            c(n_cols - 3, n_cols - 2) # Keep outcome and treatment
          tmp[, -keep] <- '*'
          out[[k]][[i]][[counter]] <- tmp

          counter <- counter + 1
          if (counter == n_unit_MGs + 1) {
            break
          }
        }
      }
    }
    if (!multiple) {
      out[[k]] %<>% lapply(`[[`, 1)
    }
  }
  if (n_df == 1) {
    out <- out[[1]]
  }
  return(out)
}

#' Compute the treatment effect of a given unit.
#'
#' \code{te_of_unit} Computes the treatment effect of unit \code{unit} as equal
#' to the CATE of its matched group.
#'
#' @param FLAME_out An object returned by running \code{FLAME_bit}.
#' @param unit An integer between 1 and \code{nrow(FLAME_out)} for whom the
#'   treatment effect is desired.
#' @export
te_of_unit <- function(unit, FLAME_out) {
  MGs <- FLAME_out$MGs
  for (i in 1:length(MGs)) {
    if (unit %in% MGs[[i]]) {
      return(FLAME_out$CATE[i])
    }
  }
  stop('Unit not matched')
}

#' @export
CATE <- function(units, FLAME_out) {
  if (any(!(units %in% 1:nrow(FLAME_out$data)))) {
    stop('Unit not in dataset')
  }
  matched_units <- units[units %in% which(FLAME_out$data$matched)]
  if (length(matched_units) == 0) {
    stop('None of the supplied units were matched by FLAME')
  }

  MGs <- FLAME_out$MGs
  sizes <- NULL
  weighted_CATE_sum <- 0
  for (unit in matched_units) {
    for (i in 1:length(MGs)) {
      if (unit %in% MGs[[i]]) {
        MG_size <- length(MGs[[i]])
        sizes %<>% c(MG_size)
        weighted_CATE_sum %<>% add(FLAME_out$CATE[i] * MG_size)
        break
      }
    }
  }
  return(weighted_CATE_sum / sum(sizes))
}

#' Compute the ATE of a matched dataset.
#'
#' \code{ATE} Computes the average treatment effect (ATE) of a matched dataset
#' via the difference of the weighted treated and the weighted control outcomes
#' A unit's weight is the number of times it was matched.
#'
#' @param FLAME_out An object returned by running \code{FLAME_bit}.
#' @export
ATE <- function(FLAME_out) {
  weight <- FLAME_out$data$weight
  outcomes <- FLAME_out$data$outcome

  control <- FLAME_out$data$treated == 0
  treated <- FLAME_out$data$treated == 1

  ATE <-
  sum(outcomes[treated] * weight[treated]) / sum(weight[treated]) -
    sum(outcomes[control] * weight[control]) / sum(weight[control])
  return(ATE)
}

#' Compute the ATT of a matched dataset.
#'
#' \code{ATT} Computes the average treatment effect on the treated (ATT) of a matched dataset.
#'
#' For each treated unit, its counterfactual outcome is estimated via the mean
#' outcome of control units in its matched group. This value is then averaged
#' across all treated units.
#' @param FLAME_out An object returned by running \code{FLAME_bit}.
#' @export
ATT <- function(FLAME_out) {
  matched_treated <- with(FLAME_out$data, which(matched & treated == 1))
  controls <- which(FLAME_out$data$treated == 0)
  outcomes <- FLAME_out$data$outcome
  MGs <- FLAME_out$MGs
  sum_TT <- 0
  for (j in matched_treated) {
    for (i in 1:length(MGs)) {
      MG <- MGs[[i]]
      if (j %in% MG) {
        controls <- MG[MG %in% controls]
        sum_TT %<>% add(outcomes[j] - mean(outcomes[controls]))
        break
      }
    }
  }
  return(sum_TT / length(matched_treated))
}
