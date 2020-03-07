#' Matched groups of units.
#'
#' Returns the matched groups of units, as well as the covariates,
#' treatment, and outcome of the units in each matched group.
#'
#' \code{MG} returns the covariates, treatments, and outcomes of the units in
#' the matched groups of \code{units}. Multiple matched groups per unit may be
#' returned if their creation was allowed when originally running
#' \code{\link{FLAME}} with \code{replace = TRUE}. Additionally, if
#' \code{\link{FLAME}} was run with \code{missing_data = 2} to generate
#' \code{FLAME_out}, then \code{MG} will return matched group information for
#' all \code{missing_data_imputations} imputations.
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
#'   \code{\link{FLAME}} that generated \code{FLAME_out} specified \code{replace
#'   = TRUE} or not. If it did, a unit may have several matched groups,
#'   corresponding to several sets of covariates it was matched on. In this
#'   case, \code{multiple = FALSE} requests that only the matched group with
#'   matches on the greatest number of covariates be returned and \code{multiple
#'   = TRUE} requests that all matched groups involving that unit be returned.
#'   If the call to \code{\link{FLAME}} specified \code{replace = FALSE}, a unit
#'   only has one matched group that will be returned. Regardless of the value
#'   of \code{replace}, if \code{multiple = TRUE}, each entry of the returned
#'   list will be a list of one or more data frames (matched groups). And if
#'   \code{multiple = FALSE}, each entry will be a single data frame (matched
#'   group).
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

group_CATE <- function(units, FLAME_out, multiple = FALSE) {
  if (is.null(names(FLAME_out))) { # Is a list of data frames
    n_df <- length(FLAME_out)
  }
  else { # Is a single data frame
    FLAME_out <- list(FLAME_out)
    n_df <- 1
  }

  out <- vector('list', length = n_df)

  for (k in 1:n_df) {
    MGs <- FLAME_out[[k]]$MGs # MGs made by FLAME
    CATEs <- FLAME_out[[k]]$CATE

    out[[k]] <- vector('list', length = length(units))

    out[[k]] <- lapply(units, function(x) {
      CATEs[which(sapply(MGs, function(y) x %in% y))]
    })
    # for (i in seq_along(units)) {
    #
    #   # To store CATEs of unit
    #   in_MG <- which(sapply(MGs, function(x) units[i] %in% x))
    #   out[[k]][[i]] <- CATEs[in_MG]
    # }
    # if (!multiple) { # Check this for replace = FALSE, multiple = FALSE
    #   out[[k]] %<>% lapply(`[[`, 1)
    # }
  }
  if (n_df == 1) {
    out <- out[[1]]
  }
  return(out)
}

#' Compute the treatment effect of a given unit.
#'
#' \code{CATE} Computes the treatment effect of unit \code{unit} as equal
#' to the CATE of its matched group.
#'
#' @param FLAME_out An object returned by running \code{FLAME_bit}.
#' @param units A vector of integers between 1 and \code{nrow(FLAME_out$data)}
#'   whose matched groups are desired.
#' @export

CATE <- function(units, FLAME_out) {

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

  if (is.null(names(FLAME_out))) { # Is a list of data frames
    n_df <- length(FLAME_out)
  }
  else { # Is a single data frame
    FLAME_out <- list(FLAME_out)
    n_df <- 1
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
  if (is.null(names(FLAME_out))) { # Is a list of data frames
    n_df <- length(FLAME_out)
  }
  else { # Is a single data frame
    FLAME_out <- list(FLAME_out)
    n_df <- 1
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
      weight_sum <- weight_sum + MG_weight

      mean_control_outcome <- mean(outcomes[MG_controls])

      for (k in seq_along(MG_treated)) {
        weighted_TT_sum <-
          weighted_TT_sum + MG_weight * (outcomes[MG_treated[k]] - mean_control_outcome)
      }
    }

    ATT <- weighted_TT_sum / weight_sum
    out[i] <- ATT
  }
  return(out)
}
