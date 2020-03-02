#' Matched groups of units.
#'
#' Returns the matched groups of units \code{units}, as well as the covariates,
#' treatment, and outcome of the units in each matched group.
#'
#' \code{MG} returns a list, each entry of which contains information about the
#' matched group(s) of the corresponding unit. If \code{\link{FLAME}} was run with
#' argument \code{repeats = FALSE} to generate \code{flame_obj}, \code{MG}
#' should be run with the default \code{multiple = FALSE} and will return
#' information about the matched group of each unit. If \code{link{FLAME}} was run
#' with argument \code{repeats = TRUE}, then \code{MG} may be run with either
#' \code{multiple = FALSE} or \code{multiple = TRUE}. In the first case, only
#' one matched group, that with matches on the greatest number of covariates,
#' will be returned. In the second case, all matched groups involving the
#' relevant unit will be returned.
#'
#' @seealso \code{\link{FLAME}}
#'
#' @param units A vector of integers between 1 and \code{nrow(flame_obj$data)}
#'   whose matched groups are desired.
#' @param flame_obj An object returned by running \code{\link{FLAME}}.
#' @param multiple A logical scalar. If \code{FALSE} (default), then will only
#'   return each unit's first matched group; else, will return all of them. See
#'   below for details.
#'
#' @return A list of length \code{length(units)}. Each entry is a data frame (if
#'   \code{multiple = FALSE}) or a list of data frames (if \code{multiple =
#'   TRUE}). For a given entry, these data frames are subsets of \code{data}
#'   passed to \code{\link{FLAME}} to generate \code{flame_obj}, whose rows
#'   correspond to the units in the matched group(s) of that entry.
#'
#'   Each entry thus always corresponds to the matched group(s) of the
#'   corresponding unit in \code{units}, but the structure of the entries
#'   depends on the value of \code{multiple} and on whether or not the call to
#'   \code{\link{FLAME}} that generated \code{flame_obj} specified \code{repeats =
#'   TRUE} or not. If it did, a unit may have several matched
#'   groups, corresponding to several sets of covariates it was matched on. In
#'   this case, \code{multiple = FALSE} requests that only the matched group
#'   with matches on the greatest number of covariates be returned and
#'   \code{multiple = TRUE} requests that all matched groups involving that unit
#'   be returned. If it did not, a unit only has one matched group
#'   that will be returned.
#'
#'   The starred entries (*) in the returned data frames have the same meaning
#'   as in \code{flame_obj$data}, except for if both \code{multiple = TRUE} and
#'   \code{repeats = TRUE}. In this case, if \emph{all} units do not match on a given
#'   covariate, all entries of that covariate will be starred, even though a
#'   subset of the units may have matched on them.
#'
#'
#' @export
#'
MG <- function(units, flame_obj, multiple = FALSE) {
  n_cols <- length(colnames(flame_obj$data))
  cov_names <-
    colnames(flame_obj$data)[1:(n_cols - 4)]
  MGs <- flame_obj$MGs # MGs made by FLAME

  out <- vector('list', length = length(units))

  for (i in seq_along(units)) {
    unit <- units[i]

    # Number of MGs to return for unit
    # Only > 1 if multiple = TRUE and matched multiple times
    n_unit_MGs <- ifelse(multiple, flame_obj$data$weight[unit], 1)

    # To store MGs of unit
    out[[i]] <- vector('list', length = n_unit_MGs)
    counter <- 1
    for (j in 1:length(MGs)) {
      if (unit %in% MGs[[j]]) {
        tmp <-
          flame_obj$data[MGs[[j]], ] %>%
          dplyr::select(-c(matched, weight))
        keep <-
          match(names(flame_obj$matched_on[[j]]), cov_names) %>%
          c(n_cols - 3, n_cols - 2) # Keep outcome and treatment
        tmp[, -keep] <- '*'
        out[[i]][[counter]] <- tmp

        counter <- counter + 1
        if (counter == n_unit_MGs + 1) {
          break
        }
      }
    }
  }
  if (!multiple) {
    out %<>% lapply(`[[`, 1)
  }
  return(out)
}

#' Compute the treatment effect of a given unit.
#'
#' \code{te_of_unit} Computes the treatment effect of unit \code{unit} as equal
#' to the CATE of its matched group.
#'
#' @param flame_obj An object returned by running \code{FLAME_bit}.
#' @param unit An integer between 1 and \code{nrow(flame_obj)} for whom the
#'   treatment effect is desired.
#' @export
te_of_unit <- function(unit, flame_obj) {
  MGs <- flame_obj$MGs
  for (i in 1:length(MGs)) {
    if (unit %in% MGs[[i]]) {
      return(flame_obj$CATE[i])
    }
  }
  stop('Unit not matched')
}

#' @export
CATE <- function(units, flame_obj) {
  if (any(!(units %in% 1:nrow(flame_obj$data)))) {
    stop('Unit not in dataset')
  }
  matched_units <- units[units %in% which(flame_obj$data$matched)]
  if (length(matched_units) == 0) {
    stop('None of the supplied units were matched by FLAME')
  }

  MGs <- flame_obj$MGs
  sizes <- NULL
  weighted_CATE_sum <- 0
  for (unit in matched_units) {
    for (i in 1:length(MGs)) {
      if (unit %in% MGs[[i]]) {
        MG_size <- length(MGs[[i]])
        sizes %<>% c(MG_size)
        weighted_CATE_sum %<>% add(flame_obj$CATE[i] * MG_size)
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
#' @param flame_obj An object returned by running \code{FLAME_bit}.
#' @export
ATE <- function(flame_obj) {
  weight <- flame_obj$data$weight
  outcomes <- flame_obj$data$outcome

  control <- flame_obj$data$treated == 0
  treated <- flame_obj$data$treated == 1

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
#' @param flame_obj An object returned by running \code{FLAME_bit}.
#' @export
ATT <- function(flame_obj) {
  matched_treated <- with(flame_obj$data, which(matched & treated == 1))
  controls <- which(flame_obj$data$treated == 0)
  outcomes <- flame_obj$data$outcome
  MGs <- flame_obj$MGs
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
