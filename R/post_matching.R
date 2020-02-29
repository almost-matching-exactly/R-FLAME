#' Returns the matched group of a given unit.
#'
#' \code{mmg_of_unit} Returns the matched group of unit \code{unit}, as well as
#' the covariates of the units in the matched group.
#'
#' @param flame_obj An object returned by running \code{FLAME_bit}.
#' @param unit An integer between 1 and \code{nrow(flame_obj)} for whom the
#'   treatment effect is desired.
#' @param output_style If 1, only returns the values of covariates these units
#'   were matched on. If 0, returns all covariate values.
#' @export
mmg_of_unit <- function(unit, flame_obj, output_style = 1) {
  MGs <- flame_obj$MGs
  for (i in 1:length(MGs)) {
    if (unit %in% MGs[[i]]) {
      if (output_style == 1) {
        return(flame_obj$data[MGs[[i]], ])
      }
      return(flame_obj$original_data[MGs[[i]], ])
    }
  }
  stop('Unit not matched')
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
  weights <- flame_obj$data$weights
  outcomes <- flame_obj$data$outcome

  control <- flame_obj$data$treated == 0
  treated <- flame_obj$data$treated == 1

  ATE <-
  sum(outcomes[treated] * weights[treated]) / sum(weights[treated]) -
    sum(outcomes[control] * weights[control]) / sum(weights[control])
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
