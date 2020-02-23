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
te_of_unit <- function(unit, flame_obj) {
  MGs <- flame_obj$MGs
  for (i in 1:length(MGs)) {
    if (unit %in% MGs[[i]]) {
      return(flame_obj$CATE[i])
    }
  }
  stop('Unit not matched')
}

#' Compute the ATE of a matched dataset.
#'
#' \code{ATE} Computes the ATE of a matched dataset via a weighted (by size)
#' average of the CATEs of each matched group in the matched dataset.
#'
#' @param flame_obj An object returned by running \code{FLAME_bit}.
ATE <- function(flame_obj) {
  CATE <- flame_obj$CATE
  MGs <- flame_obj$MGs
  size <- sapply(MGs, length)
  return(sum(size * CATE) / sum(size))
}
