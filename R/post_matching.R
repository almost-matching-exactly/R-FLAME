#' Matched Groups
#'
#' \code{MG} returns the matched groups of the supplied units.
#'
#' The \code{units} argument refers to units with respect to
#' \code{rownames(ame_out$data)}. Typically, this will also correspond to the
#' indexing of the data (i.e. passing \code{units = 3} will return the matched
#' group of the 3rd unit in the matching data). However, if a separate holdout
#' set was not passed to the matching algorithm or if the original matching data
#' had rownames other than \code{1:nrow(data)}, then this is not the case.
#'
#' The \code{multiple} argument toggles whether only a unit's main matched group
#' (MMG) or all matched groups a unit is part of should be returned. A unit's
#' MMG contains its highest quality matches (that is, the units with which it
#' first matched in the sequence of considered covariate sets). If the original
#' call that generated \code{ame_out} specified \code{replace = FALSE} then
#' units only are part of one matched group (which is also their MMG) and
#' \code{multiple} must be set to \code{FALSE}.
#'
#' @param units A vector of units whose matched groups are desired.
#' @param ame_out An object of class \code{ame}.
#' @param multiple A logical scalar. If \code{FALSE} (default), then \code{MG}
#'   will only return the main matched group for each unit. See below for
#'   details. Cannot be set to \code{TRUE} if \code{ame_out} was generated
#'   without replacement.
#' @param id_only A logical scalar. If \code{TRUE}, then only the IDs of the
#'   units in each matched group are returned, and not their treatment, outcome,
#'   or covariate information.
#' @param index_only Defunct. Use `id_only` instead.
#'
#' @return
#'
#'   A list of length \code{length(units)}, each entry of which corresponds to a
#'   different unit in \code{units}. For matched units, if \code{multiple =
#'   FALSE}, each entry is 1. a data frame containing the treatment and outcome
#'   information of members of the matched group, along with covariates they
#'   were matched on if \code{id_only = FALSE} or 2. a vector of the IDs of
#'   matched units if \code{id_only = TRUE} . If \code{multiple = TRUE}, each
#'   entry of the returned list is a list containing the previously described
#'   information, but with each entry corresponding to a different matched
#'   group. In either case, entries corresponding to unmatched units are
#'   \code{NULL}.
#' @examples
#' \dontrun{
#' data <- gen_data()
#' holdout <- gen_data()
#' FLAME_out <- FLAME(data = data, holdout = holdout, replace = TRUE)
#'
#' # Only the main matched group of unit 1
#' MG(1, FLAME_out, multiple = F)
#'
#' # All matched groups of unit 1
#' MG(1, FLAME_out, multiple = T)
#' }
#' @export
MG <- function(units, ame_out, multiple = FALSE,
               id_only = FALSE, index_only) {

  if (!missing(index_only)) {
    stop('Argument `index_only` is defunct and will be removed in a later',
         'release; please use `id_only` instead.')
  }

  if (multiple && !ame_out$info$replacement) {
    stop(paste('Multiple matched groups cannot be queried if',
               ame_out$info$algo, 'was run without replacement.'))
  }

  if (any(!(units %in% as.numeric(rownames(ame_out$data))))) {
    stop('Supplied a unit not in the matched data.')
  }

  # In the case that the rownames of the original data frame were not 1:n
  # or that a holdout set was not explicitly passed to the algo
  units <- match(units, rownames(ame_out$data))

  if (id_only & !multiple) {
    return(lapply(ame_out$MGs[units], function(z) {
      if (is.null(z)) {
        NULL
      }
      else {
        rownames(ame_out$data)[z]
      }
    }))
  }

  outcome_name <- ame_out$info$outcome
  treated_name <- ame_out$info$treatment

  col_names <- colnames(ame_out$data)

  cov_inds <-
    which(!(col_names %in%
              c(treated_name, outcome_name, 'weight', 'matched', 'MG', 'CATE')))

  all_MGs <- vector('list', length(units))
  all_MGs <- lapply(all_MGs, function(x) x <- vector('list', 1))

  for (i in seq_along(units)) {
    if (is.null(ame_out$MGs[[units[i]]])) {
      next
    }

    if (!multiple) {
      MGs_to_return <- units[i]
    }
    else {
      in_MG <-
        which(vapply(ame_out$MGs, function(x) units[i] %in% x, logical(1)))
      MGs_to_return <- in_MG[!duplicated(ame_out$MGs[in_MG])]
    }

    for (j in seq_along(MGs_to_return)) {

      if (id_only) {
        all_MGs[[i]][[j]] <-
          rownames(ame_out$data)[ame_out$MGs[[MGs_to_return[j]]]]
        next
      }

      MG_all_cols <-
        ame_out$data[ame_out$MGs[[MGs_to_return[j]]], , drop = FALSE]

      # Can't just iterate one covariate at a time because for FLAME the
      # following might happen: 1. two units don't match exactly, 2. a covariate
      # they have identical values of is dropped, 3. the covariate they differ
      # on is now also dropped. In this case, the units technically were not
      # matched on that first dropped covariate. You could argue it should be
      # shown anyway, but at least for now, for consistency with the python
      # code, we'll ignore it.
      for (cov_set in ame_out$cov_sets) {
        matched_on <- setdiff(col_names[cov_inds], cov_set)
        data <- data.matrix(MG_all_cols[, matched_on, drop = FALSE])
        tmp <- colSums(sweep(data, 2, data[1, ]) ^ 2)
        if (any(is.na(tmp))) {
          next
        }
        if (all(tmp == 0)) {
          all_MGs[[i]][[j]] <-
            MG_all_cols[, c(matched_on, treated_name, outcome_name)]
          break
        }
      }
    }
  }

  if (!multiple) {
    all_MGs <- lapply(all_MGs, `[[`, 1)
  }

  return(all_MGs)
}

#' Conditional Average Treatment Effects
#'
#' \code{CATE} returns an estimate of the conditional average treatment effect
#' for the subgroup defined by \code{units}.
#'
#'This function returns CATE estimates and the estimated variances of such
#'estimates for \code{units} of interest. The CATE of a given unit is estimated
#'by the difference between the average treated and the average control outcome
#'in that unit's main matched group. As shown by Morucci 2021, under standard
#'regularity conditions, such an estimator is asymptotically normal and unbiased
#'for the true CATE, with a variance that can be estimated by the sum of the
#'variance of treated and control outcomes in the matched group, each normalized
#'by the number of treated and control units in the matched group, respectively.
#'Note that CATEs cannot be estimated for unmatched units and estimator
#'variances cannot be computed for units whose main matched group only contains
#'a single treated or control unit. Note also that these CATE estimates differ
#'from those that are used to compute average treatment effects in
#'\code{print.ame} and \code{summary.ame} and from those that will be returned
#'in \code{ame_out$data$CATE} if \code{estimate_CATEs = TRUE}. For a treated
#'(control) unit \eqn{i}, the latter estimate the treated (control) outcome
#'conditioned on \eqn{X = x_i} simply as \eqn{Y_i}, and do not average across
#'other treated (control) units in the matched group as is done here. This
#'averaging is necessary in order to compute variance estimates. The different
#'estimates can always be manually compared, though they are the same in
#'expectation (assuming mean 0 noise) and we expect them to be similar in
#'practice, in the absence of large noise.
#'
#' Lastly, note that the \code{units} argument refers to units with respect to
#' \code{rownames(ame_out$data)}. Typically, this will also correspond to the
#' indexing of the data (i.e. passing \code{units = 3} will return the matched
#' group of the 3rd unit in the matching data). However, if a separate holdout
#' set was not passed to the matching algorithm or if the original matching data
#' had rownames other than \code{1:nrow(data)}, then this is not the case.
#'
#' @seealso \code{\link{FLAME}}, \code{\link{DAME}}
#'
#' @param units A vector of units whose CATE estimates are desired.
#' @param ame_out An object of class \code{ame}.
#'
#' @return A matrix whose columns correspond to CATE estimates and their
#'   variances and whose rows correspond to queried units. \code{NA}'s therein
#'   correspond to inestimable quantities.
#'
#' @examples
#' \dontrun{
#' data <- gen_data()
#' holdout <- gen_data()
#' FLAME_out <- FLAME(data = data, holdout = holdout)
#' CATE(1:5, FLAME_out)
#' }
#' @export
CATE <- function(units, ame_out) {


  if (any(!(units %in% as.numeric(rownames(ame_out$data))))) {
    stop('Supplied a unit not in the matched data.')
  }

  # In the case that the rownames of the original data frame were not 1:n
  # or that a holdout set was not explicitly passed to the algo
  units <- match(units, rownames(ame_out$data))

  if (ame_out$info$outcome_type != 'continuous') {
    stop('CATEs and their variance estimates can only be computed for ',
         'continuous outcomes.')
  }

  CATE_mat <- matrix(nrow = length(units), ncol = 2,
                     dimnames = list(units, c('Mean', 'Variance')))

  Tr <- ame_out$data[[ame_out$info$treatment]]
  Y <- ame_out$data[[ame_out$info$outcome]]
  n <- length(Y)

  for (unit in units) {
    mg <- ame_out$MGs[[unit]]
    if (is.null(mg)) {
      next
    }

    control <- mg[Tr[mg] == 0]
    treated <- mg[Tr[mg] == 1]

    CATE_mat[as.character(unit), 'Mean'] <- mean(Y[treated]) - mean(Y[control])
    CATE_mat[as.character(unit), 'Variance'] <-
      var(Y[treated]) / length(treated) +
      var(Y[control]) / length(control)
  }

  return(CATE_mat)
}

#' Average Treatment Effect estimates
#'
#' These functions are deprecated and will be made defunct at a later release.
#' See \code{summary.ame} for average treatment effects estimates and their
#' variance.
#'
#' \code{ATE}, \code{ATT}, and \code{ATC} estimate the average treatment effect
#' (ATE), average treatment effect on the treated (ATT), and average treatment
#' effect on the controls (ATC), respectively, of a matched dataset.
#'
#' The ATE is estimated as the average CATE estimate across all matched units in
#' the data, while the ATT and ATC average only across matched treated or
#' matched control units, respectively.
#' @seealso \code{\link{CATE}}
#' @param ame_out An object of class \code{ame}.
#' @export
ATE <- function(ame_out) {
  .Deprecated(msg = paste('`ATE` is now deprecated, as ATE, ATT, and ATC mean',
                          'and variance estimates are now computed in the',
                          '`summary.ame` method'))
  if (ame_out$info$estimate_CATEs &&
      ame_out$info$outcome_type == 'continuous') {
    return(mean(ame_out$data$CATE, na.rm = TRUE))
  }
  return(get_average_effects(ame_out)['All', 'Mean'])
}

#' @rdname ATE
#' @export
ATT <- function(ame_out) {
  .Deprecated(msg = paste('`ATT` is now deprecated, as ATE, ATT, and ATC mean',
                          'and variance estimates are now computed in the',
                          '`summary.ame` method'))
  if (ame_out$info$estimate_CATEs &&
      ame_out$info$outcome_type == 'continuous') {
    return(mean(ame_out$data$CATE[ame_out$data[[ame_out$info$treatment]] == 1],
                na.rm = TRUE))
  }
  return(get_average_effects(ame_out)['Treated', 'Mean'])
}


#' @rdname ATE
#' @export
ATC <- function(ame_out) {
  .Deprecated(msg = paste('`ATC` is now deprecated, as ATE, ATT, and ATC mean',
                          'and variance estimates are now computed in the',
                          '`summary.ame` method'))

  if (ame_out$info$estimate_CATEs &&
      ame_out$info$outcome_type == 'continuous') {
   return(mean(ame_out$data$CATE[ame_out$data[[ame_out$info$treatment]] == 0],
               na.rm = TRUE))
  }
  return(get_average_effects(ame_out)['Control', 'Mean'])
}
