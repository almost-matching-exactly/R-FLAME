#' Matched Groups
#'
#' \code{MG} returns the matched groups of the supplied units.
#'
#' \code{MG} returns the treatment, outcome, and covariates matched on for all
#' units in the relevant matched groups. If only the indices of units in the
#' matched groups are desired, \code{index_only} can be set to \code{TRUE}.
#'
#' The \code{multiple} argument toggles whether only a unit's main matched group
#' (MMG) or all matched groups a unit is part of should be returned. A unit's
#' MMG contains its highest quality matches. Note that if the original call that
#' generated \code{AME_out} specified \code{replace = FALSE} then units only are
#' part of one matched group (their MMG) and this argument has no effect.
#' Otherwise, \code{multiple = TRUE} requests that \code{MG} returns all matched
#' groups and \code{multiple = FALSE} that only MMGs are returned.
#'
#' @seealso \code{\link{FLAME}}, \code{\link{DAME}}
#'
#' @param units A vector of indices for the units whose matched groups
#' are desired.
#' @param AME_out An object of class \code{AME}.
#' @param multiple A logical scalar. If \code{FALSE} (default), then \code{MG}
#'   will only return the main matched group for each unit. See below for
#'   details. Should not be set to \code{TRUE} if \code{AME_out} was generated
#'   without replacement.
#' @param index_only A logical scalar. If \code{TRUE} then only the indices of
#' units in each matched group are returned.
#'
#' @return
#'
#'   A list of length \code{length(units)}. For matched units, each entry is a
#'   data frame (if \code{multiple = FALSE}) or a list of data frames (if
#'   \code{multiple = TRUE}). For a given entry, these data frames are subsets
#'   of \code{data} passed to \code{\link{FLAME}} to generate \code{FLAME_out},
#'   whose rows correspond to the units in the matched group(s) of that entry
#'   and whose columns are the respective values of treatment, outcome, and
#'   covariates matched on. For unmatched units, the relevant list entries are
#'   \code{NULL}.
#'
#' @export
MG <- function(units, AME_out, multiple = FALSE, index_only = FALSE) {

  if (any(!(units %in% as.numeric(rownames(AME_out$data))))) {
    stop('Supplied a unit not in `AME_out$data`.')
  }

  if (index_only) {
    return(AME_out$MGs[units])
  }

  col_names <- colnames(AME_out$data)

  cov_inds <-
    which(!(col_names %in%
              c('treated', 'outcome', 'weight', 'matched', 'original_ind')))

  all_MGs <- vector('list', length(units))
  all_MGs <- lapply(all_MGs, function(x) x <- vector('list', 1))

  for (i in seq_along(units)) {
    if (is.null(AME_out$MGs[[units[i]]])) {
      next
    }

    if (!multiple) {
      MGs_to_return <- units[i]
    }
    else {
      in_MG <- which(sapply(AME_out$MGs, function(x) units[i] %in% x))
      MGs_to_return <- in_MG[!duplicated(AME_out$MGs[in_MG])]
    }

    for (j in seq_along(MGs_to_return)) {
      MG <- AME_out$data[AME_out$MGs[[MGs_to_return[j]]], ]
      matched_inds <-
        which(sapply(MG[, cov_inds], function(col) !('*' %in% col)))

      keep_inds <-
        c(matched_inds, which(col_names %in% c('treated', 'outcome')))
      all_MGs[[i]][[j]] <- MG[, keep_inds]
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
#' This function returns a single CATE estimate for a subgroup of interest, as
#' defined by the units in \code{units}. This is computed by averaging the CATE
#' estimates of all specified units, weighted by the weights of their main
#' matched groups. The weight of a main matched group is defined as the sum of
#' the inverse weights of the constituent units. Unmatched units are omitted
#' from the computation. To retrieve CATE estimates for individual units or
#' matched groups, use \code{AME_out$CATE}.
#'
#' @seealso \code{\link{FLAME}}, \code{\link{DAME}}
#'
#' @param units A vector of indices for the units whose CATEs
#'   are desired.
#' @param AME_out An object of class \code{AME}.
#'
#' @return A scalar numeric with the CATE estimate for \code{units}.
#'
#' @examples
#' data <- gen_data()
#' holdout <- gen_data()
#' FLAME_out <- FLAME(data = data, holdout = holdout)
#' CATE(which(data$X1 == 1), FLAME_out)
#' @export
CATE <- function(units, AME_out) {
  if (!('outcome' %in% colnames(AME_out$data))) {
    stop('`AME_out$data` does not contain an outcome; cannot estimate CATE.')
  }

  if (is.factor(AME_out$data$outcome)) {
    stop("Cannot estimate treatment effects with a categorical outcome.")
  }

  if (any(!(units %in% as.numeric(rownames(AME_out$data))))) {
    stop('Supplied a unit not in `AME_out$data`.')
  }

  # Filter out unmatched units
  units <- units[AME_out$data$matched[units]]
  if (length(units) == 0) {
    stop('None of the supplied units were matched; a CATE cannot be estimated.')
  }

  weights <- sapply(AME_out$MGs[units], function(x) sum(AME_out$data$weight[x]))
  CATEs <- AME_out$CATE[units]

  return(weighted.mean(CATEs, weights))
}

#' ATE of a matched dataset
#'
#' \code{ATE} estimates the average treatment effect (ATE) of a matched dataset.
#'
#' The ATE is estimated as the CATE estimate across all units in the data.
#' @seealso \code{\link{CATE}}, \code{\link{ATC}}, \code{\link{ATT}}
#' @param AME_out An object of class \code{AME}.
#' @export
ATE <- function(AME_out) {

  if (!('outcome' %in% colnames(AME_out$data))) {
    stop('`AME_out$data` does not contain an outcome; cannot estimate ATE.')
  }

  if (is.factor(AME_out$data$outcome)) {
    stop("Cannot estimate treatment effects with a categorical outcome.")
  }

  weights <- sapply(1:nrow(AME_out$data), function(x) {
    if (is.null(AME_out$MGs[x])) {
      return(0)
    }
    return(sum(1 / AME_out$data$weight[AME_out$MGs[[x]]]))
  })

  return(weighted.mean(AME_out$CATE, weights))
}

#' ATT of a matched dataset
#'
#' \code{ATT} estimates the average treatment effect on the treated (ATT) of a
#' matched dataset.
#'
#' The ATT is estimated as the CATE estimate across all treated units in the
#' data.
#' @seealso \code{\link{CATE}}, \code{\link{ATE}}, \code{\link{ATC}}
#' @param AME_out An object of class \code{AME}.
#' @export
ATT <- function(AME_out) {
  return(CATE(which(AME_out$data$treated == 0), AME_out))
  # if (!('outcome' %in% colnames(AME_out$data))) {
  #   stop('`AME_out$data` does not contain an outcome; cannot estimate ATT.')
  # }
  #
  # if (is.factor(AME_out$data$outcome)) {
  #   stop("Cannot estimate treatment effects with a categorical outcome.")
  # }
  #
  # weights <- sapply(1:nrow(AME_out$data), function(x) {
  #   if (is.null(AME_out$MGs[x]) | AME_out$data$treated[x] == 0) {
  #     return(0)
  #   }
  #   controls <- intersect(AME_out$MGs[[x]], AME_out$data$treated == 0)
  #   return(sum(1 / AME_out$data$weight[controls]))
  # })
  #
  # all_controls <- which(AME_out$data$treated == 0)
  # sum_weights <- 0
  # weighted_sum_CATEs <- 0
  #
  # for (i in 1:nrow(AME_out$data)) {
  #   if (AME_out$data$treated[i] == 0 || !AME_out$data$matched[i]) {
  #     next
  #   }
  #   MMG <- AME_out$MGs[[i]]
  #
  #   controls <- intersect(MMG, all_controls)
  #
  #   weight <- sum(1 / AME_out$data$weight[controls])
  #   sum_weights <- sum_weights + weight
  #
  #   weighted_sum_CATEs <-
  #     weighted_sum_CATEs +
  #     weight * (AME_out$data$outcome[i] - mean(AME_out$data$outcome[controls]))
  # }
  # return(weighted_sum_CATEs / sum_weights)
}

#' ATC of a matched dataset
#'
#' \code{ATC} estimates the average treatment effect on the controls (ATC) of a
#' matched dataset.
#'
#' The ATC is estimated as the CATE estimate across all control units in the
#' data.
#' @seealso \code{\link{CATE}}, \code{\link{ATE}}, \code{\link{ATT}}
#' @param AME_out An object of class \code{AME}.
#' @export
ATC <- function(AME_out) {
  return(CATE(which(AME_out$data$treated == 0), AME_out))
  # treated <- which(AME_out$data$treated == 1)
  # control <- which(AME_out$data$treated == 0)
  # AME_out$data$treated[treated] <- 0
  # AME_out$data$treated[control] <- 1
  # return(-ATT(AME_out))
}
