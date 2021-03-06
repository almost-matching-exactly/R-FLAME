% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/post_matching.R
\name{MG}
\alias{MG}
\title{Matched Groups}
\usage{
MG(units, ame_out, multiple = FALSE, id_only = FALSE, index_only)
}
\arguments{
\item{units}{A vector of units whose matched groups are desired.}

\item{ame_out}{An object of class \code{ame}.}

\item{multiple}{A logical scalar. If \code{FALSE} (default), then \code{MG}
will only return the main matched group for each unit. See below for
details. Cannot be set to \code{TRUE} if \code{ame_out} was generated
without replacement.}

\item{id_only}{A logical scalar. If \code{TRUE}, then only the IDs of the
units in each matched group are returned, and not their treatment, outcome,
or covariate information.}

\item{index_only}{Defunct. Use `id_only` instead.}
}
\value{
A list of length \code{length(units)}, each entry of which corresponds to a
  different unit in \code{units}. For matched units, if \code{multiple =
  FALSE}, each entry is 1. a data frame containing the treatment and outcome
  information of members of the matched group, along with covariates they
  were matched on if \code{id_only = FALSE} or 2. a vector of the IDs of
  matched units if \code{id_only = TRUE} . If \code{multiple = TRUE}, each
  entry of the returned list is a list containing the previously described
  information, but with each entry corresponding to a different matched
  group. In either case, entries corresponding to unmatched units are
  \code{NULL}.
}
\description{
\code{MG} returns the matched groups of the supplied units.
}
\details{
The \code{units} argument refers to units with respect to
\code{rownames(ame_out$data)}. Typically, this will also correspond to the
indexing of the data (i.e. passing \code{units = 3} will return the matched
group of the 3rd unit in the matching data). However, if a separate holdout
set was not passed to the matching algorithm or if the original matching data
had rownames other than \code{1:nrow(data)}, then this is not the case.

The \code{multiple} argument toggles whether only a unit's main matched group
(MMG) or all matched groups a unit is part of should be returned. A unit's
MMG contains its highest quality matches (that is, the units with which it
first matched in the sequence of considered covariate sets). If the original
call that generated \code{ame_out} specified \code{replace = FALSE} then
units only are part of one matched group (which is also their MMG) and
\code{multiple} must be set to \code{FALSE}.
}
\examples{
\dontrun{
data <- gen_data()
holdout <- gen_data()
FLAME_out <- FLAME(data = data, holdout = holdout, replace = TRUE)

# Only the main matched group of unit 1
MG(1, FLAME_out, multiple = F)

# All matched groups of unit 1
MG(1, FLAME_out, multiple = T)
}
}
