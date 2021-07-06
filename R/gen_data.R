#' Generate Toy Data for Matching
#'
#' \code{gen_data} generates toy data that can be used to explore FLAME and DAME
#' functionality.
#'
#' \code{gen_data} simulates data in the format accepted by \code{\link{FLAME}}
#' and \code{link{DAME}}. Covariates \eqn{X_i} and treatment \eqn{T} are both
#' independently generated according to a Bernoulli(0.5) distribution. The
#' outcome \eqn{Y} is generated according to \eqn{Y = 15X_1 - 10X_2 + 5X_3 + 5T
#' + \epsilon}, where \eqn{\epsilon \sim N(0, I_n)}. Thus, the value of \code{p}
#' must be at least 3 and any additional covariates beyond \eqn{X_1, X_2, X_3}
#' are irrelevant.
#'
#' @param n Number of units desired in the data set. Defaults to 250.
#' @param p Number of covariates in the data set. Must be greater than 2.
#'   Defaults to 5.
#' @param write A logical scalar. If \code{TRUE}, the resulting data is stored
#'   as a .csv file as specified by arguments \code{path} and \code{filename}.
#'   Defaults to \code{FALSE}.
#' @param path The path to the location where the data should be written if
#'   \code{write = TRUE}. Defaults to \code{getwd()}.
#' @param filename The name of the file to which the data should be written if
#'   \code{write = TRUE}. Defaults to AME.csv.
#'
#' @return A data frame that may be passed to \code{\link{FLAME}} or
#'   \code{\link{DAME}}. Covariates are numeric, treatment is binary numeric and
#'   outcome is numeric.
#' @export
gen_data <- function(n = 250, p = 5,
                     write = FALSE, path = getwd(), filename = 'AME.csv') {
  if (p <= 2) {
    stop('`p` must be greater than 2')
  }
  TE <- 5

  covs <- matrix(rbinom(n * p, 1, prob = 0.5), nrow = n)

  covs <-
    sample(1:4, size = n * p, replace = TRUE, prob = c(0.2, 0.3, 0.4, 0.1))
  covs <- matrix(covs, nrow = n)

  treated <- rbinom(n, 1, prob = 0.5)

  outcome <-
    15 * covs[, 1] - 10 * covs[, 2] + 5 * covs[, 3] +
    TE * treated +
    rnorm(n)

  data <- data.frame(covs, outcome = outcome, treated = treated)
  if (write) {
    write.csv(data, file = paste0(path, '/', filename),
              row.names = FALSE)
  }
  return(data)
}

gen_data2 <- function(n = 250, p_miss = 0) {
  # Data generation with non-numeric factor names to test that functionality
  p <- 5
  data <- data.frame(Gender = sample(c('M', 'F'), size = n, replace = TRUE),
                        Race = sample(c('Black', 'Asian',
                                        'White', 'Native American'),
                                      size = n, replace = TRUE,
                                      prob = c(.15, .10, .7, .05)),
                        Latino = sample(c('Yes', 'No', 'Not Recorded'),
                                        size = n, replace = TRUE,
                                        prob = c(0.1, 0.88, 0.02)),
                        NumPriors = sample(0:4,
                                           size = n, replace = TRUE,
                                           prob = c(0.6, 0.2, 0.1, 0.05, 0.05)),
                        Education = sample(c('Below HS', 'HS', 'College',
                                             'Grad', 'Not Recorded'),
                                           size = n, replace = TRUE,
                                           prob = c(0.05, 0.3, 0.5, 0.1, 0.05)),
                        treated = sample(c(TRUE, FALSE),
                                         size = n, replace = TRUE))
  data$outcome <-
    3 * (data$Race %in% c('White', 'Black')) +
    1 * (data$Latino == 'Yes') +
    2 * (data$NumPriors == 0) + 2 * (data$NumPriors > 2) +
    5.5 * data$treated +
    3 * (data$Education %in% c('HS', 'College')) +
    rnorm(n)

  if (p_miss > 0) {
    for (j in 1:p) {
      data[[j]][sample(1:n, size = round(p_miss * n))] <- NA
    }
  }
  return(data)
}

gen_missing_data <- function(n = 250, p = 3,
                             write = FALSE, percent_missing = 0.05) {
  # Data generation to test missing data functionality
  data <- gen_data(n = n, p = p, write = write, filename = 'missing_data.csv')
  covs <- as.matrix(data[, 1:p])
  size <- n * p
  inds <- sample(1:size, size = round(percent_missing * size), replace = FALSE)
  covs[inds] <- NA
  data[, 1:p] <- covs
  data[, seq_len(ncol(covs))] <- lapply(data[, seq_len(ncol(covs))], as.factor)
  return(data)
}
