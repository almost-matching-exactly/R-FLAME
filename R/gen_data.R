#' Generate Toy Data for Matching
#'
#' \code{gen_data} generates toy data that can be used to explore FLAME's
#' functionality.
#'
#' \code{gen_data} simulates data in the format accepted by \code{FLAME}.
#'   Covariates \eqn{X_i} and treatment \eqn{T} are both independently generated according
#'   to a Bernoulli(0.5) distribution. The outcome \eqn{Y} is generated according to
#'   \eqn{Y = 15X_1 - 10X_2 + 5X_3 + 5T + \epsilon}, where
#'   \eqn{\epsilon \sim N(0, I_n)}. Thus, the value of
#'   \code{p} must be at least 3 and any additional covariates beyond
#'   \eqn{X_1, X_2, X_3} are irrelevant.
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
#'   \code{write = TRUE}. Defaults to FLAME.csv.
#' @export
gen_data <- function(n = 250, p = 5, write = FALSE, path = getwd(), filename = 'FLAME.csv') {
  if (p <= 2) {
    stop('p must be greater than 2')
  }
  TE <- 5
  covs <-
    rbinom(n * p, 1, prob = 0.5) %>%
    matrix(nrow = n)

  treated <- rbinom(n, 1, prob = 0.5)

  outcome <-
    (15 * covs[, 1] - 10 * covs[, 2] + 5 * covs[, 3]) %>%
    magrittr::add(rnorm(n)) %>%
    magrittr::add(TE * treated)

  data <- data.frame(covs, outcome = outcome, treated = treated)
  data[, 1:ncol(covs)] %<>% lapply(as.factor)
  if (write) {
    write.csv(data, file = paste0(path, '/', filename),
              row.names = FALSE)
  }
  return(data)
}

gen_mixed_data <- function(n = 200, p = 8) {
  TE <- 5
  # covs <-
  #   rbinom(n * (p - 3), 1, prob = 0.5) %>%
  #   matrix(nrow = n) %>%
  #   cbind(rnorm(n, 5, 7),
  #         rnorm(n, 0, 1),
  #         rnorm(n, -10, 1))

  covs <-
    cbind(rnorm(n, 5, 7),
          rnorm(n, 0, 1),
          rbinom(n, 1, 0.3),
          rbinom(n, 1, 0.7))


  treated <- rbinom(n, 1, prob = 0.5)

  outcome <-
    outcome <- 10 * covs[, 1] + 5 * covs[, 3] %>%
    # (15 * covs[, 7] - 10 * covs[, 2] + 5 * covs[, 3] - 2.5 * covs[, 4] + 3 * covs[, 6]) %>%
    magrittr::add(rnorm(n)) %>%
    magrittr::add(TE * treated)
  data <- data.frame(covs, outcome = outcome, treated = treated)
  # data[, 1:(ncol(covs) - 3)] %<>% lapply(as.factor)
  return(data)
}

gen_missing_data <- function(n = 250, p = 3, write = FALSE, percent_missing = 0.05) {
  data <- gen_data(n = n, p = p, write = write, filename = 'missing_data.csv')
  # data <- gen_mixed_data(n = n, p = p)

  covs <- as.matrix(data[, 1:p])
  size <- n * p
  inds <- sample(1:size, size = round(percent_missing * size), replace = FALSE)
  covs[inds] <- NA
  data[, 1:p] <- covs
  data[, 1:ncol(covs)] %<>% lapply(as.factor)
  return(data)
}
