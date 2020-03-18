#' @export
gen_data <- function(n = 50, p = 5, write = FALSE, path = getwd(), filename = 'FLAME.csv') {
  TE <- 5
  # covs <-
  #   rbinom(n * p, 1, prob = 0.5) %>%
  #   matrix(nrow = n)
  covs <-
    apply(rmultinom(n * p, size = 1, prob = c(0.3, 0.2, 0.4, 0.1)) == 1,
          2,
          which) %>%
    matrix(nrow = n)
  treated <- rbinom(n, 1, prob = 0.5)

  # outcome <-
  #   1 * (covs[, 1] == 1 | ((covs[, 2] + covs[, 3] + treated * covs[, 4]) >= 2))

  outcome <-
    (15 * covs[, 1] - 10 * covs[, 2] + 5 * covs[, 3] - 2.5 * covs[, 4]) %>%
    magrittr::add(rnorm(n)) %>%
    magrittr::add(TE * treated)
  # outcome <-
  #   (15 * (covs[, 1] == 3) - 10 * (covs[, 2] == 1) + 5 * (covs[, 3] == 2)) %>%
  #   magrittr::add(rnorm(n)) %>%
  #   magrittr::add(TE * treated)
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
          rnorm(n, -10, 1))


  treated <- rbinom(n, 1, prob = 0.5)

  outcome <-
    outcome <- 10 * covs[, 1] + 5 * covs[, 2] %>%
    # (15 * covs[, 7] - 10 * covs[, 2] + 5 * covs[, 3] - 2.5 * covs[, 4] + 3 * covs[, 6]) %>%
    magrittr::add(rnorm(n)) %>%
    magrittr::add(TE * treated)
  data <- data.frame(covs, outcome = outcome, treated = treated)
  # data[, 1:(ncol(covs) - 3)] %<>% lapply(as.factor)
  return(data)
}

gen_missing_data <- function(n = 50, p = 5, write = FALSE, percent_missing = 0.05) {
  data <- gen_data(n = n, p = p, write = write)
  covs <- as.matrix(data[, 1:p])
  size <- n * p
  inds <- sample(1:size, size = round(percent_missing * size), replace = FALSE)
  covs[inds] <- NA
  data[, 1:p] <- covs
  data[, 1:ncol(covs)] %<>% lapply(as.factor)
  return(data)
}
