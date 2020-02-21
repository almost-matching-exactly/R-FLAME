gen_data <- function(n, p) {
  covs <-
    rbinom(n * p, 1, prob = 0.5) %>%
    matrix(nrow = n)
  treated <- rbinom(n, 1, prob = 0.5)
  outcome <-
    3 * (covs[, 1] + covs[, 2]) - 5 * covs[, 3] %>%
    add(rnorm(n))
  data <- data.frame(covs, outcome = outcome, treated = treated)
  return(data)
}
