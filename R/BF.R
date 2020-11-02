get_BF <- function(cov_to_drop, data, replace, covs) {

  # Calculate number of units eligible to be matched
  if (replace) {
    n_control <- sum(data$treated[!data$missing] == 0)
    n_treated <- sum(data$treated[!data$missing] == 1)
  }
  else {
    n_control <- sum(data$treated[!data$matched & !data$missing] == 0)
    n_treated <- sum(data$treated[!data$matched & !data$missing] == 1)
  }

  match_out <- exact_match_bit(data, setdiff(covs, cov_to_drop), replace)
  newly_matched <- match_out$newly_matched

  # Newly matched
  n_control_matched <- sum(data$treated[newly_matched] == 0)
  n_treated_matched <- sum(data$treated[newly_matched] == 1)

  # All matched units; for stopping rule purposes
  all_unmatched <-
    setdiff(1:nrow(data), union(newly_matched, which(data$matched)))

  n_control_unmatched <- sum(all_unmatched %in% which(data$treated == 0))
  n_treated_unmatched <- sum(all_unmatched %in% which(data$treated == 1))

  prop_control_unmatched <- n_control_unmatched / n_control
  prop_treated_unmatched <- n_treated_unmatched / n_treated

  # Is this if_else really necessary? We should catch it earlier
  if (n_control == 0 | n_treated == 0) {
    BF <- 0
  }
  else {
    BF <- n_control_matched / n_control + n_treated_matched / n_treated
  }

  return(list(BF = BF,
              prop_unmatched =
                list(control = prop_control_unmatched,
                     treated = prop_treated_unmatched)))
}
