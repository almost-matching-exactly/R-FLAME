gen_new_active_sets <- function(s, delta) {
  k <- length(s)
  Z <- list()

  # 3
  if (length(delta) == 0) {
    return (Z)
  }

  delta_k <- list()
  counter <- 1
  for (i in seq_along(delta)) {
    if (length(delta[[i]]) == k) {
      delta_k[[counter]] <- delta[[i]]
      counter <- counter + 1
    }
  }
  delta_k[[counter]] <- s

  # 4, 5
  supp <- table(unlist(delta_k))

  # 6
  omega <- setdiff(strtoi(names(supp[supp >= k])), s)

  # 7
  counter <- 1
  if (all(supp[match(s, names(supp))] >= k)) {
    # 8
    for (a in omega) {
      # 9
      # Necessary to sort?
      r <- sort(c(s, a))
      # 10
      if (all(combn(r, k, simplify = FALSE) %in% delta_k)) {
        # 11
        Z[[counter]] <- r
        counter <- counter + 1
      }
    }
  }
  return(Z)
}
