check_all_matches <- function(flout, data) {
  max_covs <- length(flout$matched_on[[1]])
  min_covs <- length(flout$matched_on[[length(flout$matched_on)]])
  matched <- NULL
  for (n_covs in seq(max_covs, min_covs)) {
    message('n_covs = ', n_covs)
    check_matches(flout, n_covs, data, matched)
    MG_exclude <- which(sapply(flout$matched_on, function(x) length(x) >= n_covs))
    matched <- unlist(flout$MGs[MG_exclude])
  }
  return(TRUE)
}
check_matches <- function(flout, n_covs, data, matched) {
  MGs <- flout$MGs
  matched_on <- flout$matched_on
  to_eliminate <- NULL
  for (i in 1:length(matched_on)) {
    if (length(matched_on[[i]]) != n_covs) {
      to_eliminate <- c(to_eliminate, i)
    }
    else {
      covs <-
          match(names(matched_on[[i]]), colnames(data)[1:(length(colnames(data)) - 2)])
    }
  }
  # browser()
  MGs[to_eliminate] <- NULL
  for (i in 1:nrow(data)) { # For each unit
    if (i %in% matched) {
      next
    }
    ever_matched <- FALSE # Have we found them to be matched?
    matches <- NULL # Who have we found them to be matched with?
    to_match <- data[i, covs] # Covariate values we're matching to
    for (j in 1:nrow(data)) { # For everyone else in the data
      if (j == i | j %in% matched) { # If they've already been matched, carry on
        next
      }
      if (all(data[j, covs] == to_match)) { # If unit j matches unit i
        ever_matched <- TRUE # i has been matched
        matches <- c(matches, j) # Add j to the list of units i has been matched to
        # matched <- c(matched, j) # Add j to the list of units that has ever been matchedd
      }
    }
    # Now we have a list, matches, of all the people i matched to

    if (length(unique(union(data$treated[i], data$treated[matches]))) != 2) {
      ever_matched <- FALSE
    }

    if (!ever_matched) {
      for (k in 1:length(MGs)) {
        if (i %in% MGs[[k]]) { # If i found
          print(i)
          print(matches)
          stop('uh oh A')
          break
        }
      }
      next # if no discrepancy, carry on
    }

    # Valid match
    matches <- c(matches, i) # add him to his own MG
    # matched <- c(matched, i) # add him to the list of matched people

    # if we matched i
    ever_found <- FALSE # did FLAME ever match i
    for (k in 1:length(MGs)) {
      if (i %in% MGs[[k]]) { #if FLAME did match i
        ever_found <- TRUE
        if (!all(sort(MGs[[k]]) == sort(matches))) { # If MGs don't agree
          ever_found <- FALSE # discrepancy
          print(i)
          print(matches)
          stop('uh oh C')
          break
        }
        else {
          # If they matched on all the same covariates np
        }
      }
    }
    if (!ever_found) { # If discrepancy, cut
      print(i)
      print(matches)
      stop('uh oh D')
    }
  }
  return(TRUE)
}

# data <- gen_data()
# holdout <- gen_data()
# treatment_column_name <- 'treated'
# outcome_column_name <- 'outcome'
# out <- organize_data(data, holdout, treatment_column_name, outcome_column_name)
#
# f_a <- function(x) {
#   c(data, holdout, covs, n_covs, n_levels, cov_names, sorting_order) %<-% x
#   list(data, holdout, covs, n_covs, n_levels, cov_names, sorting_order)
# }
#
# f_b <- function(out) {
#   data <- out[[1]]
#   holdout <- out[[2]]
#   covs <- out[[3]]
#   n_covs <- out[[4]]
#   n_levels <- out[[5]]
#   cov_names <- out[[6]]
#   sorting_order <- out[[7]]
#   list(data, holdout, covs, n_covs, n_levels, cov_names, sorting_order)
# }

# eta <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5)
# max_depth <- c(2, 3, 4, 6, 8)
# alpha <- c(0.01, 0.1, 0.5, 1, 5)
# nrounds <- c(5, 10, 50, 100, 200)
# subsample <- c(0.1, 0.3, 0.5, 0.75, 1)
#
# f_a <- function() {
#   expand.grid(eta, max_depth, alpha, nrounds, subsample) %>%
#     `colnames<-`(c('eta', 'max_depth', 'alpha', 'nrounds', 'subsample')) %>%
#     return()
# }
#
# f_b <- function() {
#   grd <- expand.grid(eta, max_depth, alpha, nrounds, subsample)
#   colnames(grd) <- c('eta', 'max_depth', 'alpha', 'nrounds', 'subsample')
#   return(grd)
# }
#
# microbenchmark::microbenchmark(f_a(), f_b(), check = 'identical')


# microbenchmark::microbenchmark()
# n <- c(500, 1000, 2500, 5000, 10000, 20000)
# p <- c(25, 30, 35, 40)
# n_reps <- 3
#
# avg_times <- array(dim = c(length(p), length(n), 2))
# sd_times <- array(dim = c(length(p), length(n), 2))
#
# for (i in seq_along(n)) {
#   for (j in seq_along(p)) {
#     data <- gen_data(n[i], p[j])
#     holdout <- gen_data(n[i], p[j])
#     reps <- matrix(nrow = 2, ncol = n_reps)
#     for (k in 1:n_reps) {
#       reps[1, k] <-
#         system.time(flout <- FLAME_bit(data = data,
#                                        holdout = holdout,
#                                        verbose = 0, opt = 0))[3]
#       reps[2, k] <-
#         system.time(flout <- FLAME_bit(data = data,
#                                        holdout = holdout,
#                                        verbose = 0, opt = 1))[3]
#     }
#     avg_times[j, i, ] <- rowMeans(reps)
#     sd_times[j, i, ] <- apply(reps, 1, sd)
#     print(avg_times[j, i, ])
#     saveRDS(avg_times, 'avgs.rds')
#     saveRDS(sd_times, 'sd.rds')
#     message(paste0('Done with n = ', n[i], ' and p = ', p[j]))
#   }
# }
#
# saved_times <- readRDS('avgs.rds')
# saved_times <-
#   rbind(saved_times[,,1], saved_times[,,2]) %>%
#   as.data.frame() %>%
#   dplyr::mutate(FLAME = rep(c('Old', 'New'), each = nrow(.) / 2))
#
# ggdata <-
#   reshape2::melt(saved_times, id.vars = 'FLAME') %>%
#   dplyr::mutate(p = rep(rep(p, 2), 6),
#          n = rep(n, each = 8))
#
# saved_sd_times <- readRDS('sd.rds')
# saved_sd_times <-
#   rbind(saved_sd_times[,,1], saved_sd_times[,,2]) %>%
#   as.data.frame() %>%
#   dplyr::mutate(FLAME = rep(c('Old', 'New'), each = nrow(.) / 2)) %>%
#   reshape2::melt(id.vars = 'FLAME') %>%
#   dplyr::mutate(p = rep(rep(p, 2), 6),
#                 n = rep(n, each = 8))
#
# ggdata %<>%
#   dplyr::mutate(sd = saved_sd_times$value)
#
# ggdata %<>%
#   dplyr::mutate(ymin = value - sd,
#                 ymax = value + sd)
#
# ggplot(ggdata, aes(color = FLAME, x = n, y = value, ymin = ymin, ymax = ymax)) +
#   geom_point() +
#   geom_line() +
#   geom_pointrange() +
#   facet_grid(. ~ p) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   labs(y = 'Runtime (sec)',
#        title = 'Runtime Scaling as a Function of p')
#
# # new_times <- NULL
# # for (i in seq_along(n)) {
# #   reps <- matrix(ncol = 3, nrow = 2)
# #   for (j in 1:3) {
# #     data <- gen_data(n[i], p)
# #     holdout <- gen_data(n[i], p)
# #     reps[1, j] <-
# #       system.time(flout <- FLAME_bit(data = data,
# #                                      holdout = holdout,
# #                                      verbose = 0, opt = 0))[3]
# #     reps[2, j] <-
# #       system.time(flout_new <- FLAME_bit_new(data = data,
# #                                              holdout = holdout,
# #                                              verbose = 0, opt = 1))[3]
# # if (!all(flout$CATE == flout_new$CATE)) {
# #   stop('Unequal methods')
# # }
# #   }
# #   message('new n value')
# #   mean_times <- rowMeans(reps)
# #   sd_times <- apply(reps, 1, sd)
# #   new_times %<>% rbind(cbind(rep(n[i], 2), mean_times, sd_times))
# # }
# # new_times %<>%
# #   cbind(p) %>%
# #   as.data.frame() %>%
# #   `colnames<-`(c('n_obs', 'avg', 'sd', 'p')) %>%
# #   mutate(estimator = rep(c('Old', 'New'), length(n)),
# #          ymin = avg - sd, ymax = avg + sd)
# #
# # # times %<>%
# # #   cbind(p) %>%
# # #   as.data.frame() %>%
# # #   `colnames<-`(c('n_obs', 'avg', 'sd', 'p')) %>%
# # #   mutate(estimator = rep(c('Old', 'New'), length(n)),
# # #          ymin = avg - sd, ymax = avg + sd)
# #
# # times %<>% rbind(new_times)
# #
# # ggplot(times, aes(x = n_obs, y = avg, ymin = ymin, ymax = ymax,
# #                     color = factor(p), shape = factor(estimator))) +
# #   geom_pointrange()
# #
# #
# #
# #
# #
# # # aggregate_table <- function(list) {
# # #   tab = table(as.character(list))
# # #   tab = unclass(tab)
# # #   name = names(tab)
# # #   list_val = as.character(list)
# # #   return(as.vector(tab[sapply(list_val, function(x) which(name ==  x))]))
# # # }
# # #
# # # opt1 <- function(vec) {
# # #   tab = table(as.character(vec))
# # #   tab = unclass(tab)
# # #   name = names(tab)
# # #   list_val = as.character(vec)
# # #   return(as.vector(tab[match(vec, name)]))
# # # }
# # # #
# # #
# # #
# # # #
# # # bu <- sample(1:100, 5000, replace = T)
# # # #
# # # microbenchmark(aggregate_table(bu), opt1(bu), check = 'identical')
# # # microbenchmark({table(bu);as.numeric(names(bu))})
