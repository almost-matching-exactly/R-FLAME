
# n <- c(500, 1000, 2500, 5000, 10000)
# p <- 30
# new_times <- NULL
# for (i in seq_along(n)) {
#   reps <- matrix(ncol = length(n), nrow = 2)
#   for (j in 1:5) {
#     data <- gen_data(n[i], p)
#     holdout <- gen_data(n[i], p)
#     reps[1, j] <-
#       system.time(flout <- FLAME_bit(data = data,
#                                      holdout = holdout,
#                                      verbose = 0, opt = 0))[3]
#     reps[2, j] <-
#       system.time(flout_new <- FLAME_bit_new(data = data,
#                                              holdout = holdout,
#                                              verbose = 0, opt = 1))[3]
#     if (!all(flout$CATE == flout_new$CATE)) {
#       stop('Unequal methods')
#     }
#   }
#   message('new n value')
#   mean_times <- rowMeans(reps)
#   sd_times <- apply(reps, 1, sd)
#   new_times %<>% rbind(cbind(rep(n[i], 2), mean_times, sd_times))
# }
# new_times %<>%
#   cbind(p) %>%
#   as.data.frame() %>%
#   `colnames<-`(c('n_obs', 'avg', 'sd', 'p')) %>%
#   mutate(estimator = rep(c('Old', 'New'), length(n)),
#          ymin = avg - sd, ymax = avg + sd)
#
# # times %<>%
# #   cbind(p) %>%
# #   as.data.frame() %>%
# #   `colnames<-`(c('n_obs', 'avg', 'sd', 'p')) %>%
# #   mutate(estimator = rep(c('Old', 'New'), length(n)),
# #          ymin = avg - sd, ymax = avg + sd)
#
# times %<>% rbind(new_times)
#
# ggplot(times, aes(x = n_obs, y = avg, ymin = ymin, ymax = ymax,
#                     color = factor(p), shape = factor(estimator))) +
#   geom_pointrange()
#
#
#
#
#
# # aggregate_table <- function(list) {
# #   tab = table(as.character(list))
# #   tab = unclass(tab)
# #   name = names(tab)
# #   list_val = as.character(list)
# #   return(as.vector(tab[sapply(list_val, function(x) which(name ==  x))]))
# # }
# #
# # opt1 <- function(vec) {
# #   tab = table(as.character(vec))
# #   tab = unclass(tab)
# #   name = names(tab)
# #   list_val = as.character(vec)
# #   return(as.vector(tab[match(vec, name)]))
# # }
# # #
# #
# #
# # #
# # bu <- sample(1:100, 5000, replace = T)
# # #
# # microbenchmark(aggregate_table(bu), opt1(bu), check = 'identical')
# # microbenchmark({table(bu);as.numeric(names(bu))})
