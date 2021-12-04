# make_MGs_dev <- function(MGs, valid_matches, match_vals,
#                      matched, newly_matched, data, info) {
#
#   if (info$estimate_CATEs && info$outcome_type == 'continuous') {
#     Tr <- data$treated
#     Y <- data$outcome
#   }
#
#   # b_u values for those first matched on this cov set
#   newly_matched_vals <- match_vals[match(newly_matched, valid_matches)]
#
#   MG_ids <- match(as.character(newly_matched_vals),
#                   as.character(unique(newly_matched_vals)))
#
#   MG_counter <- max(data$MG)
#
#   for (i in seq_along(newly_matched)) {
#     new_MG <- valid_matches[match_vals == newly_matched_vals[i]]
#     MGs[[newly_matched[i]]] <- new_MG
#
#     if (info$estimate_CATEs && info$outcome_type == 'continuous') {
#       if (Tr[newly_matched[i]] == 1) {
#         data$CATE[newly_matched[i]] <-
#           Y[newly_matched[i]] - mean(Y[new_MG[Tr[new_MG] == 0]])
#       }
#       else {
#         data$CATE[newly_matched[i]] <-
#           mean(Y[new_MG[Tr[new_MG] == 1]]) - Y[newly_matched[i]]
#       }
#     }
#
#     # data$MG[new_MG] <- MG_counter + MG_ids[i]
#   }
#   return(list(MGs, data))
# }
#
# make_MGs_dev2 <- function(MGs, valid_matches, match_vals,
#                          matched, newly_matched, data, info) {
#
#   if (info$estimate_CATEs && info$outcome_type == 'continuous') {
#     Tr <- data$treated
#     Y <- data$outcome
#   }
#
#   # b_u values for those first matched on this cov set
#   newly_matched_vals <- match_vals[match(newly_matched, valid_matches)]
#
#   MG_ids <- match(as.character(newly_matched_vals),
#                   as.character(unique(newly_matched_vals)))
#
#   MG_counter <- max(data$MG)
#
#   oo <- outer(match_vals, newly_matched_vals, FUN = `==`)
#   browser()
#   MGs[newly_matched] <- apply(oo, 2, function(x) valid_matches[x])
#
#     # data$MG[new_MG] <- MG_counter + MG_ids[i]
#
#   return(MGs)
# }
#
# library(profvis)
# library(dplyr)
# n <- 5000
# p <- 10
#
# data <-
#   gen_data(n, p) %>%
#   mutate(matched = FALSE, missing = FALSE, MG = 0)
#
#
# replace <- FALSE
# covs <- 1:p
# MGs <- vector('list', n)
#
# info <- list(outcome_type = 'continuous',
#              treatment = 'treated',
#              algo = 'FLAME',
#              outcome = 'outcome',
#              replacement = replace,
#              estimate_CATEs = FALSE,
#              missing_data = 'none',
#              missing_holdout = 'none')
#
# match_out <- exact_match_bit(data, covs, replace)
# valid_matches <- match_out$valid_matches
# match_vals <- match_out$match_vals
# newly_matched <- match_out$newly_matched
# matched <- match_out$matched
#
# make_MGs_dev2(MGs, valid_matches, match_vals, matched, newly_matched, data, info)
#
# # profvis(make_MGs_dev(MGs, valid_matches, match_vals, matched, newly_matched, data, info))
#
# # microbenchmark::microbenchmark(make_MGs_dev(MGs, valid_matches, match_vals, matched, newly_matched, data, info),
# #                                make_MGs_dev2(MGs, valid_matches, match_vals, matched, newly_matched, data, info),
# #                                times = 10, check = 'equal')
#
# # newly_matched_vals <- match_vals[match(newly_matched, valid_matches)]
# # valid_matches[match_vals == newly_matched_vals[1]]
# # valid_matches[match_vals == newly_matched_vals[2]]
#
