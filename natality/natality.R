# library(readr)
# library(dplyr)
# library(magrittr)
#
# replace_missing_indicator <- function(col, missing_indicator) {
#   col[col == missing_indicator] <- NA
#   return(col)
# }
#
# df <- read.csv('~/Downloads/natl2010.csv', header = TRUE, fileEncoding = 'ASCII', nrows = 500000)
#
# df %<>%
#   filter(dplural == 1, gestrec3 == 2, dbwt != 9999,
#          cig_1 != 99, cig_2 != 99, cig_3 != 99) %>%
#   mutate(smokes10 = (cig_1 >= 10 & cig_2 >= 10 & cig_3 >= 10),
#          nonsmoker = (cig_1 == 0 & cig_2 == 0 & cig_3 == 0)) %>%
#   filter(smokes10 == 1 | nonsmoker == 1) %>%
#   select(c('fagerec11', 'meduc', 'dob_wk', 'mager9',
#          'fbrace', 'mar','rf_diab', 'rf_cesar',
#          'tbo_rec', 'rf_phyp', 'urf_chyper', 'lbo_rec',
#          'rf_ppterm', 'sex', 'mbrace', 'bfacil3', 'smokes10', 'dbwt')) %>%
#   mutate(bfacil3 = replace_missing_indicator(bfacil3, 9),
#          mar = replace_missing_indicator(mar, 9),
#          meduc = replace_missing_indicator(meduc, 9),
#          fagerec11 = replace_missing_indicator(fagerec11, 11),
#          fbrace = replace_missing_indicator(fbrace, 99),
#          lbo_rec = replace_missing_indicator(lbo_rec, 9),
#          tbo_rec = replace_missing_indicator(tbo_rec, 9),
#          urf_chyper = replace_missing_indicator(urf_chyper, 9),
#          rf_diab = replace_missing_indicator(rf_diab, 'U'),
#          rf_phyp = replace_missing_indicator(rf_phyp, 'U'),
#          rf_ppterm = replace_missing_indicator(rf_ppterm, 'U'),
#          rf_cesar = replace_missing_indicator(rf_cesar, 'U'))
#
# outcome <- 'dbwt'
# treatment <- 'smokes10'
# df[, !(colnames(df) %in% c(outcome, treatment))] %<>% lapply(as.factor)
#
# # write.csv(df, file = 'natl2010_cleaned.csv')
# my_PE <- function(X, Y) {
#   df <- as.data.frame(cbind(X, Y = Y))
#   return(lm(Y ~ ., df)$fitted.values)
# }
#
# system.time(flout <- FLAME(df, holdout = 0.25, treated_column_name = 'smokes10',
#                outcome_column_name = 'dbwt', verbose = 3,
#                missing_holdout = 'drop', missing_data = 'drop',
#                estimate_CATEs = TRUE,
#                PE_method = my_PE))


system.time(m_nn_pscore_10 <- matchit(smokes10 ~ fagerec11 + meduc + dob_wk + mager9 + fbrace +
                    mar + rf_diab + rf_cesar + tbo_rec + rf_phyp + urf_chyper +
                    lbo_rec + rf_ppterm + sex + mbrace + bfacil3, data = df,
                  method = "nearest", distance = "glm", ratio = 10))

system.time(m_cem <- matchit(smokes10 ~ fagerec11 + meduc + dob_wk + mager9 + fbrace +
                   mar + rf_diab + rf_cesar + tbo_rec + rf_phyp + urf_chyper +
                   lbo_rec + rf_ppterm + sex + mbrace + bfacil3, data = df,
                 method = 'cem'))


# saveRDS(flout, 'natality_out_500k_lm.rds')
#
# #
# # MG_size <- sapply(flout$MGs, length)
# # matched <- MG_size > 0
# # plot(MG_size[matched], flout$CATE[matched], xlim = c(0, 1000))
# # abline(h = 0, col = 'red')
