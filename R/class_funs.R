print.ame <- function(x, ...) {

  df <- x$data
  n_matches <- c(table(df$treated, df$matched))

  n_total <- c(sum(n_matches[c(1, 3)]), sum(n_matches[c(2, 4)]))

  cat(paste0('Sample Sizes:\n'))
  cat(format('', width = 10), format('Control Treated'))
  cat('\n')
  cat(format('All', width = 10), format(n_total, width = 7))
  cat('\n')
  cat(format('Matched', width = 10), format(n_matches[3:4], width = 7))
  cat('\n')
  cat(format('Unmatched', width = 10), format(n_matches[1:2], width = 7))
  cat('\n')


  TEs <- c(ATE = ATE(x), ATC = ATC(x), ATT = ATT(x))
  cat(format(c('ATE', 'ATC', 'ATT'), width = 5, justify = 'c'))
  cat('\n')
  cat(format(TEs, digits = 2, width = 5, justify = 'c'))
  cat('\n')


  cat('Matched Groups:\n')
  matched <- !sapply(x$MGs, is.null)

  MGs <- x$MGs[!duplicated(x$MGs) & matched]
  median_MG_size <- median(sapply(MGs, length))
  cat('  Median size:', median_MG_size)
  cat('\n')

  outcome_name <- 'outcome'
  treated_name <- 'treated'

  cov_inds <- !(colnames(df) %in% c(outcome_name, treated_name,
                                    'original_ind', 'matched', 'weight'))

  covs <- df[df$matched, cov_inds]
  n_covs_matched_on <- apply(covs, 1, function(x) sum(x != '*'))

  cat('  Median number of covariates matched on:',
      median(n_covs_matched_on))

  return(invisible(x))
}

plot.ame <- function(x, ...) {
  plot(sapply(x$MGs, length), x$CATE,
       xlab = 'Size of Matched Group',
       ylab = 'Estimated Conditional Average Treatment Effect')

  # Could also consider plotting ATE
  abline(h = 0, lty = 2)

  plot(density(x$CATE, na.rm = TRUE),
       xlab = c('Estimated Conditional Average Treatment Effect'),
       ylab = '',
       main = '')
  abline(v = ATE(flout), lty = 2)
  legend('topright', legend = c('Estimated ATE'), lty = 2)
}
