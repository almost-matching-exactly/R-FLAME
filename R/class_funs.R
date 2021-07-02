n_covs_matched_on <- function(x, unit, cov_inds) {
  # Faster than for loop for moderate p
  if (is.null(x$MGs[[unit]])) {
    return(0)
  }

  data <- data.matrix(x$data[x$MGs[[unit]], cov_inds])

  # na.rm takes care of `missing_data = 'keep'`
  return(sum(colSums(sweep(data, 2, data[1, ]) ^ 2) == 0, na.rm = TRUE))
}

get_average_effects <- function(x) {

  unmatched <- !x$data$matched
  new_ind <- seq_along(unmatched) - cumsum(unmatched)

  x$MGs <- x$MGs[!unmatched]
  x$MGs <- lapply(x$MGs, function(z) new_ind[z])

  n <- length(x$MGs)
  Tr <- x$data[[x$info$treatment]][!unmatched]
  Y <- x$data[[x$info$outcome]][!unmatched]

  K <- numeric(n)

  for (i in 1:n) {
    Tr_i <- Tr[i]

    for (j in 1:n) {
      if (Tr[j] == Tr_i) {
        next
      }
      MG <- x$MGs[[j]]
      if (!(i %in% MG)) {
        next
      }
      opp_sign <- sum(Tr[MG] == Tr_i)
      K[i] <- K[i] + 1 / opp_sign
    }
  }

  n1 <- sum(Tr == 1)
  n0 <- sum(Tr == 0)

  ATE <- sum((2 * Tr - 1) * (1 + K) * Y) / n
  ATT <- sum((Tr - (1 - Tr) * K) * Y) / n1
  ATC <- sum((Tr * K - (1 - Tr)) * Y) / n0

  cond_var <- 0
  cond_var_t <- 0
  cond_var_c <- 0

  for (i in 1:n) {
    cond_var_tmp <- 0
    cond_var_tmp_t <- 0
    cond_var_tmp_c <- 0

    Tr_i <- Tr[i]
    Y_i <- Y[i]

    MG <- x$MGs[[i]]
    MG <- MG[Tr[MG] != Tr_i]

    if (Tr_i == 1) {
      cond_var <- cond_var + mean((Y_i - Y[MG] - ATE) ^ 2)
      cond_var_t <- cond_var_t + mean((Y_i - Y[MG] - ATT) ^ 2)
    }
    else {
      cond_var <- cond_var + mean((Y[MG] - Y_i - ATE) ^ 2)
      cond_var_c <- cond_var_c + mean((Y[MG] - Y_i - ATT) ^ 2)
    }
  }

  cond_var <- cond_var / (2 * n)
  cond_var_t <- cond_var_t / (2 * n1)
  cond_var_c <- cond_var_c / (2 * n0)

  V_sample_sate <- sum(cond_var * (1 + K) ^ 2) / n ^ 2
  V_sample_satt <- sum(cond_var_t * (Tr - (1 - Tr) * K) ^ 2) / n1 ^ 2
  V_sample_satc <- sum(cond_var_c * (Tr * K - (1 - Tr)) ^ 2) / n0 ^ 2

  return(matrix(c(ATE, ATT, ATC, V_sample_sate, V_sample_satt, V_sample_satc),
                ncol = 2,
                dimnames = list(c('All', 'Treated', 'Control'),
                                c('Mean', 'Variance'))))
}

#' @param x An object of class \code{ame}, returned by a call to
#'   \code{\link{FLAME}} or \code{\link{DAME}}.
#' @param digits Number of significant digits for printing the average treatment
#' effect.
#' @param ... Additional arguments to be passed to other methods.
#' @rdname AME
#' @export
print.ame <- function(x, digits = getOption("digits"), ...) {

  df <- x$data

  algo <- x$info$algo

  outcome_type <- x$info$outcome_type

  replacement <- x$info$replacement

  n_iters <- length(x$cov_sets)

  n_matched <- sum(df$matched)
  n_total <- nrow(df)

  cat('An object of class `ame`:\n',
      '  ', algo, ' ran for ', n_iters, ' iterations, matching ',
      n_matched, ' out of ', n_total, ' units ',
      ifelse(replacement, 'with', 'without'), ' replacement.\n', sep = '')


  if (outcome_type == 'continuous') {
    cat('  The average treatment effect of treatment `', x$info$treatment,
        '` on outcome `', x$info$outcome, '` is estimated to be ',
        round(ifelse(x$info$estimate_CATEs && outcome_type == 'continuous',
                      mean(x$data$CATE, na.rm = TRUE),
                      get_average_effects(x)['All', 'Mean']),
               digits = digits),
        '.\n', sep = '')
  }

  if (x$info$missing_data == 'drop') {
    missing_data_message <-
      'Units with missingness in the matching data were not matched.'
  }
  else if (x$info$missing_data == 'impute') {
    missing_data_message <-
      'Missing values in the matching data were imputed by MICE.'
  }
  else if (x$info$missing_data == 'keep') {
    missing_data_message <-
      'Missing values in the matching data were not matched on.'
  }
  else if (x$info$missing_data == 'none') {
    missing_data_message <- NULL
  }

  if (x$info$missing_holdout == 'drop') {
    missing_holdout_message <-
      'Units with missingness in the holdout data were not used to compute PE.'
  }
  else if (x$info$missing_holdout == 'impute') {
    missing_holdout_message <-
      'Missing values in the holdout data were imputed by MICE.'
  }
  else if (x$info$missing_holdout == 'none') {
    missing_holdout_message <- NULL
  }

  if (!is.null(missing_data_message)) {
    cat('  ')
    cat(missing_data_message)
    cat('\n')
  }
  if (!is.null(missing_holdout_message)) {
    cat('  ')
    cat(missing_holdout_message)
    cat('\n')
  }

  return(invisible(x))
}

#' Summarize the output of FLAME or DAME
#'
#' These methods create and print objects of class \code{summary.ame} containing
#' information on the numbers of units matched by the AME algorithm, matched
#' groups formed, and, if applicable, average treatment effects.
#'
#' The average treatment effect (ATE) is estimated as the average CATE estimate
#' across all matched units in the data, while the average treatment effect on
#' the treated (ATT) and average treatment effect on controls (ATC) average only
#' across matched treated or matched control units, respectively. Variances of
#' these estimates are computed as in Abadie, Drukker, Herr, and Imbens (The
#' Stata Journal, 2004) assuming constant treatment effect and homoscedasticity.
#' Note that the implemented estimator is \strong{not} =asymptotically normal and
#' so in particular, asymptotically valid confidence intervals or hypothesis
#' tests cannot be conducted on its basis. In the future, the estimation
#' procedure will be changed to employ the nonparametric regression bias
#' adjustment estimator of Abadie and Imbens 2011 which is asymptotically
#' normal.
#'
#' @return A list of type \code{summary.ame} with the following entries:
#' \describe{
#' \item{MG}{
#'   A list with the number and median size of matched groups formed.
#'   Additionally, two of the highest quality matched groups formed. Quality
#'   is determined first by number of covariates matched on and second by
#'   matched group size.
#'  }
#' \item{n_matches}{
#'   A matrix detailing the number of treated and control units matched.
#'  }
#' \item{TEs}{
#'   If the matching data had a continuous outcome, estimates of the ATE, ATT,
#'   and ATC and the corresponding variance of the estimates.
#'  }
#' }

#' @name summary.ame
NULL
#> NULL

#' @param object An object of class \code{ame}, returned by a call to
#'   \code{\link{FLAME}} or \code{\link{DAME}}
#' @param ... Additional arguments to be passed on to other methods
#' @rdname summary.ame
#' @export
summary.ame <- function(object, ...) {
  df <- object$data
  matched_df <- droplevels(df[df$matched, ]) # anywhere else we need to do this?

  outcome_name <- object$info$outcome
  treated_name <- object$info$treatment

  cov_inds <-
    which(!(colnames(df) %in%
              c(outcome_name, treated_name,
                'matched', 'weight', 'MG', 'CATE')))

  n_matched_on <-
    vapply(seq_len(nrow(df)),
           function(z) n_covs_matched_on(object, z, cov_inds),
           numeric(1))


  if (all(df$matched)) {
    n_match_mat <- matrix(c(sum(df[[treated_name]] == 0),
                            sum(df[[treated_name]] == 1),
                            sum(df[[treated_name]] == 0),
                            sum(df[[treated_name]] == 1),
                            0, 0),
                          byrow = TRUE,
                          ncol = 2,
                          dimnames = list(c('All', 'Matched', 'Unmatched'),
                                          c('Control', 'Treated')))
  }
  else {
    n_matches <- c(table(df[[treated_name]], df$matched))

    n_total <- c(sum(n_matches[c(1, 3)]), sum(n_matches[c(2, 4)]))

    n_match_mat <- matrix(c(n_total, n_matches[3:4], n_matches[1:2]),
                          byrow = TRUE,
                          ncol = 2,
                          dimnames = list(c('All', 'Matched', 'Unmatched'),
                                          c('Control', 'Treated')))
  }

  if (object$info$outcome_type == 'continuous') {
    average_effects <- get_average_effects(object)
  }

  cov_treatment_outcome_inds <-
    !(colnames(df) %in% c('matched', 'weight', 'MG', 'CATE'))

  max_matched_on <- max(n_matched_on)

  matched_on_max <- which(n_matched_on == max_matched_on)
  MG_matched_on_max <- object$MGs[matched_on_max]
  MG_matched_on_max <- MG_matched_on_max[!duplicated(MG_matched_on_max)]
  if (length(MG_matched_on_max) == 1) {
    highest_quality <- MG_matched_on_max[[1]][1]
  }
  else {
    units_matched_on_max <-
      vapply(MG_matched_on_max[!duplicated(MG_matched_on_max)],
             `[`,
             FUN.VALUE = numeric(1),
             1)

    MG_matched_on_max_sizes <- vapply(MG_matched_on_max, length, numeric(1))

    sorted_quality <-
      sort(MG_matched_on_max_sizes, decreasing = TRUE, index.return = TRUE)

    highest_quality <-
      rownames(object$data)[units_matched_on_max[sorted_quality$ix[1:2]]]
  }

  MG_size <- vapply(object$MGs, length, FUN.VALUE = numeric(1))
  MG_number <- sum(MG_size > 0)
  MG_median_size <- median(MG_size[MG_size > 0])

  summary_obj <- list(MG = list(`number` = MG_number,
                                `median_size` = MG_median_size,
                                `highest_quality` = highest_quality),
                      n_matches = n_match_mat)
  if (object$info$outcome_type == 'continuous') {
    summary_obj <- c(summary_obj, list(TEs = average_effects))
  }

  class(summary_obj) <- 'summary.ame'
  return(summary_obj)
}

#' Print a summary of FLAME or DAME
#'
#' @param x An object of class \code{summary.ame}, returned by a call to
#'   \code{\link{summary.ame}}
#' @param digits Number of significant digits for printing the average treatment
#' effect estimates and their variances.
#' @param ... Additional arguments to be passed on to other methods.
#' @rdname summary.ame
#' @export
print.summary.ame <- function(x, digits = 3, ...) {

  max_meanlen <- 7
  max_varlen <- 8
  lablen <- 13

  if ('TEs' %in% names(x)) {

    ATE_meanstr <- format(x$TEs['All', 1], digits = digits, justify = 'right')
    ATE_varstr <- format(x$TEs['All', 2], digits = digits, justify = 'right')
    ATT_meanstr <- format(x$TEs['Treated', 1], digits = digits, justify = 'right')
    ATT_varstr <- format(x$TEs['Treated', 2], digits = digits, justify = 'right')
    ATC_meanstr <- format(x$TEs['Control', 1], digits = digits, justify = 'right')
    ATC_varstr <- format(x$TEs['Control', 2], digits = digits, justify = 'right')

    max_meanlen <- max(max_meanlen, nchar(c(ATE_meanstr, ATT_meanstr, ATC_meanstr)))
    max_varlen <- max(max_varlen, nchar(c(ATE_varstr, ATT_varstr, ATC_varstr)))
  }

  total_line_len <- sum(max_meanlen, max_varlen, 2, lablen)

  cat('Number of Units:\n')

  cat(format('', width = lablen), format('Control', width = max_meanlen, justify = 'right'),
      format('Treated', width = max_varlen, justify = 'right'))
  cat('\n')

  cat(format('  All', width = lablen, justify = 'left'),
      format(x$n_matches['All', 'Control'], width = max_meanlen, justify = 'right'),
      format(x$n_matches['All', 'Treated'], width = max_varlen, justify = 'right'),
      '\n')

  cat(format('  Matched', width = lablen, justify = 'left'),
      format(x$n_matches['Matched', 'Control'], width = max_meanlen, justify = 'right'),
      format(x$n_matches['Matched', 'Treated'], width = max_varlen, justify = 'right'),
      '\n')

  cat(format('  Unmatched', width = lablen, justify = 'left'),
      format(x$n_matches['Unmatched', 'Control'], width = max_meanlen, justify = 'right'),
      format(x$n_matches['Unmatched', 'Treated'], width = max_varlen, justify = 'right'),
      '\n')

  if ('TEs' %in% names(x)) {

    cat('\nAverage Treatment Effects:\n')
    cat(format('', width = lablen),
        format('Mean', width = max_meanlen, justify = 'right'),
        format('Variance', width = max_varlen, justify = 'right'),
        '\n')

    cat(format('  All', width = lablen),
        format(x$TEs['All', 1], digits = digits, width = max_meanlen, justify = 'right'),
        format(x$TEs['All', 2], digits = digits, width = max_varlen, justify = 'right'),
        '\n')

    cat(format('  Treated', width = lablen),
        format(x$TEs['Treated', 1],
               digits = digits, width = max_meanlen, justify = 'right'),
        format(x$TEs['Treated', 2],
               digits = digits, width = max_varlen, justify = 'right'),
        '\n')

    cat(format('  Control', width = lablen),
        format(x$TEs['Control', 1],
               digits = digits, width = max_meanlen, justify = 'right'),
        format(x$TEs['Control', 2],
               digits = digits, width = max_varlen, justify = 'right'),
        '\n')
  }

  cat('\nMatched Groups:\n')
  cat(format('  Number', width = lablen),
      format(x$MG$number, width = total_line_len - lablen, justify = 'right'),
      '\n')

  cat(format('  Median size', width = lablen),
      format(x$MG$median_size, width = total_line_len - lablen, justify = 'right'),
      '\n')
  cat('  Highest quality:', format(ifelse(length(x$MG$highest_quality) == 1,
                                   as.character(x$MG$highest_quality),
                                   paste(x$MG$highest_quality,
                                         collapse =  ' and ')),
      width = total_line_len - 18, justify = 'right'))
  cat('\n')
  return(invisible(x))
}

#' Plot a summary of FLAME or DAME
#'
#' Plot information about numbers of covariates matched on, CATE estimates, and
#' covariate set dropping order after a call to \code{FLAME} or \code{DAME}.
#'
#' \code{plot.ame} displays four plots by default. The first contains
#' information on the number of covariates that matched groups were formed on,
#' and thereby gives some indication of the quality of matched groups across the
#' matched data. The second plots matched group sizes against CATEs, which can
#' be useful for determining whether higher quality matched groups yield
#' different treatment effect estimates than lower quality ones. The third plots
#' a density estimate of the estimated CATE distribution. The fourth displays a
#' heatmap showing which covariates were dropped (shown in black) throughout the
#' matching procedure.
#'
#' @param x An object of class \code{ame}, returned by a call to
#'   \code{\link{FLAME}} or \code{link{DAME}}.
#' @param which_plots A vector describing which plots should be displayed. See
#' details.
#' @param ... Additional arguments to passed on to other methods.
#' @export
plot.ame <- function(x, which_plots = c(1, 2, 3, 4), ...) {

  if (min(which_plots) <= 0 | max(which_plots) >= 5) {
    stop('Please supply an integer 1 through 4 for `which_plots`.')
  }

  # Worry about memory? Should we always do x$data[x$data$matched, ]?
  df <- x$data
  df <- df[df$matched, ]

  n_plots <- length(which_plots)
  n_plotted <- 0
  first_plot <- min(which_plots)

  outcome_name <- x$info$outcome
  treated_name <- x$info$treatment

  Y <- x$data[[outcome_name]]
  Tr <- x$data[[treated_name]]

  cov_inds <-
    !(colnames(df) %in% c(outcome_name, treated_name,
                          'matched', 'weight', 'MG', 'CATE'))

  MGs <- x$MGs

  n_MGs <- length(MGs)
  MG_size <- vapply(x$MGs, length, FUN.VALUE = numeric(1))


  if (!x$info$estimate_CATEs) {
    CATEs <- numeric(n_MGs)
    for (i in seq_len(nrow(x$data))) {
      MG <- MGs[[i]]
      if (is.null(MG)) {
        CATEs[i] <- NA
        next
      }
      if (Tr[i] == 1) {
        CATEs[i] <- Y[i] - mean(Y[MG[Tr[MG] == 0]])
      }
      else {
        CATEs[i] <- mean(Y[MG[Tr[MG] == 1]] - Y[i])
      }
    }
  }
  else {
    CATEs <- x$data$CATE
  }

  MG_size <- MG_size[x$data$matched]
  CATEs <- CATEs[x$data$matched]
  ATE <- mean(CATEs)

  n_covs_matched <-
    vapply(seq_len(nrow(df)),
           function(z) n_covs_matched_on(x, z, cov_inds),
           numeric(1))

  # Number of covariates matched on
  if (1 %in% which_plots) {
    barplot(table(n_covs_matched[n_covs_matched > 0]),
            xlab = 'Number of Covariates Matched On',
            ylab = 'Number of Units')
    n_plotted <- n_plotted + 1
  }

  if (n_plotted == n_plots) {
    return(invisible(x))
  }

  # MG Size vs. CATE
  if (2 %in% which_plots) {
    if (interactive() & first_plot != 2) {
      readline(prompt="Press <enter> to view next plot")
    }
    plot(MG_size, CATEs,
         xlab = 'Size of Matched Group',
         ylab = 'Estimated Conditional Average Treatment Effect')

    if (min(CATEs) < 0 && 0 < max(CATEs)) {
      include_null <- TRUE
    }
    else {
      include_null <- FALSE
    }


    abline(h = ATE, lty = 2)
    if (include_null) {
      abline(h = 0, lty = 3)
    }

    if (include_null) {
      legend('topright',
             legend = c('Estimated ATE', 'Null Effect'), lty = c(2, 3))
    }
    else {
      legend('topright', legend = c('Estimated ATE'), lty = c(2))
    }

    n_plotted <- n_plotted + 1
  }

  if (n_plotted == n_plots) {
    return(invisible(x))
  }

  # CATE Density
  if (3 %in% which_plots) {
    if (interactive() & first_plot != 3) {
      readline(prompt="Press <enter> to view next plot")
    }
    dens <- density(CATEs, na.rm = TRUE)
    if (0 > min(dens$x) && 0 < max(dens$x)) {
      include_null <- TRUE
    }
    else {
      include_null <- FALSE
    }

    plot(dens,
         xlab = c('Estimated Conditional Average Treatment Effect'),
         ylab = '', main = '',
         zero.line = FALSE)

    abline(v = ATE, lty = 2)
    if (include_null) {
      abline(v = 0, lty = 3)
    }

    if (include_null) {
      legend('topright',
             legend = c('Estimated ATE', 'Null Effect'), lty = c(2, 3))
    }
    else {
      legend('topright', legend = c('Estimated ATE'), lty = c(2))
    }

    n_plotted <- n_plotted + 1
  }

  if (n_plotted == n_plots) {
    return(invisible(x))
  }

  # Dropped covariate sets
  if (4 %in% which_plots) {
    if (interactive() & first_plot != 4) {
      readline(prompt="Press <enter> to view next plot")
    }
    cov_sets <- x$cov_sets
    covs_dropped <- matrix(1, nrow = sum(cov_inds), ncol = length(x$cov_sets),
                           dimnames = list(colnames(df)[cov_inds],
                                           seq_along(x$cov_sets)))
    for (i in seq_along(x$cov_sets)) {
      covs_dropped[x$cov_sets[[i]], i] <- 0
    }

    image(z = t(covs_dropped), col = c('black', 'white'),
          xaxt = 'n', yaxt = 'n',
          xlab = 'Iteration', ylab = 'Variables Dropped')

    axis(side = 1, at = seq(0, 1, length.out = ncol(covs_dropped)),
         labels = seq_len(ncol(covs_dropped)))
    axis(side = 2, at = seq(0, 1, length.out = nrow(covs_dropped)),
         labels = rownames(covs_dropped), las = 1)

    return(invisible(x))
  }
}
