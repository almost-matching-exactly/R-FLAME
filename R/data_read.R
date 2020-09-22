read_data <- function(data, holdout,
                      treated_column_name, outcome_column_name) {
  if (is.character(data)) {
    tryCatch(
      error = function(cnd) {
        stop('Cannot read `data` .csv file from working directory')
      },
      data <- read.csv(data, header = TRUE)
    )
    cov_inds <-
      which(!(colnames(data) %in% c(treated_column_name, outcome_column_name)))
    for (cov in cov_inds) {
      data[[cov]] <- as.factor(data[[cov]])
    }
  }

  if (is.character(holdout)) {
    tryCatch(
      error = function(cnd) {
        stop('Cannot read `holdout` .csv file from working directory')
      },
      holdout <- read.csv(holdout, header = TRUE)
    )
    cov_inds <- which(!(colnames(holdout) %in%
                        c(treated_column_name, outcome_column_name)))

    for (cov in cov_inds) {
      holdout[[cov]] <- as.factor(holdout[[cov]])
    }
  }

  if (is.numeric(holdout) & length(holdout) == 1) {
      holdout_inds <- sample(1:nrow(data), size = round(holdout * nrow(data)))
      holdout <- data[holdout_inds, ]
      data <- data[-holdout_inds, ]
  }

  return(list(data = data,
              holdout = holdout))
}
