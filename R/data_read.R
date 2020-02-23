read_data <- function(data, holdout) {
  if (is.character(data)) {
    tryCatch(
      error = function(cnd) {
        stop('Cannot read data .csv file from current directory')
      },
      data <- read.csv(data, header = TRUE)
    )
  }

  if (is.character(holdout)) {
    tryCatch(
      error = function(cnd) {
        stop('Cannot read holdout .csv file from current directory')
      },
      holdout <- read.csv(holdout, header = TRUE)
    )
  }

  if (is.logical(holdout)) {
      holdout_inds <- sample(1:nrow(data), size = round(0.1 * nrow(data)))
      holdout <- data[holdout_inds, ]
      data <- data[-holdout_inds, ]
  }

  return(list(data = data,
              holdout = holdout))
}
