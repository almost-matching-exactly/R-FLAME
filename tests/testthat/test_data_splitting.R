n <- 250
p <- 5
data <- gen_data(5 * n, p)
weights <- runif(p)

test_that("pre-splitting data doesn't affect post matching", {

  holdout_inds <- sample(seq_len(nrow(data)), round(0.1 * nrow(data)), TRUE)
  holdout <- data[holdout_inds, ]
  data <- data[-holdout_inds, ]

  flout <- FLAME(data, holdout, weights = weights)
  for (i in seq_len(nrow(flout$data))) {
    i <- as.numeric(i)
    inds <- flout$MGs[[i]]
    tmp <-
      MG(as.numeric(rownames(flout$data))[i], flout, index_only = TRUE)[[1]]
    if (is.null(inds)) {
      expect_null(tmp)
    }
    else {
      expect_equivalent(sort(inds), sort(tmp))
    }
  }
})

test_that("auto-splitting data doesn't affect post matching", {
  flout <- FLAME(data, weights = weights)
  for (i in seq_len(nrow(flout$data))) {
    i <- as.numeric(i)
    inds <- flout$MGs[[i]]
    tmp <-
      MG(as.numeric(rownames(flout$data))[i], flout, index_only = TRUE)[[1]]
    if (is.null(inds)) {
      expect_null(tmp)
    }
    else {
      expect_equivalent(sort(inds), sort(tmp))
    }
  }
})

test_that("scrambled rownames and auto-splitting data
          doesn't affect post matching", {

  rownames(data) <- sample(nrow(data))

  flout <- FLAME(data, weights = weights)
  for (i in seq_len(nrow(flout$data))) {
    i <- as.numeric(i)
    inds <- flout$MGs[[i]]
    tmp <-
      MG(as.numeric(rownames(flout$data))[i], flout, index_only = TRUE)[[1]]
    if (is.null(inds)) {
      expect_null(tmp)
    }
    else {
      expect_equivalent(sort(inds), sort(tmp))
    }
  }
})

test_that("can't require data splitting without outcome", {
  data$outcome <- NULL
  expect_error(flout <- FLAME(data, weights = weights))
})
