# Thanks to Dirk Eddelbuettel:
#  https://stackoverflow.com/questions/36166288/skip-tests-on-cran-but-run-locally
if (length(strsplit(packageDescription("FLAME")$Version, "\\.")[[1]]) > 3) {
  skipping_on_cran <- FALSE
} else {
  skipping_on_cran <- TRUE
}

if (!skipping_on_cran) {
  n <- 500
  p <- 4
  data <- gen_data(n, p)
  holdout <- gen_data(n, p)

  test_that("ridge works with continuous outcome", {
    flout <- FLAME(data, holdout, PE_method = 'ridge')
    expect_true(TRUE)
  })

  test_that("XGBoost works with continuous outcome", {
    flout <- FLAME(data, holdout, PE_method = 'xgb')
    expect_true(TRUE)
  })

  data$outcome <-
    exp(scale(data$outcome, scale = F)) /
    (1 + exp(scale(data$outcome, scale = F)))

  data$outcome <- round(data$outcome)

  holdout$outcome <-
    exp(scale(holdout$outcome, scale = F)) /
    (1 + exp(scale(holdout$outcome, scale = F)))

  holdout$outcome <- round(holdout$outcome)

  test_that("ridge works with binary continuous outcome", {
    flout <- FLAME(data, holdout, PE_method = 'ridge')
    expect_true(TRUE)
  })

  test_that("XGBoost works with binary continuous outcome", {
    flout <- FLAME(data, holdout, PE_method = 'xgb')
    expect_true(TRUE)
  })

  data$outcome <- factor(data$outcome)
  holdout$outcome <- factor(holdout$outcome)

  test_that("ridge works with binary factor outcome", {
    flout <- FLAME(data, holdout, PE_method = 'ridge')
    expect_null(flout$CATE)
  })

  test_that("XGBoost works with binary factor outcome", {
    flout <- FLAME(data, holdout, PE_method = 'xgb')
    expect_null(flout$CATE)
  })

  data$outcome <- factor(sample(c('Green', 'White', 'Red'), n, TRUE))
  holdout$outcome <- factor(sample(c('Green', 'White', 'Red'), n, TRUE))

  test_that("ridge works with multiclass outcomes", {
    flout <- FLAME(data, holdout, PE_method = 'ridge')
    expect_null(flout$CATE)
  })

  test_that("XGBoost works with multiclass outcomes", {
    flout <- FLAME(data, holdout, PE_method = 'xgb')
    expect_null(flout$CATE)
  })

  ######
  # test_that("independent of outcome levels", {
  #   data <- gen_data(n, p)
  #   holdout <- gen_data(n, p)
  #
  #
  #   data$outcome <-
  #     exp(scale(data$outcome, scale = F)) /
  #     (1 + exp(scale(data$outcome, scale = F)))
  #
  #   data$outcome <- round(data$outcome)
  #
  #   holdout$outcome <-
  #     exp(scale(holdout$outcome, scale = F)) /
  #     (1 + exp(scale(holdout$outcome, scale = F)))
  #
  #   holdout$outcome <- round(holdout$outcome)
  #
  #   flout <- FLAME(data, holdout)
  #
  #   data$outcome <- factor(data$outcome)
  #   holdout$outcome <- factor(holdout$outcome)
  #
  #   flout_factor <- FLAME(data, holdout)
  #   if (identical(flout$dropped, flout_factor$dropped)) {
  #     expect_identical(flout$MGs, flout_factor$MGs)
  #   }
  #
  #   levels(data$outcome) <- c('A', 'B')
  #   levels(holdout$outcome) <- c('A', 'B')
  #
  #   flout_factor2 <- FLAME(data, holdout)
  #   if (identical(flout$dropped, flout_factor2$dropped)) {
  #     expect_identical(flout$MGs, flout_factor2$MGs)
  #   }
  # })
}
