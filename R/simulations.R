library(ggplot2)

y_dgp_FLAME <- function(X, Z, eps, p_rel, U = 1, max_nonlinear = 5) {
  n <- nrow(X)

  s <- sample(c(-1, 1), size = p_rel, replace = TRUE)
  alpha <- rnorm(p_rel, 10 *s, 1)
  baseline <- X[, 1:p_rel] %*% alpha

  beta <- rnorm(p_rel, 1.5, sqrt(0.15))
  linear_treated <- X[, 1:p_rel] %*% beta * Z

  interactions <- combn(1:max_nonlinear, 2)
  nonlinear_treated <- rep(0, n)
  for (i in 1:ncol(interactions)) {
    nonlinear_treated %<>%
      magrittr::add(X[, interactions[1, i]] * X[, interactions[2, i]])
  }
  nonlinear_treated <- nonlinear_treated * Z * U

  Y <- baseline + linear_treated + nonlinear_treated + eps

  return(Y)
}

# FLAME Simulation 6.1.1
sim1 <- FALSE
set.seed(910916)
n_control <- 10000
n_treated <- 10000
n <- n_control + n_treated
p <- 10

treatment <- rep(0, n)
treatment[sample(1:n, size = n / 2)] <- 1
X <-
  rbinom(n * p, size = 1, prob = 0.5) %>%
  matrix(nrow = n, ncol = p)

eps <- rnorm(n, 0, sqrt(0.1))

Y0 <- y_dgp_FLAME(X, rep(0, n), eps, p_rel = p, U = 10, max_nonlinear = 5)
Y1 <- y_dgp_FLAME(X, rep(1, n), eps, p_rel = p, U = 10, max_nonlinear = 5)
HTE <- Y1 - Y0
Y <-  Y1 * treatment + Y0 * (1 - treatment)

HTEs <- HTE[1:(n / 2), ]

df <-
  cbind(X, treatment, Y) %>%
  as.data.frame()

holdout <- df[(n / 2 + 1):n, ]
df <- df[1:(n / 2), ]

colnames(df) <- c(colnames(df)[1:p], 'treated', 'outcome')
colnames(holdout) <- c(colnames(df)[1:p], 'treated', 'outcome')
df[, 1:p] %<>% lapply(as.factor)
holdout[, 1:p] %<>% lapply(as.factor)

if (sim1) {
  flout <- FLAME(df, holdout, verbose = 3)
  CATEs <- group_CATE(1:(n / 2), flout)
  unmatched <- which(!flout$data$matched)
  error <- abs(HTEs[-unmatched] - unlist(CATEs))

  ggdata <- data.frame(true = HTEs[-unmatched],
                       predicted = unlist(CATEs))

  ggplot(ggdata, aes(x = true, y = predicted)) +
    geom_point() +
    geom_abline(slope = 1, color = 'red')
}

# FLAME Simulation 6.1.2
sim2 <- FALSE
set.seed(910916)
n_control <- 10000
n_treated <- 10000
n <- n_control + n_treated
p <- 30
p_rel <- 10

treatment <- rep(0, n)
treatment[sample(1:n, size = n / 2)] <- 1
X_rel <-
  rbinom(n * p_rel, size = 1, prob = 0.5) %>%
  matrix(nrow = n, ncol = p_rel)

X <- X_rel

for (i in 1:(p - p_rel)) {
  X <- cbind(X, matrix(rbinom(n, size = 1, prob = 0.1 + 0.8 * treatment), nrow = n))
}

eps <- rnorm(n, 0, sqrt(0.1))

Y0 <- y_dgp_FLAME(X, rep(0, n), eps, p_rel = p_rel, U = 1, max_nonlinear = 5)
Y1 <- y_dgp_FLAME(X, rep(1, n), eps, p_rel = p_rel, U = 1, max_nonlinear = 5)
HTE <- Y1 - Y0
Y <-  Y1 * treatment + Y0 * (1 - treatment)

HTEs <- HTE[1:(n / 2), ]

df <-
  cbind(X, treatment, Y) %>%
  as.data.frame()

holdout <- df[(n / 2 + 1):n, ]
df <- df[1:(n / 2), ]

colnames(df) <- c(colnames(df)[1:p], 'treated', 'outcome')
colnames(holdout) <- c(colnames(df)[1:p], 'treated', 'outcome')
df[, 1:p] %<>% lapply(as.factor)
holdout[, 1:p] %<>% lapply(as.factor)
if (sim2) {
  flout <- FLAME(df, holdout, verbose = 3)
  CATEs <- group_CATE(1:(n / 2), flout)
  unmatched <- which(!flout$data$matched)
  error <- abs(HTEs[-unmatched] - unlist(CATEs))

  ggdata <- data.frame(true = HTEs[-unmatched],
                       predicted = unlist(CATEs))

  ggplot(ggdata, aes(x = true, y = predicted)) +
    geom_point() +
    geom_abline(slope = 1, color = 'red')
}
