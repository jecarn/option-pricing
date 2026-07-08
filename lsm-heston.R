# American call for max of J uncorrelated assets

# Libraries
library(matrixStats)
library(Rfast)

# Parameters
vol <- 0.2 # volatility (sigma)
S_0 <- 110 # initial stock price
r <- 0.05 # risk-free rate
d <- 0.1 # proportional dividend rate
K <- 100 # strike price
T <- 3 # time to maturity (years)

# N <- 5 # number of time steps
# R <- 20 # number of replications / asset
# J <- 5 # number of assets

N <- 252 # number of time steps
R <- 10000 # number of replications / asset
J <- 5 # number of assets

cor_matrix <- matrix(
  c(1, 0, 0, 0, 0,
    0, 1, 0, 0, 0,
    0, 0, 1, 0, 0,
    0, 0, 0, 1, 0,
    0, 0, 0, 0, 1),
  nrow = J, ncol = J, byrow = TRUE
)

# svd(cor_matrix)

start_time <- Sys.time()

# Model values
dt <- T / N # time increment

# Path simulations
Z <- matrix(rnorm(N * R * J), nrow = R * J, ncol = N) # matrix of standard normal RVs
dW <- Z * sqrt(dt) # Brownian increments
W <- rowCumsums(dW) # Brownian motion at each time
deterministic_matrix <- matrix((r - d - vol^2 / 2) * seq(dt, T, dt), nrow = R * J, ncol = N, byrow = TRUE) # deterministic component of GMB
X <- S_0 * exp(deterministic_matrix) * exp(vol * W) # simulated GBMs

path_maxes <- matrix(0, nrow = R, ncol = N)

for (i in 0 : (R - 1)){
  path_maxes[i + 1, ] <- colMaxs(X[(1 + i * J) : (J * (i + 1)), ], value = TRUE)
}

exercise_matrix <- ifelse(path_maxes - K <= 0, 0, path_maxes - K) # matrix of immediate exercise values

cf_matrix <- matrix(c(rep(0, R * (N - 1)), exercise_matrix[, N]), nrow = R, ncol = N) # cashflow matrix

matrix_disc <- function(mx, r, n, dt){ # function to discount a matrix (indexed by column)
  disc_vec <- exp(-r * dt * (1 : n))
  return(t(t(mx) * disc_vec))
}

hermite_basis <- function(x, n){ # takes a vector x and returns the first n polynomials as a matrix
  rows <- length(x)
  He <- matrix(1, nrow = rows, ncol = n + 1)
  He[, 2] <- x
  
  for (i in 1 : (n - 1)){
    He[, i + 2] <- x * He[, i + 1] - i * He[, i]
  }
  return(He[, -1])
}

for (j in rev(seq(N - 1))){ # iterating over time steps
  path_regr_filter <- ifelse(path_maxes[, j] > K, TRUE, FALSE) # filtering for positive exercise values
  regr_filter <- rep(path_regr_filter, each = J)
  
  # x <- X[, j, drop = FALSE][regr_filter] / K # basis variable
  x <- matrix(X[regr_filter, j], nrow = sum(path_regr_filter), ncol = J, byrow = TRUE)
  x_sorted <- rowSort(x, descending = TRUE)
  
  n <- 5 # number of Hermite basis polynomials
  
  x_matrix <- matrix(0, ncol = 3 * J + n - 2, nrow = sum(path_regr_filter))
  
  x_matrix[, 1 : n] <- hermite_basis(x_sorted[, 1], n) # adding Hermite basis polynomials for max prices
  x_matrix[, (n + 1) : (n + J - 1)] <- x_sorted[, -1] # adding 2nd, 3rd, ... highest prices
  x_matrix[, (n + J) : (n + 2 * (J - 1))] <- x_sorted[, -1]^2 # adding squares of 2nd, 3rd, ...
  x_matrix[, (n + 2 * J - 1) : (n + 3 * J - 3)] <- x_sorted[, 1 : (J - 1)] * x_sorted[, 2 : J] # adding 1st * 2nd, 2nd * 3rd, ...
  x_matrix[, 3 * J + n - 2] <- rowprods(x_sorted)
  
  y <- rowSums(matrix_disc(cf_matrix[path_regr_filter, (j + 1) : N], r, N - j, dt)) # includes discounted future CFs
  
  ls_model <- lm( # as in L&S
    y ~ x_matrix
  )
  
  exp_cont <- fitted(ls_model) # expected continuation value
  imm_exer <- exercise_matrix[path_regr_filter, j] # immediate exercise value
  
  exer_filter <- ifelse(exp_cont < imm_exer, TRUE, FALSE) # filtering for cases of greater immediate exercise value
  
  comb_filter <- which(path_regr_filter)[exer_filter] # combining regression and exercise filters
  
  cf_matrix[comb_filter, (j + 1) : N] = 0 # setting future CFs to zero for early exercise paths
  cf_matrix[comb_filter, j] = imm_exer[exer_filter] # setting current CF to immediate exercise value
}

disc_cf_matrix <- matrix_disc(cf_matrix, r, N, dt) # discounting all CFs

max(sum(disc_cf_matrix) / R, S_0 - K) # call value

end_time <- Sys.time()

perf_time <- end_time - start_time
perf_time

lsm <- X[, N]

plot(density(X[, N]))

