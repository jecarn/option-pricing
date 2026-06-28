# Libraries
library(matrixStats)

# Parameters
vol <- sqrt(0.1) # volatility (sigma)
S_0 <- 100 # initial stock price
r <- 0.06 # risk-free rate
K <- 110 # strike price
T <- 1 # time to maturity (years)

N <- 50 # number of time steps
R <- 50000 # number of replications

start_time <- Sys.time()

# Model values
dt <- T / N # time increment

# Path simulations
Z <- matrix(rnorm(N * R), nrow = R, ncol = N) # matrix of standard normal RVs
dW <- Z * sqrt(dt) # Brownian increments
W <- rowCumsums(dW) # Brownian motion at each time
deterministic_matrix <- matrix((r - vol^2 / 2) * seq(dt, T, dt), nrow = R, ncol = N, byrow = TRUE) # deterministic component of GMB
X <- S_0 * exp(deterministic_matrix) * exp(vol * W) # simulated GBMs

exercise_matrix <- ifelse(K - X <= 0, 0, K - X) # matrix of immediate exercise values

cf_matrix <- matrix(c(rep(0, R * (N - 1)), exercise_matrix[, N]), nrow = R, ncol = N) # cashflow matrix

matrix_disc <- function(mx, r, n, dt){ # function to discount a matrix (indexed by column)
  disc_vec <- exp(-r * dt * (1 : n))
  return(t(t(mx) * disc_vec))
}

for (j in rev(seq(N - 1))){ # iterating over time steps
  regr_filter <- ifelse(X[, j] < K, TRUE, FALSE) # filtering for positive exercise values
  
  x <- X[, j][regr_filter] / K # basis variable
  y <- rowSums(matrix_disc(cf_matrix[regr_filter, (j + 1) : N], r, N - j, dt)) # includes discounted future CFs
  
  ls_model <- lm( # regressing on 3 Laguerre polynomials
    y ~ I(exp(-x / 2)) + I(exp(-x / 2) * (1 - x)) + I(exp(-x / 2) * (1 - 2 * x + x^2 / 2))
  )
  
  exp_cont <- fitted(ls_model) # expected continuation value
  imm_exer <- exercise_matrix[regr_filter, j] # immediate exercise value
  
  exer_filter <- ifelse(exp_cont < imm_exer, TRUE, FALSE) # filtering for cases of greater immediate exercise value
  
  comb_filter <- which(regr_filter)[exer_filter] # combining regression and exercise filters
  
  cf_matrix[comb_filter, (j + 1) : N] = 0 # setting future CFs to zero for early exercise paths
  cf_matrix[comb_filter, j] = imm_exer[exer_filter] # setting current CF to immediate exercise value
}

disc_cf_matrix <- matrix_disc(cf_matrix, r, N, dt) # discounting all CFs

max(sum(disc_cf_matrix) / R, K - S_0) # put value

end_time <- Sys.time()

perf_time <- end_time - start_time
perf_time
