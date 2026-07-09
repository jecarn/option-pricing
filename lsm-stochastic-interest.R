# Libraries
library(matrixStats)

# Option parameters
vol <- sqrt(0.1)
S_0 <- 100 # initial stock price
K <- 110 # strike price
T <- 1 # time to maturity (years)

# CIR model parameters
a <- 0.5 # mean-reversion rate
b <- 0.06 # mean risk-free rate
sigma <- 0.15 # variance of risk-free rate
r_0 <- 0.06 # initial risk-free rate

# N <- 252 # number of time steps
# R <- 100000 # number of replications

N <- 252 # number of time steps
R <- 50000 # number of replications

start_time <- Sys.time()

# Model values
dt <- T / N # time increment

# Path simulations
Z_1 <- matrix(rnorm(N * R), nrow = R, ncol = N) # matrix of standard normal RVs
Z_2 <- matrix(rnorm(N * R), nrow = R, ncol = N)
dW_S <- Z_1 * sqrt(dt) # Brownian increments
dW_r <- Z_2 * sqrt(dt)

r <- matrix(r_0, ncol = (N + 1), nrow = R)
S <- matrix(S_0, ncol = (N + 1), nrow = R)

for (n in 1 : N){
  r_n <- r[, n]
  r[, n + 1] <- pmax(r_n + a * (b - r_n) * dt + sigma * sqrt(r_n) * dW_r[, n], 0)
  
  S_n <- S[, n]
  # S[, n + 1] <- S_n + r_n * S_n * dt + S_n * vol * dW_S[, n]
  S[, n + 1] <- S_n * exp((r_n - vol^2 / 2) * dt + vol * dW_S[, n]) # alternative discretization
}

X <- S[, -1] # asset prices
r <- r[, -1] # risk-free rate paths

exercise_matrix <- ifelse(K - X <= 0, 0, K - X) # matrix of immediate exercise values

cf_matrix <- matrix(c(rep(0, R * (N - 1)), exercise_matrix[, N]), nrow = R, ncol = N) # cashflow matrix

# matrix_disc <- function(mx, r, n, dt){ # function to discount a matrix (indexed by column)
#   disc_vec <- exp(-r * dt * (1 : n))
#   return(t(t(mx) * disc_vec))
# }

for (j in rev(seq(N - 1))){ # iterating over time steps
  regr_filter <- ifelse(X[, j] < K, TRUE, FALSE) # filtering for positive exercise values
  
  # if (!(sum(regr_filter) == 0)){
  #   
  # }
  
  x <- X[, j][regr_filter] / K # basis variable
  
  avg_rf_matrix <- t(t(rowCumsums(r[regr_filter, (j + 1) : N, drop = FALSE])) / c(1 : (N - j)))
  
  disc_matrix <- t(t(avg_rf_matrix) * dt * (1 : (N - j)))
  
  y <- rowSums(cf_matrix[regr_filter, (j + 1) : N, drop = FALSE] * exp(-disc_matrix)) # includes discounted future CFs
  
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

# disc_cf_matrix <- matrix_disc(cf_matrix, r, N, dt) # discounting all CFs (DOESNT WORK BECAUSE r IS NOT A CONSTANT)

avg_rf_matrix <- t(t(rowCumsums(r)) / c(1 : N))

disc_matrix <- t(t(avg_rf_matrix) * dt * (1 : N))

disc_cf_matrix <- cf_matrix * exp(-disc_matrix)

max(sum(disc_cf_matrix) / R, K - S_0) # American put option value

end_time <- Sys.time()

perf_time <- end_time - start_time
perf_time

# DIAGNOSTICS

# exp(-r * T) * mean(payoff(X[, N], K)) # European put option price
# exp(-r * T) * mean(pmax(X[, N] - K, 0)) # European call option price

lsm_cir <- X[, N]

plot(density(X[, N]))

plot(density(lsm))
lines(density(lsm_heston))
lines(density(lsm_cir))
