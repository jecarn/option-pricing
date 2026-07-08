# Libraries
library(matrixStats)

# Option parameters
vol <- 0.2 # volatility (sigma)
S_0 <- 100 # initial stock price
r <- 0.06 # risk-free rate
K <- 110 # strike price
T <- 1 # time to maturity (years)

# Merton model parameters
lambda <- 1.5 # jump intensity
mu_j <- -0.1 # mean of log jump size
sigma_j <- 0.15 # volatility of log jump size

N <- 252 # number of time steps
R <- 50000 # number of replications

start_time <- Sys.time()

# Model values
dt <- T / N # time increment
k_bar <- exp(mu_j + sigma_j^2 / 2) - 1 # expected jump compensation

# Path simulations
Z <- matrix(rnorm(N * R), nrow = R, ncol = N) # matrix of standard normal RVs
dW <- Z * sqrt(dt) # Brownian increments

# Jump simulations
dq <- matrix(ifelse(rpois(N * R, lambda * dt) > 0, 1, 0), nrow = R, ncol = N, byrow = TRUE) # Poisson occurances dq_t
k <- matrix(exp(rnorm(N * R, mu_j, sigma_j)) - 1, nrow = R, ncol = N, byrow = TRUE) # random normal jump sizes

S <- matrix(S_0, ncol = (N + 1), nrow = R)

for (n in 1 : N){
  S_n <- S[, n]
  S[, n + 1] <- S_n + (r - lambda * k_bar) * S_n * dt + vol * S_n * dW[, n] + S_n * k[, n] * dq[, n]
  # S[, n + 1] <- S_n * exp((r - V_n / 2) * dt + sqrt(V_n) * dW_S[, n]) # alternative discretization
}

X <- S[, -1] # asset prices

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
  
  ls_model <- lm( # regressing on 3 weighted Laguerre polynomials
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

max(sum(disc_cf_matrix) / R, K - S_0) # American put option value

end_time <- Sys.time()

perf_time <- end_time - start_time
perf_time

# DIAGNOSTICS

exp(-r * T) * mean(payoff(X[, N], K)) # European put option price
exp(-r * T) * mean(pmax(X[, N] - K, 0)) # European call option price

# lsm_heston <- X[, N]

plot(density(X[, N]))

mean(X[, N])

