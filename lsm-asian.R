# American Asian call option

# Libraries
library(matrixStats)
library(Rfast)

# Parameters
vol <- 0.2 # volatility (sigma)
S_0 <- 100 # initial stock price
A_0 <- 100 # initial average price
r <- 0.06 # risk-free rate
K <- 100 # strike price
T <- 2 # time to maturity (years)
L <- 0.25 # lookback period (years)
E <- 0.25 # exercise start (years)

N <- 250 # number of time steps
R <- 100000 # number of replications

start_time <- Sys.time()

# Model values
dt <- T / N # time increment
l <- max(floor(L / dt), 1) # lookback time steps
e <- ceiling(E / dt) # exercise start time step

# Path simulations
Z <- matrix(rnorm(N * R), nrow = R, ncol = N) # matrix of standard normal RVs
dW <- Z * sqrt(dt) # Brownian increments
W <- rowCumsums(dW) # Brownian motion at each time
deterministic_matrix <- matrix((r - vol^2 / 2) * seq(dt, T, dt), nrow = R, ncol = N, byrow = TRUE) # deterministic component of GMB
S <- S_0 * exp(deterministic_matrix) * exp(vol * W) # simulated GBMs

X <- matrix(nrow = R, ncol = (l + 1 + N))

X[, 1 : l] <- A_0

X[, (l + 1)] <- S_0

X[, (l + 2) : (l + 1 + N)] <- S

A <- X

for (n in (2 : (l + 1 + N))){
  A[, n] <- rowmeans(X[, 1 : n])
}

A_exc <- A[, (l + 2) : (l + 1 + N)]

exercise_matrix <- ifelse(A_exc - K <= 0, 0, A_exc - K) # matrix of immediate exercise values

cf_matrix <- matrix(c(rep(0, R * (N - 1)), exercise_matrix[, N]), nrow = R, ncol = N) # cashflow matrix

matrix_disc <- function(mx, r, n, dt){ # function to discount a matrix (indexed by column)
  disc_vec <- exp(-r * dt * (1 : n))
  return(t(t(mx) * disc_vec))
}

for (j in seq(N - 1, e)){ # iterating over time steps
  regr_filter <- ifelse(A_exc[, j] > K, TRUE, FALSE) # filtering for positive exercise values
  
  x_spot <- S[, j][regr_filter] / K # stock price basis variable
  x_avg <- A_exc[, j][regr_filter] / K # stock average basis variable
  y <- rowSums(matrix_disc(cf_matrix[regr_filter, (j + 1) : N, drop = FALSE], r, N - j, dt)) # includes discounted future CFs
  
  if (sum(regr_filter != 0) & length(y) != 0){
    ls_model <- lm( # regressing on 3 Laguerre polynomials
      y ~ I(1 - x_spot) + I(1 - 2 * x_spot + x_spot^2 / 2) +
        I(1 - x_avg) + I(1 - 2 * x_avg + x_avg^2 / 2) +
        I((1 - x_spot) * (1 - x_avg)) +
        I((1 - x_spot) * (1 - 1 * x_avg + x_avg^2 / 2)) + I((1 - x_avg) * (1 - 1 * x_spot + x_spot^2 / 2))
    )
    
    exp_cont <- fitted(ls_model) # expected continuation value
    imm_exer <- exercise_matrix[regr_filter, j] # immediate exercise value
    
    exer_filter <- ifelse(exp_cont < imm_exer, TRUE, FALSE) # filtering for cases of greater immediate exercise value
    
    comb_filter <- which(regr_filter)[exer_filter] # combining regression and exercise filters
    
    cf_matrix[comb_filter, (j + 1) : N] = 0 # setting future CFs to zero for early exercise paths
    cf_matrix[comb_filter, j] = imm_exer[exer_filter] # setting current CF to immediate exercise value
  }
}

disc_cf_matrix <- matrix_disc(cf_matrix, r, N, dt) # discounting all CFs

max(sum(disc_cf_matrix) / R, S_0 - K) # call value

end_time <- Sys.time()

perf_time <- end_time - start_time
perf_time

plot(density(A_exc[, N]))

exp(-r * T) * mean(pmax(K - A_exc[, N], 0)) # European put option price
exp(-r * T) * mean(pmax(A_exc[, N] - K, 0)) # European call option price

# exp(-r * T) * mean(S[, N])
#exp(-r * T) * mean(A[, N])

# exp(-r * T) * mean(pmax(A_T - K, 0))

