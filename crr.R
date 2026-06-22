# Parameters
vol <- sqrt(0.2)
S_0 <- 100
K <- 100
r <- 0.05
T <- 1

N <- 3

start_time <- Sys.time()

# Calculations
dt <- T / N
u <- exp(vol * dt)
d <- 1 / u

p <- (exp(r * dt) - d) / (u - d)
q <- 1 - p

# Functions
prices_at_t <- function(S_0, u, d, n){
  prices_vec <- c()
  for (j in seq(0, n)){
    prices_vec <- c(prices_vec, u^j * d^(n - j) * S_0)
  }
  return(prices_vec)
}

payoff <- function(S, K){
  return(ifelse(K - S <= 0, 0, K - S))
}

exp_disc_payoff <- function(p, q, r, dt, future_payoff){
  m <- length(future_payoff)
  payoff_vec <- c()
  for (j in seq(m - 1)){
    payoff_vec <- c(payoff_vec, exp(-r * dt) * (q * future_payoff[j] + p * future_payoff[j + 1]))
  }
  return(payoff_vec)
}
prices_at_t(S_0, u, d, 1)
future_payoff <- payoff(prices_at_t(S_0, u, d, N), K)

# exp_disc_payoff(p, q, r, dt, future_payoff)

for (j in rev(seq(N - 1))){
  exp_disc_poff <- exp_disc_payoff(p, q, r, dt, future_payoff)
  current_prices <- prices_at_t(S_0, u, d, j)
  future_payoff <- pmax(exp_disc_poff, payoff(current_prices, K))
}

exp_disc_payoff(p, q, r, dt, future_payoff)

end_time <- Sys.time()

perf_time <- end_time - start_time

