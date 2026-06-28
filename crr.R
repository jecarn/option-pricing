# Parameters
vol <- sqrt(0.1) # volatility (sigma)
S_0 <- 100 # initial stock price
r <- 0.06 # risk-free rate
K <- 110 # strike price
T <- 1 # time to maturity (years)

N <- 200 # number of time steps

start_time <- Sys.time()

# Model values
dt <- T / N # time increment
u <- exp(vol * sqrt(dt)) # upturn factor
d <- 1 / u # downturn factor

p <- (exp(r * dt) - d) / (u - d) # upturn probability
q <- 1 - p # downturn probability

# Functions
prices_at_t <- function(S_0, u, d, n){ # returns the vector of prices at time step n
  prices_vec <- c()
  for (j in seq(0, n)){
    prices_vec <- c(prices_vec, u^j * d^(n - j) * S_0)
  }
  return(prices_vec)
}

payoff <- function(S, K){ # returns the non-negative payoff
  return(ifelse(K - S <= 0, 0, K - S))
}

exp_disc_payoff <- function(p, q, r, dt, future_payoff){ # discounts expected future payoffs
  m <- length(future_payoff)
  payoff_vec <- c()
  for (j in seq(m - 1)){
    payoff_vec <- c(payoff_vec, exp(-r * dt) * (q * future_payoff[j] + p * future_payoff[j + 1]))
  }
  return(payoff_vec)
}

# Calculations
future_payoff <- payoff(prices_at_t(S_0, u, d, N), K)

for (j in rev(seq(N - 1))){ # iterating over time steps
  exp_disc_poff <- exp_disc_payoff(p, q, r, dt, future_payoff)
  current_prices <- prices_at_t(S_0, u, d, j)
  future_payoff <- pmax(exp_disc_poff, payoff(current_prices, K))
}

exp_disc_payoff(p, q, r, dt, future_payoff) # option value

end_time <- Sys.time()

perf_time <- end_time - start_time
perf_time

