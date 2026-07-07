# Parameters
vol <- 0.2 # volatility (sigma)
S_0 <- 100 # initial stock price
r <- 0.06 # risk-free rate
K <- 100 # strike price
T <- 2 # time to maturity (years)

d_1 <- (log(S_0 / K) + (r + vol^2 / 2) * T) / (vol * sqrt(T))
d_2 <- d_1 - vol * sqrt(T)

C <- S_0 * pnorm(d_1) - K * exp(-r * T) * pnorm(d_2) # European call option price
P <- C + K * exp(-r * T) - S_0 # European put option price

C
P