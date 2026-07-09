# https://arxiv.org/pdf/1502.02963

# Option parameters
S_0 <- 100 # initial stock price
r <- 0.06 # risk-free rate
K <- 110 # strike price
T <- 1 # time to maturity (years)

# Heston model parameters
rho <- -0.7 # correlation between asset and variance Brownian motions
kappa <- 2 # mean-reversion rate
theta <- 0.1 # mean variance
v_0 <- 0.11 # initial variance
xi <- 0.4 # volatility of volatility

# Heston characteristic function
heston_char_func <- function(w, S_0, K, r, T, kappa, theta, xi, rho, v_0){
  alpha <- -w^2 / 2 - 1i * w / 2
  beta <- kappa - rho * xi * 1i * w
  gamma <- xi^2 / 2
  h <- sqrt(beta^2 - 4 * alpha * gamma)
  
  r_p <- (beta + h) / xi^2
  r_m <- (beta - h) / xi^2
  g <- r_m / r_p
  
  C <- kappa * (r_m * T - 2 * log((1 - g * exp(-h * T)) / (1 - g)) / xi^2)
  D <- r_m * (1 - exp(-h * T)) / (1 - g * exp(-h * T))
  
  return(exp(C * theta + D * v_0 + 1i * w * log(S_0 * exp(r * T))))
}

# Heston call and put prices
heston_call <- function(S_0, K, r, T, kappa, theta, xi, rho, v_0){
  Psi <- function(w){heston_char_func(w, S_0 = S_0, K = K, r = r, T = T, kappa = kappa, theta = theta, xi = xi, rho = rho, v_0 = v_0)}
  integrand_1 <- function(w){Re((exp(-1i * w * log(K)) * Psi(w - 1i)) / (1i * w * Psi(-1i)))}
  integrand_2 <- function(w){Re((exp(-1i * w * log(K)) * Psi(w)) / (1i * w))}
  
  integral_1 <- integrate(
    integrand_1,
    lower = 0, upper = Inf
  )
  integral_2 <- integrate(
    integrand_2,
    lower = 0, upper = Inf
  )
  
  Pi_1 <- 1/2 + integral_1$value / pi
  Pi_2 <- 1/2 + integral_2$value / pi
  
  return(S_0 * Pi_1 - exp(-r * T) * K * Pi_2)
}

heston_put <- function(S_0, K, r, T, kappa, theta, xi, rho, v_0){
  C <- heston_call(S_0, K, r, T, kappa, theta, xi, rho, v_0)
  return(C + exp(-r * T) * K - S_0) # calculated using the put-call parity
}

heston_call(S_0, K, r, T, kappa, theta, xi, rho, v_0)
heston_put(S_0, K, r, T, kappa, theta, xi, rho, v_0)

