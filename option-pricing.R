library(matrixStats)
library(ggplot2)

Asset <- function(S_0, r, vol, T){
  structure(
    list(
      params = list(S_0 = S_0, r = r, vol = vol, T = T),
      models = list(Heston = 0, Merton = 0),
      layers = list(),
      diagnostics = list(start_time = Sys.time())
    ),
    class = "Asset"
  )
}

Simulation <- function(R, N){
  structure(
    list(R = as.numeric(R), N = as.numeric(N)),
    class = "SimulationLayer"
  )
}

TerminalSimulation <- function(R){
  structure(
    list(R = as.numeric(R)),
    class = "TerminalSimulationLayer"
  )
}

HestonStochasticVol <- function(kappa, theta, xi, rho){
  structure(
    list(
      kappa = as.numeric(kappa),
      theta = as.numeric(theta),
      xi = as.numeric(xi),
      rho = as.numeric(rho)
    ),
    class = "HestonLayer"
  )
}

MertonJumpDiffusion <- function(lambda, gamma, delta){
  structure(
    list(
      lambda = as.numeric(lambda),
      gamma = as.numeric(gamma),
      delta = as.numeric(delta)
    ),
    class = "MertonLayer"
  )
}

PriceDensity <- function(){
  structure(
    list(),
    class = "DensityLayer"
  )
}

PricePlot <- function(n_plots = 15){
  structure(
    list(n_plots = as.numeric(n_plots)),
    class = "PricePlotLayer"
  )
}

EuropeanVanillaMC <- function(K){
  structure(
    list(K = as.numeric(K)),
    class = "EuroVanillaMCLayer"
  )
}

EuropeanVanillaBS <- function(K){
  structure(
    list(K = as.numeric(K)),
    class = "EuroVanillaBSLayer"
  )
}

AmericanVanillaCRR <- function(K, N = 250){
  structure(
    list(K = as.numeric(K), N = as.numeric(N)),
    class = "AmerVanillaCRRLayer"
  )
}

AmericanVanillaLSM <- function(K){
  structure(
    list(K = as.numeric(K)),
    class = "AmerVanillaLSMLayer"
  )
}

TimePerformance <- function(){
  structure(
    list(),
    class = "TimePerformanceLayer"
  )
}

`+.Asset` <- function(object, layer){
  UseMethod("add_layer", layer)
}

add_layer.HestonLayer <- function(object, layer){
  object$models$Heston = 1
  object$layers <- c(object$layers, layer)
  object
}

add_layer.MertonLayer <- function(object, layer){
  object$models$Merton = 1
  object$layers <- c(object$layers, layer)
  object
}

add_layer.SimulationLayer <- function(object, layer){
  object$layers <- c(object$layers, layer)
  
  dt <- object$params$T / object$layers$N
  
  # Asset Brownian increments
  Z_1 <- matrix(rnorm(layer$R * layer$N), nrow = layer$R, ncol = layer$N)
  dW_S <- Z_1 * sqrt(dt)
  
  # Heston model
  if (object$models$Heston == 1){
    Z_2 <- matrix(rnorm(layer$R * layer$N), nrow = layer$R, ncol = layer$N)
    dW_V <- object$layers$rho * Z_1 * sqrt(dt) + sqrt(1 - object$layers$rho^2) * Z_2 * sqrt(dt)
    
    kappa <- object$layers$kappa
    theta <- object$layers$theta
    xi <- object$layers$xi
  } else {
    dW_V <- matrix(0, nrow = layer$R, ncol = layer$N)
    
    kappa <- 0
    theta <- 0
    xi <- 0
  }
  
  # Merton model
  if (object$models$Merton == 1){
    k_bar <- exp(object$layers$gamma + object$layers$delta^2 / 2) - 1 # expected jump compensation
    
    dq <- matrix(
      ifelse(rpois(layer$R * layer$N, object$layers$lambda * dt) > 0, 1, 0),
      nrow = layer$R, ncol = layer$N, byrow = TRUE
    ) # Poisson occurrences dq_t
    
    k <- matrix(
      exp(rnorm(layer$R * layer$N, object$layers$gamma, object$layers$delta)) - 1,
      nrow = layer$R, ncol = layer$N, byrow = TRUE
    ) # random normal jump sizes
    
    lambda <- object$layers$lambda
  } else {
    lambda <- 0
    k_bar <- 0
    dq <- matrix(0, nrow = layer$R, ncol = layer$N)
    k <- matrix(0, nrow = layer$R, ncol = layer$N)
  }
  
  V <- matrix(object$params$vol^2, ncol = (layer$N + 1), nrow = layer$R) # variance matrix
  S <- matrix(object$params$S_0, ncol = (layer$N + 1), nrow = layer$R) # asset price path matrix
  
  for (n in 1 : layer$N){
    V_n <- pmax(V[, n], 0)
    V[, n + 1] <- V_n + kappa * (theta - pmax(V_n, 0)) * dt + xi * sqrt(pmax(V_n, 0)) * dW_V[, n]
    # V[, n + 1] <- pmax(V_n + kappa * (theta - V_n) * dt + xi * sqrt(V_n) * dW_V[, n], 0) # incorrect
    
    S_n <- S[, n]
    S[, n + 1] <- S_n +
      (object$params$r - lambda * k_bar) * S_n * dt + # drift
      S_n * sqrt(V_n) * dW_S[, n] + # diffusion
      S_n * k[, n] * dq[, n] # jump
    # S[, n + 1] <- S_n * exp((r - V_n / 2) * dt + sqrt(V_n) * dW_S[, n]) # alternative discretization
  }
  
  object$asset_paths <- S
  object$vol_paths <- V
  
  object
}

add_layer.TerminalSimulationLayer <- function(object, layer){
  object$layers <- c(object$layers, layer)
  
  # Retrieving parameters
  S_0 <- object$params$S_0
  vol <- object$params$vol
  r <- object$params$r
  T <- object$params$T
  
  R <- layer$R
  
  # Diffusion
  Z <- rnorm(R)
  
  # Merton model
  if (object$models$Merton == 1){
    lambda <- object$layers$lambda
    gamma <- object$layers$gamma
    delta <- object$layers$delta
    
    k_bar <- exp(gamma + delta^2 / 2) - 1 # expected jump compensation
    
    q <- rpois(layer$R, lambda * T)
    
    jumps_mx <- cbind(q, 1, exp(matrix(rnorm(R * max(q), gamma, delta), nrow = R, ncol = max(q))))
    
    jump_func <- function(mx_row, output){
      q <- mx_row[1]
      return(prod(mx_row[2 : (q + 2)]))
    }
    
    J <- apply(jumps_mx, 1, jump_func)
    
  } else {
    k_bar <- 0
    J <- 1
  }
  
  S <- S_0 * exp((r - lambda * k_bar - vol^2 / 2) * T + vol * Z) * J # terminal asset price vector
  
  object$asset_paths <- cbind(S_0, S)
  
  object
}

add_layer.DensityLayer <- function(object, layer){
  df <- data.frame(x = object$asset_paths[, object$layers$N + 1])
  
  density_plot <- ggplot(df, aes(x = x)) +
    geom_vline(xintercept = object$params$S_0, linetype = "dashed", color = "darkgray") +
    geom_density(color = "grey15", fill = "cornflowerblue", alpha = 0.4) +
    labs(x = "Terminal Price", y = "Density") +
    theme_bw()
  
  print(density_plot)
  
  object
}

add_layer.EuroVanillaMCLayer <- function(object, layer){
  object$layers <- c(object$layers, layer)
  
  # Retrieving parameters
  S_0 <- object$params$S_0
  vol <- object$params$vol
  r <- object$params$r
  T <- object$params$T
  
  K <- layer$K
  
  terminal_prices <- object$asset_paths[, ncol(object$asset_paths)]
  
  object$euro_vanilla_mc_call <- exp(-r * T) * mean(pmax(terminal_prices - K, 0))
  # object$euro_vanilla_mc_put <- exp(-object$params$r * object$params$T) * mean(pmax(layer$K - terminal_prices, 0))
  object$euro_vanilla_mc_put <- object$euro_vanilla_mc_call + exp(-r * T) * K - S_0
  
  object
}

add_layer.EuroVanillaBSLayer <- function(object, layer){
  object$layers <- c(object$layers, layer)
  
  # Retrieving parameters
  S_0 <- object$params$S_0
  vol <- object$params$vol
  r <- object$params$r
  T <- object$params$T
  
  K <- layer$K
  
  d_1 <- (log(S_0 / K) + (r + vol^2 / 2) * T) / (vol * sqrt(T))
  d_2 <- d_1 - vol * sqrt(T)
  
  object$euro_vanilla_bs_call <- S_0 * pnorm(d_1) - K * exp(-r * T) * pnorm(d_2)
  object$euro_vanilla_bs_put <- object$euro_vanilla_bs_call + K * exp(-r * T) - S_0
  
  object
}

add_layer.AmerVanillaCRRLayer <- function(object, layer){
  object$layers <- c(object$layers, layer)
  
  # Retrieving parameters
  S_0 <- object$params$S_0
  vol <- object$params$vol
  r <- object$params$r
  T <- object$params$T
  
  K <- layer$K
  N <- layer$N
  
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
  
  exp_disc_payoff <- function(p, q, r, dt, future_payoff){ # discounts expected future payoffs
    m <- length(future_payoff)
    payoff_vec <- c()
    for (j in seq(m - 1)){
      payoff_vec <- c(payoff_vec, exp(-r * dt) * (q * future_payoff[j] + p * future_payoff[j + 1]))
    }
    return(payoff_vec)
  }
  
  # Calculations
  future_payoff <- pmax(K - prices_at_t(S_0, u, d, N), 0)
  
  for (j in rev(seq(N - 1))){ # iterating over time steps
    exp_disc_poff <- exp_disc_payoff(p, q, r, dt, future_payoff)
    current_prices <- prices_at_t(S_0, u, d, j)
    future_payoff <- pmax(exp_disc_poff, pmax(K - current_prices, 0))
  }
  
  object$amer_vanilla_crr_put <- exp_disc_payoff(p, q, r, dt, future_payoff)
  
  object
}

add_layer.AmerVanillaLSMLayer <- function(object, layer){
  object$layers <- c(object$layers, layer)
  
  # Retrieving parameters
  S_0 <- object$params$S_0
  vol <- object$params$vol
  r <- object$params$r
  T <- object$params$T
  
  K <- layer$K
  
  R <- object$layers$R
  N <- object$layers$N
  
  dt <- T / N # time increment
  
  X <- object$asset_paths[, -1]
  
  exercise_matrix <- pmax(K - X, 0) # matrix of immediate exercise values
  
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
  
  object$amer_vanilla_lsm_put <- max(sum(disc_cf_matrix) / R, K - S_0) # put value
  
  object
}

add_layer.TimePerformanceLayer <- function(object, layer){
  object$layers <- c(object$layers, layer)
  
  print(Sys.time() - object$diagnostics$start_time)
  
  object
}

add_layer.PricePlotLayer <- function(object, layer){
  object$layers <- c(object$layers, layer)
  
  price_paths <- t(object$asset_paths[1 : layer$n_plots, ])
  
  matplot(
    price_paths,
    xlab = "Timestep (n)", ylab = "Asset Price",
    type = "l", lwd = 1.5, lty = 1
  )
  
  object
}

a <- Asset(S_0 = 100, r = 0.05, vol = 0.2, T = 1) +
  Merton(lambda = 1, gamma = 0,  delta = 0.1) +
  Simulation(R = 10000, N = 250) +
  PricePlot() +
  PriceDensity() +
  TimePerformance()


a <- Asset(S_0 = 100, r = 0.05, vol = 0.2, T = 1) +
  Simulation(R = 50000, N = 250) +
  AmericanVanillaCRR(K = 100) +
  AmericanVanillaLSM(K = 100) +
  TimePerformance()

a <- Asset(S_0 = 100, r = 0.05, vol = 0.2, T = 1) +
  Merton(lambda = 1, gamma = 0,  delta = 0.1) +
  Simulation(R = 10000, N = 250) +
  PricePlot() +
  TimePerformance()

a




a <- Asset(S_0 = 100, r = 0.05, vol = 0.2, T = 1) +
  Heston(kappa = 1, theta = 0.2, xi = 0.5, rho = -0.7) +
  Merton(lambda = 1, gamma = 0,  delta = 0.1) +
  Simulation(R = 100000, N = 1) +
  PriceDensity() +
  EuropeanVanillaMC(K = 100) +
  EuropeanVanillaBS(K = 100)

a <- Asset(S_0 = 100, r = 0.05, vol = 0.2, T = 1) +
  Simulation(R = 100000, N = 50) +
  PriceDensity() +
  EuropeanVanillaMC(K = 100) +
  EuropeanVanillaBS(K = 100)

a <- Asset(S_0 = 100, r = 0.05, vol = 0.2, T = 1) +
  Merton(lambda = 1, gamma = 0,  delta = 0.1) +
  TerminalSimulation(R = 1000) +
  EuropeanVanillaBS(K = 100) +
  EuropeanVanillaMC(K = 100) +
  AmericanVanillaCRR(K = 100)






