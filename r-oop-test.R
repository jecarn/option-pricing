asset <- function(S_0, K, r, vol, T){
  self <- list(S_0 = S_0, K = K, r = r, vol = vol, T = T)
  class(self) <- "Asset"
}

a <- asset(100, 100, 0.05, 0.2, 1)
a
a$r


library(matrixStats)


Asset <- function(S_0, r, vol, T){
  structure(
    list(
      params = list(S_0 = S_0, r = r, vol = vol, T = T),
      models = list(Heston = 0, Merton = 0),
      layers = list()
    ),
    class = "Asset"
  )
}

Simulation <- function(R, N = 100){
  structure(
    list(R = as.numeric(R), N = as.numeric(N)),
    class = "SimulationLayer"
  )
}

Heston <- function(kappa, theta, xi, rho){
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

Merton <- function(lambda, gamma, delta){
  structure(
    list(
      lambda = as.numeric(lambda),
      gamma = as.numeric(gamma),
      delta = as.numeric(delta)
    ),
    class = "MertonLayer"
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
  object$layer <- c(object$layers, layer)
  object
}

add_layer.SimulationLayer <- function(object, layer){
  object$layers <- c(object$layers, layer)
  
  # Asset Brownian increments
  Z_1 <- matrix(rnorm(layer$R * layer$N), nrow = layer$R, ncol = layer$N)
  dW_S <- Z_1 * sqrt(dt)
  
  # Volatility Brownian increments
  if (object$models$Heston == 1){
    Z_2 <- matrix(rnorm(layer$R * layer$N), nrow = layer$R, ncol = layer$N)
    dW_V <- object$layers$rho * Z_1 * sqrt(dt) + sqrt(1 - object$layers$rho^2) * Z_2 * sqrt(dt)
  } else {
    V <- object$params$vol
  }
  
  V <- matrix(v_0, ncol = (N + 1), nrow = R)
  S <- matrix(object$params$S_0, ncol = (object$layers$N + 1), nrow = object$layers$R)
  
  for (n in 1 : object$layers$N){
    V_n <- V[, n]
    V[, n + 1] <- pmax(V_n + kappa * (theta - V_n) * dt + xi * sqrt(V_n) * dW_V[, n], 0)
    
    S_n <- S[, n]
    S[, n + 1] <- S_n + r * S_n * dt + S_n * sqrt(V_n) * dW_S[, n]
    # S[, n + 1] <- S_n * exp((r - V_n / 2) * dt + sqrt(V_n) * dW_S[, n]) # alternative discretization
  }
  
  X <- S[, -1] # asset prices
  V <- V[, -1] # variance paths
  
  
  object$vol_paths <- matrix(
    rnorm(layer$R * layer$N), nrow = layer$R, ncol = layer$N
  )
  object$asset_paths <- matrix()
  
  object
}

a <- Asset(100, 0.05, 0.2, 1)
a <- a + Simulation(R = 5, N = 100)

a <- a + Heston(1, 0.2, 0.5, -0.7)
a <- a + Merton(1, 0, 0.1)


