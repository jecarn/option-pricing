asset <- function(S_0, K, r, vol, T){
  self <- list(S_0 = S_0, K = K, r = r, vol = vol, T = T)
  class(self) <- "Asset"
}

a <- asset(100, 100, 0.05, 0.2, 1)
a
a$r


library(matrixStats)
library(ggplot2)


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

PriceDensity <- function(){
  structure(
    list(),
    class = "DensityLayer"
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
  
  V <- matrix(object$params$vol, ncol = (layer$N + 1), nrow = layer$R)
  S <- matrix(object$params$S_0, ncol = (layer$N + 1), nrow = layer$R)
  
  for (n in 1 : layer$N){
    V_n <- V[, n]
    V[, n + 1] <- pmax(V_n + kappa * (theta - V_n) * dt + xi * sqrt(V_n) * dW_V[, n], 0)
    
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

add_layer.DensityLayer <- function(object, layer){
  df <- data.frame(x = object$asset_paths[, object$layers$N + 1])
  
  density_plot <- ggplot(df, aes(x = x)) +
    geom_vline(xintercept = object$params$S_0, linetype = "dashed", color = "darkgray") +
    geom_density(color = "darkgray", fill = "cornflowerblue", alpha = 0.4) +
    labs(x = "Terminal Price", y = "Density", title = "Simulation Price Density") +
    theme_minimal()
  
  print(density_plot)
  
  object
}

a <- Asset(S_0 = 100, r = 0.05, vol = 0.2, T = 1) +
  Heston(kappa = 1, theta = 0.2, xi = 0.5, rho = -0.7) +
  Merton(lambda = 1, gamma = 0,  delta = 0.1) +
  Simulation(R = 10000, N = 250) +
  PriceDensity()

a <- Asset(S_0 = 100, r = 0.05, vol = 0.2, T = 1)
a <- a + Heston(1, 0.2, 0.5, -0.7)
a <- a + Merton(1, 0, 0.1)
a <- a + Simulation(R = 1000, N = 10)
a <- a + PriceDensity()


ggplot(data.frame(x = a$asset_paths[, a$layers$N + 1]), aes(x = x)) +
  geom_vline(xintercept = a$params$S_0, linetype = "dashed", color = "darkgray") +
  geom_density(color = "darkgray", fill = "cornflowerblue", alpha = 0.4) +
  labs(x = "Terminal Price", y = "Density", title = "Simulated Price Density") +
  theme_minimal()

