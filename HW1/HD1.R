library(MASS)
library(rjags)
library(tools)

setwd("~/Github/STAE004F/hw1")

set.seed(1)

## Part (c)

# Hyperparams
Q     <- 53
n     <- 41
k     <- 2
sigma <- 2.51

# Read in data
data <- read.table("y_hw1.dat")
colnames(data) <- paste0("y", 1:Q)
y <- do.call(c, unname(as.vector(data)))

# Construct Z
Z <- matrix(nrow = n * Q, ncol = Q + 1)
Z[,1] <- 1

for (i in 2:(Q+1)) {
  Z[,i] <- c(
    rep(0, (i - 2) * n),
    rep(1, n),
    rep(0, (Q - i + 1) * n) 
  )
}

# tau vector
delta <- 1e-3
tau_vec <- seq(0.001, 3, by = delta)

log_tau_post <- sapply(tau_vec, function(tau) {
  # Construct Qx
  Qx <- diag(tau^-2, Q + 1)
  Qx[1, 1] <- 100^-2
  
  # Construct Q_{x \mid y}
  Qxy      <- diag(1, nrow = Q + 1, ncol = Q + 1)
  Qxy[,1]  <- Qxy[1,] <- 1
  Qxy[1,1] <- Q
  Qxy <- sigma^-2 * n * Qxy + Qx
  
  mu_xy <- solve(Qxy, sigma^-2 * t(Z) %*% y)
  log_tau <- -Q * log(tau) - 1/2 * log(det(Qxy)) + (sigma^-2/2 * t(mu_xy) %*% t(Z) %*% y - 1/2 * tau)
  return(log_tau)
})

log_tau_post <- log_tau_post - max(log_tau_post)
tau_post_unscaled <- exp(log_tau_post)
c <- sum(delta * tau_post_unscaled)^-1
tau_post <- c * tau_post_unscaled

pdf("tau_post.pdf")
plot(tau_vec, tau_post, type = "l", xlab = expression(tau), ylab = "Density")
dev.off()

sum(tau_post * delta, na.rm = TRUE)

## Part (e)

log_prior <- function(theta) { return(log(k^-1) - k^-1 * exp(theta) + theta) }
log_like  <- function(theta) {
  Qx <- diag(exp(theta)^-2, Q + 1)
  Qx[1, 1] <- 100^-2
  Qxy <- Qx + sigma^-2 * t(Z) %*% Z
  mu_xy <- solve(Qxy, sigma^-2 * t(Z) %*% y)

  return(-Q * theta - 1/2 * log(det(Qxy)) + (sigma^-2/2 * t(mu_xy) %*% t(Z) %*% y - 1/2 * exp(theta)))
}

log_post <- function(theta) { log_like(theta) + log_prior(theta) }

mcmc_simulation <- function(S, burnin) { 
  
  # initial values
  theta_mc <- x_mc <- rep(0, S)
  x_mc <- matrix(nrow = S, ncol = Q + 2)
  colnames(x_mc) <- c("eta", paste0("u", 1:53), "tau")
  x_mc[1,] <- c(mean(y), apply(data, 2, mean) - mean(y), tau_vec[which.max(tau_post)])
  xi <- 2.38

  for (s in 2:S) {
    theta_curr <- rnorm(1, theta_mc[s - 1], xi)
    alpha <- log_post(theta_curr) - log_post(theta_mc[s - 1])
    u <- log(runif(1))
   
    theta_mc[s] <- ifelse(u <= alpha, theta_curr, theta_mc[s - 1])

    tau <- exp(theta_mc[s])

    # Construct Qx
    Qx <- diag(tau^-2, Q + 1)
    Qx[1, 1] <- 100^-2

    Qxy <- Qx + sigma^-2 * t(Z) %*% Z
    
    R <- chol(Qxy)
    z <- mvrnorm(1, rep(0, Q + 1), diag(1, Q + 1))
    u <- solve(R, z)
    b <- sigma^-2 * t(Z) %*% y
    v <- solve(t(R), b) 
    mu_xy <- solve(R, v)
    x_ast <- mu_xy + u

    x_mc[s,] <- c(x_ast, tau)
  }

  sel_mc <- x_mc[(burnin + 1):S,]
  return(mcmc(sel_mc))

}

## Part (f)
chains <- lapply(1:4, function(i) mcmc_simulation(13000, 3000))
chains <- mcmc.list(chains)

sum_chains <- summary(chains)

# draw the posterior and prior densities of tau on the same graph
tau_prior <- dexp(tau_vec, k^-1)

pdf("tau_densities.pdf")
plot(tau_vec, tau_prior, lwd = 1.5, type = "l", xlab = expression(tau), ylab = "Density", ylim = range(0, 4.5), lty = 2)
lines(tau_vec, tau_post, lwd = 1.5)
legend("topright", inset = c(0.025, 0.025), bty = "n", lty = c(2, 1), legend = c("Prior", "Posterior"))
lines(density(runjags::combine.mcmc(chains)[,"tau"], adjust = 3), lty = 2, col = "red")
dev.off()

# present random effects in a figure (posterior median + 95% CI vs. cat. no.)
pdf("random_effects.pdf")
plot(sum_chains$quantiles[2:(Q+1), 3], ylim = range(sum_chains$quantiles[2:(Q+1),]), xlab = "Category", ylab = "Effect size")
arrows(x0 = 1:Q, y0 = sum_chains$quantiles[2:(Q+1), 1], y1 = sum_chains$quantiles[2:(Q+1), 5], code = 3, angle = 90, length = 0.05)
abline(h = 0, lty = 2)
dev.off()
