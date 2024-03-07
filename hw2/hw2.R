library(latex2exp)
library(MASS)
library(Matrix)
library(lattice)
library(rSPDE)
library(plot3D)
library(rjags)

set.seed(42)

setwd("~/Github/STAE004F/hw2")

dat <- read.table("qh_hw2.dat")
colnames(dat) <- c("id", "Q", "B", "h", "S", "T")

n <- nrow(dat)
p <- ncol(dat)

dat$A <- with(dat, B * h)
dat$P <- with(dat, B + 2 * h)
dat$R <- with(dat, A / P)

## part (a)
lm.kappa <- lm(Q ~ I(R^(2/3) * A * S^(1/2)), data = dat)
kappa    <- unname(coef(lm.kappa)[2])^-1 # estimate Manning roughness coef.

pdf("1a_Q_vs_h.pdf", width = 12)
par(mfrow = c(1, 2))
plot(Q ~ h, data = dat)
with(dat, lines(h, kappa^-1 * R^(2/3) * A * S^(1/2), lty = 2))
with(dat, plot(kappa^-1 * R^(2/3) * A * S^(1/2), Q, xlab = TeX("$\\hat{Q} = \\kappa^{-1} R^{2/3} A S^{1/2}$")))
abline(a = 0, b = 1, lty = 2)
dev.off()

## part (b)

dat$y <- with(dat, log(Q) - 2/3 * log(R) - log(A) - 1/2 * log(S))

my <- mean(dat$y)
sy <- sd(dat$y)

pdf("1b_y_vs_h.pdf")
plot(y ~ h, data = dat)
abline(h = mean(dat$y))
abline(h = my - qt(0.975, df = n - 1) * sy / sqrt(n), lty = 2)
abline(h = my + qt(0.975, df = n - 1) * sy / sqrt(n), lty = 2)
dev.off()

## part (e)
y  <- dat$y
h  <- dat$h
m  <- 15
nu <- 2.5
h_p <- seq(0.000, 0.028, length.out = m)

mu_eta <- 4.565
mu_su  <- 0.10
mu_se  <- 0.05

s_eta <- 0.215
phi_u <- 0.16

mu_x <- c(mu_eta, rep(0, n))

# generate matern covariance matrix
R <- R22 <- R21 <- matrix(0, n, n)

for (i in 1:n) {
  for (j in 1:n) {
    R[i, j] <- (1 + (sqrt(20) * abs(h[i] - h[j])) / phi_u + (20 * (h[i] - h[j])^2) / (3 * phi_u^2)) * exp(-(sqrt(20) * abs(h[i] - h[j])) / phi_u)
    R21[i, j] <- (1 + (sqrt(20) * abs(h[i] - h_p[j])) / phi_u + (20 * (h[i] - h_p[j])^2) / (3 * phi_u^2)) * exp(-(sqrt(20) * abs(h[i] - h_p[j])) / phi_u)
    R22[i, j] <- (1 + (sqrt(20) * abs(h_p[i] - h_p[j])) / phi_u + (20 * (h_p[i] - h_p[j])^2) / (3 * phi_u^2)) * exp(-(sqrt(20) * abs(h_p[i] - h_p[j])) / phi_u)
  }
}

R21 <- t(R21[1:n, 1:15])
R22 <- R22[1:15, 1:15]

log_tau_post <- function(s_u, s_e) {
  cov_y <- s_eta^2 * matrix(1, n, n) + s_u^2 * R + s_e^2 * diag(n)
  mu_y  <- rep(mu_eta, n)

  prior <- - mu_su^-1 * s_u - mu_se^-1 * s_e
  like  <- - 1/2 * log(det(cov_y)) - 1/2 * t(y - mu_y) %*% solve(cov_y) %*% (y - mu_y)
  res   <- prior + like

  return(res)
}

n_grid  <- 1.5e+2
su_grid <- seq(0.0001, 0.16, length.out = n_grid)
se_grid <- seq(0.0001, 0.03, length.out = n_grid)
log_tau_post_grid <- matrix(0, n_grid, n_grid)

for (i in 1:n_grid) {
  for (j in 1:n_grid) {
    log_tau_post_grid[i, j] <- log_tau_post(su_grid[i], se_grid[j])
  }
}

mesh_ss <- mesh(su_grid, se_grid)
su_mesh <- mesh_ss$x
se_mesh <- mesh_ss$y
dsu <- (0.16 - 0.01) / n_grid
dse <- (0.03 - 0.01) / n_grid
dA <- dsu * dse

tau_post_grid <- exp(log_tau_post_grid - max(log_tau_post_grid))

Cint <- sum(tau_post_grid * dA)

tau_post_grid <- Cint^-1 * tau_post_grid

surf3D(su_mesh, se_mesh, tau_post_grid, box = TRUE, 
       bty = "b", colvar = tau_post_grid, colkey = TRUE,
       phi = 20, theta = 120, xlab = "σ_u", ylab = "σ_ε", zlab = "π(τ | y)")

## part (g)

log_post <- function(theta) {
  theta_u <- theta[1]
  theta_e <- theta[2]

  cov_y <- s_eta^2 * matrix(1, n, n) + exp(theta_u)^2 * R + exp(theta_e)^2 * diag(n)
  mu_y  <- rep(mu_eta, n)

  prior <- - mu_su^-1 * exp(theta_u) + theta_u - mu_se^-1 * exp(theta_e) + theta_e
  like  <- - 1/2 * log(det(cov_y)) - 1/2 * t(y - mu_y) %*% solve(cov_y) %*% (y - mu_y)

  return(prior + like)
}

mcmc_simulation <- function(S, burnin) {
  N <- S + burnin
  x_mc <- matrix(nrow = N, ncol = n + m + 3)
  theta_mc <- matrix(nrow = N, ncol = 2)
  colnames(x_mc) <- c("eta", paste0("u", 1:n), paste0("u", 1:m, "p"), "s_u", "s_e") 
  colnames(theta_mc) <- c("theta_u", "theta_e")

  theta_mc[1, ] <- c(-4.437, -4.726)

  x_mc[1,] <- c(mu_eta, rep(0, n + m), 0, 0)
  xi <- 2.38^2 / 2

  for (s in 2:N) {
    theta_p <- theta_mc[s - 1, ]
    theta_c <- mvrnorm(1, theta_p, xi * diag(2))
   
    alpha <- log_post(theta_c) - log_post(theta_p)
    u <- log(runif(1))

    theta_mc[s,] <- if (u < alpha) theta_c else theta_p

    tau_c <- exp(theta_mc[s,])
    names(tau_c) <- c("s_u", "s_e")
    s_u <- tau_c[1]; s_e <- tau_c[2]

    cov_x <- bdiag(s_eta^2, s_u^2 * R)

    A <- rbind(
      s_eta^2 * t(rep(1, n)),
      s_u^2 * R
    )

    B <- solve(s_eta^2 * matrix(1, n, n) + s_u^2 * R + s_e^2 * diag(n))

    mu_xy  <- mu_x + A %*% B %*% (y - mu_eta * rep(1, n)) 
    cov_xy <- cov_x - A %*% B %*% t(A)

    x_c <- mvrnorm(1, mu_xy, cov_xy)

    
    # sample from posterior predictive density
    mu_p  <- R21 %*% solve(R) %*% x_c[-1]
    cov_p <- R22 - R21 %*% solve(R) %*% t(R21)
    y_c <- mvrnorm(1, mu_p, cov_p) + mvrnorm(1, rep(0, m), s_e^2 * diag(m))
    
    x_mc[s, ] <- c(x_c, y_c, tau_c)
  }

  return(x_mc = x_mc[(burnin + 1):N,])
}

chains <- lapply(1:4, function(i) mcmc(mcmc_simulation(10000, 3000)))
save(chains, file = "chains_raw.Rdata")

mcmc.chains <- mcmc.list(chains)
save(mcmc.chains, file = "chains.Rdata")

gelman.diag(mcmc.chains)

## part (h)
sum.chains <- summary(mcmc.chains)
eta   <- sum.chains$quantiles[1, 3]
eta_l <- sum.chains$quantiles[1, 1]
eta_u <- sum.chains$quantiles[1, 5]
uh   <- sum.chains$quantiles[2:(n + 1), 3]
uh_l <- sum.chains$quantiles[2:(n + 1), 1]
uh_u <- sum.chains$quantiles[2:(n + 1), 5]

plot(h, eta + uh, type = "l", ylim = range(eta + uh_l - 0.025, eta + uh_u + 0.025), lwd = 2)
lines(h, eta + uh_l, lty = 2)
lines(h, eta + uh_u, lty = 2)
points(h, y)

## part (i)
u_p <- sum.chains$statistics[(n + 2):(n + m + 1), 1]
u_l <- sum.chains$quantiles[(n + 2):(n + m + 1), 1]
u_u <- sum.chains$quantiles[(n + 2):(n + m + 1), 5]

plot(h, eta + uh, xlim = c(0.000, 0.13), ylim = range(4, 5), col = "blue", type = "l", ylab = expression(eta + u(h[i])))
lines(h_p, eta + u_p, col = "red")
lines(h_p, eta_l + u_l, lty = 2, col = "red")
lines(h_p, eta_u + u_u, lty = 2, col = "red")
lines(h, eta_l + uh_l, lty = 2, col = "blue")
lines(h, eta_u + uh_u, lty = 2, col = "blue")
points(h, y)
legend("bottomright", bty = "n", legend = c("Predicted", "Fitted"), col = c("red", "blue"), lty = 1)
