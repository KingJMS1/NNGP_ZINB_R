library(mvtnorm)
library(Matrix)
library(ZINB.GP)

# This file can be run via `Rscript example_usage.R`

# Generate a number of samples at a number of spatial locations, return design matrices and total number of observations, avg_obs describes how many observations there area at each spatiotemporal point
make_Vs_Vt <- function(num_spatial, num_temporal, avg_obs) {
    n_time_points <- num_temporal # Number of temporal units
    n_unit_mat <- matrix(rpois(num_spatial * n_time_points, avg_obs), nrow = num_spatial, byrow = TRUE) # sample around 2 observations per sampling unit (both space and time)

    N <- sum(n_unit_mat)
    id <- c()
    for (i in seq_len(nrow(n_unit_mat))) {
        id <- c(id, rep(i, sum(n_unit_mat[i, ])))
    }

    tp_seq <- c()
    for (j in seq_len(ncol(n_unit_mat))) {
        tp_seq <- c(tp_seq, rep(j, sum(n_unit_mat[, j])))
    }

    # spatial design matrix
    Vs <- as.matrix(sparseMatrix(i = 1:N, j = id, x = rep(1, N)))

    # temporal design matrix
    Vt <- as.matrix(sparseMatrix(i = 1:N, j = tp_seq, x = rep(1, N)))

    return(list(Vs = Vs, Vt = Vt, N = N))
}

make_spatial_effects <- function(phi_nb, phi_bin, sigma_bin_s, sigma_nb_s, coords) {
    ##########################
    # Spatial Random Effects #
    ##########################
    Ds <- as.matrix(dist(coords))
    Ks_bin <- sigma_bin_s^2 * exp(-phi_bin * Ds)
    Ks_nb <- sigma_nb_s^2 * exp(-phi_nb * Ds)
    a <- t(rmvnorm(n = 1, sigma = Ks_bin))
    c <- t(rmvnorm(n = 1, sigma = Ks_nb))

    return(list(a = a, c = c, Ds = Ds))
}

make_temporal_effects <- function(l1t, l2t, sigma1t, sigma2t, n_time_points) {
    ###########################
    # Temporal Random Effects #
    ###########################
    w <- matrix(1:n_time_points, ncol = 1)
    Dt <- as.matrix(dist(w))
    
    Kt_bin <- sigma1t^2 * exp(-Dt / (l1t^2))
    Kt_nb <- sigma2t^2 * exp(-Dt / (l2t^2))
    
    b <- t(rmvnorm(n = 1, sigma = Kt_bin))
    d <- t(rmvnorm(n = 1, sigma = Kt_nb))
    
    return(list(b = b, d = d, Dt = Dt))
}

#################
# Generate Data #
#################
num_spatial <- 30
num_temporal <- 10

# Get Spatial and temporal design matrices, and total number of observations
out <- make_Vs_Vt(num_spatial, num_temporal, 2)
Vs <- out$Vs
Vt <- out$Vt
N  <- out$N

coords <- cbind(runif(num_spatial), runif(num_spatial))
x <- rnorm(N, 0, 1)
X <- as.matrix(x) # Design matrix, can add additional covariates (e.g., race, age, gender)
X <- cbind(1, X)
p <- ncol(X)


phi_nb <- 1
phi_bin <- 2
sigma_bin_s <- 1
sigma_nb_s <- 1
out <- make_spatial_effects(phi_nb, phi_bin, sigma_bin_s, sigma_nb_s, coords)
a <- out$a
c <- out$c
Ds <- out$Ds

l1t <- 2
l2t <- 3
sigma1t <- 0.5
sigma2t <- 0.5
out <- make_temporal_effects(l1t, l2t, sigma1t, sigma2t, num_temporal)
b <- out$b
d <- out$d
Dt <- out$Dt

#################
# Fixed Effects #
#################
# Binomial Part
alpha <- c(-0.25, 0.25)

# Count Part
beta <- c(.5, -.25)


################
# Random Noise #
################
sigma_eps1s <- sigma_eps2s <- 0.05
eps1s <- t(rmvnorm(n = 1, sigma = diag(sigma_eps1s^2, nrow = num_spatial)))
eps2s <- t(rmvnorm(n = 1, sigma = diag(sigma_eps2s^2, nrow = num_spatial)))

sigma_eps1t <- true.sigma_eps2t <- sigma_eps2t <- 0.05
eps1t <- t(rmvnorm(n = 1, sigma = diag(sigma_eps1t^2, nrow = num_temporal)))
eps2t <- t(rmvnorm(n = 1, sigma = diag(sigma_eps2t^2, nrow = num_temporal)))

#######################
# Binomial Simulation #
#######################
true.phi1 <- phi1 <- Vs %*% a + Vt %*% b
eta1 <- X %*% alpha + phi1 + Vs %*% eps1s + Vt %*% eps1t

p_at_risk <- exp(eta1) / (1 + exp(eta1)) # 1-pr("structural zero")
u <- rbinom(N, 1, p_at_risk[, 1]) # at-risk indicator

#################
# NB Simulation #
#################
phi3 <- Vs %*% c + Vt %*% d
eta2 <- X[u == 1, ] %*% beta + phi3[u == 1, ] + Vs[u == 1, ] %*% eps2s + Vt[u == 1, ] %*% eps2t # Linear predictor for count part
N1 <- sum(u == 1)

r <- 1 # NB dispersion
psi <- exp(eta2) / (1 + exp(eta2)) # Prob of success

mu <- r * psi / (1 - psi) # NB mean
y <- rep(0, N) # Response
y[u == 1] <- rnbinom(N1, r, mu = mu[, 1]) # If at risk, draw from NB


#################
# Run the Model #
#################
# Run for a short time for demo purposes
cat("Running Model\n")
output <- ZINB_NNGP(X, y, coords, Vs, Vt, Ds, Dt, M = 10, 500, 100, 1, TRUE)
predictions <- output$Y_pred
sim_alpha <- output$Alpha
sim_beta <- output$Beta

# Examine coefficients for regressions
cat("\nLogistic Regression Coefficients:\n")
alpha
apply(sim_alpha, 2, function(x) quantile(x, probs=c(0.05, 0.95)))

cat("\nNegative Binomial Coefficients:\n")
beta
apply(sim_beta, 2, function(x) quantile(x, probs=c(0.05, 0.95)))

# Examine how often various samples are at risk
cat("\nAt risk 'probabilities':\n")
at_risk <- output$at_risk
sim_p_at_risk <- apply(at_risk, 2, mean)
sim_p_at_risk[1:20]
