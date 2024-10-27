library(mvtnorm)
library(Matrix)

##################
# Parameter Setting #
##################
n <- 200 # Number of spatial units
n_time_points <- 20
n_unit_mat <- matrix(rpois(n * n_time_points, 2), nrow = n, byrow = TRUE) # number of obs in each sampling unit

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
Vs[1:15, 1:10]

# temporal design matrix
Vt <- as.matrix(sparseMatrix(i = 1:N, j = tp_seq, x = rep(1, N)))
Vt[1:15, 1:10]

##################
# Generate Data  #
##################
coords <- cbind(runif(n), runif(n))
x <- rnorm(N, 0, 1)
boxplot(x)
X <- cbind(1, x) # Design matrix, can add additional covariates (e.g., race, age, gender)
p <- ncol(X)


#################
# Fixed Effects #
#################
# Binomial Part
true.alpha <- alpha <- c(-0.25, 0.25)

# Count Part
true.beta <- beta <- c(.5, -.25)

#################
# Spatial Random Effects #
#################
true.phi_bin <- phi_bin <- true.phi_nb <- phi_nb <- 8
true.sigma1s <- sigma1s <- true.sigma2s <- sigma2s <- 0.5
Ds <- as.matrix(dist(coords))
Ks_bin <- sigma1s^2 * exp(-phi_bin * Ds)
Ks_nb <- sigma2s^2 * exp(-phi_nb * Ds)
true.a <- a <- t(rmvnorm(n = 1, sigma = Ks_bin))
true.c <- c <- t(rmvnorm(n = 1, sigma = Ks_nb))

#################
# Temporal Random Effects #
#################
w <- matrix(1:n_time_points, ncol = 1)
true.l1t <- l1t <- true.l2t <- l2t <- 1
true.sigma1t <- sigma1t <- true.sigma2t <- sigma2t <- 0.2
Dt <- as.matrix(dist(w))
Kt_bin <- sigma1t^2 * exp(-Dt / (l1t^2))
Kt_nb <- sigma2t^2 * exp(-Dt / (l2t^2))
true.b <- b <- t(rmvnorm(n = 1, sigma = Kt_bin))
true.d <- d <- t(rmvnorm(n = 1, sigma = Kt_nb))


#################
# Random Noise #
#################
true.sigma_eps1s <- sigma_eps1s <- true.sigma_eps2s <- sigma_eps2s <- 0.05
true.eps1s <- eps1s <- t(rmvnorm(n = 1, sigma = diag(sigma_eps1s^2, nrow = n)))
true.eps2s <- eps2s <- t(rmvnorm(n = 1, sigma = diag(sigma_eps2s^2, nrow = n)))
true.sigma_eps1t <- sigma_eps1t <- true.sigma_eps2t <- sigma_eps2t <- 0.05
true.eps1t <- eps1t <- t(rmvnorm(n = 1, sigma = diag(sigma_eps1t^2, nrow = n_time_points)))
true.eps2t <- eps2t <- t(rmvnorm(n = 1, sigma = diag(sigma_eps2t^2, nrow = n_time_points)))
summary(true.eps1t)
summary(true.eps2t)

#################
# Binomial Simulation #
#################
# combine the spatial and temporal effects
V <- cbind(Vs, Vt)
true.phi1 <- phi1 <- Vs %*% a + Vt %*% b
eta1 <- X %*% alpha + phi1 + Vs %*% eps1s + Vt %*% eps1t

plot(X %*% alpha, eta1)
boxplot(phi1[, 1])

p_at_risk <- exp(eta1) / (1 + exp(eta1)) # 1-pr("structural zero")

# plot for the binomial part
par(mfrow = c(2, 2))
boxplot(X %*% alpha, main = "x*alpha")
boxplot(phi1[, 1], main = "phi1")
boxplot(eta1[, 1], main = "eta1=X%*%alpha+phi1")
boxplot(p_at_risk[, 1], main = "pi=exp(eta1)/(1+exp(eta1))")

u <- rbinom(N, 1, p_at_risk[, 1]) # at-risk indicator
N1 <- sum(u)
pstruct0 <- 1 - mean(u) # Proportion of structural zeros


tmp <- u
tmp[u == 0] <- 0 # If in structural group then set tmp to 0
nis1 <- tapply(tmp, id, sum) # Number of at risk observations in each county

#################
# NB Simulation #
#################
true.phi3 <- phi3 <- Vs %*% c + Vt %*% d
eta2 <- X[u == 1, ] %*% beta + phi3[u == 1, ] + Vs[u == 1, ] %*% eps2s + Vt[u == 1, ] %*% eps2t # Linear predictor for count part
N1 <- sum(u == 1)

r <- 1 # NB dispersion
psi <- exp(eta2) / (1 + exp(eta2)) # Prob of success

mu <- r * psi / (1 - psi) # NB mean
y <- rep(0, N) # Response
y[u == 1] <- rnbinom(N1, r, mu = mu[, 1]) # If at risk, draw from NB
pzero <- length(y[y == 0]) / N # Proportion of zeros

par(mfrow = c(2, 2))
boxplot(X %*% beta, main = "X[u==1,]%*%beta")
boxplot(eta2[, 1], main = "eta2=X[u==1,]%*%beta+phi3[u==1]")
boxplot(psi[, 1], main = "psi=exp(eta2)/(1+exp(eta2))")
boxplot(mu[, 1], main = "mu=r*psi/(1-psi)")

tmp <- table(y) / N * 100 # convert to %s (divide by N multiply by 100)
tmp[1:5]
barplot(tmp, ylab = "Percent", xlab = "Count", col = "lightblue")

save(
    list = ls(all.names = TRUE),
    file = "./simdata/simdata_dynamic.RData",
    envir = .GlobalEnv
)
