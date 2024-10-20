estimate <- function(X, alpha, beta, Vs, Vt, a, b, c, d, eps1s, eps1t, eps2s, eps2t, r) {
    N <- nrow(X)
    eta1 <- X %*% alpha + Vs %*% a + Vt %*% b + Vs %*% eps1s + Vt %*% eps1t
    pi <- pmax(0.01, pmin(0.99, exp(eta1) / (1 + exp(eta1)))) # at-risk probability
    u <- rbinom(N, 1, pi) # at-risk indicator
    if (ncol(X) == 1) {
        eta2 <- X[u == 1, ] * beta + Vs[u == 1, ] %*% c + Vt[u == 1, ] %*% d + Vs[u == 1, ] %*% eps2s + Vt[u == 1, ] %*% eps2t # Linear predictor for count part
    } else {
        eta2 <- X[u == 1, ] %*% beta + Vs[u == 1, ] %*% c + Vt[u == 1, ] %*% d + Vs[u == 1, ] %*% eps2s + Vt[u == 1, ] %*% eps2t # Linear predictor for count part
    }
    psi <- pmax(0.01, pmin(0.99, exp(eta2) / (1 + exp(eta2)))) # Prob of success

    mu <- r * psi / (1 - psi) # NB mean
    y <- rep(0, N) # Response
    y[u == 1] <- mu # If at risk, draw from NB
    return(y)
}
