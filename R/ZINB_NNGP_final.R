#' Run the ZINB NNGP model described in https://doi.org/10.1016/j.jspi.2023.106098.
#'
#' @param X Other Predictor variables
#' @param y Zero inflated count response
#' @param coords Spatial coordinates for NNGP
#' @param Vs   Spatially varying predictor variables (e.g. one-hot indication of which location this is for varying intercept), wrapped in sparseMatrix from Matrix R package. Will be multiplied by the spatial random effects for prediction.
#' @param Vt   Temporal varying predictor variables, wrapped in sparseMatrix from Matrix R package. Will be multiplied by the temporal random effects for prediction.
#' @param Ds   Spatial distance matrix, diagonal should be 0, off diagonal is distance between elements i and j in space, inputs to the spatial NNGP kernel
#' @param Dt   Temporal distance matirx, diagonal should be 0, off diagonal is distance between elements i and j in time, inputs to the temporal GP kernel
#' @param M    How many neighbors to allow in the spatial NNGP algorithm, defaults to 10.
#' @param nsim How long to run MCMC in total, must be greater than burn.
#' @param burn How long to run MCMC before saving samples.
#' @param thin How often to save MCMC samples, default is 1, saves every iteration.
#' @param save_ypred Whether or not to output the predicted values at every iteration
#' @param print_iter How often to print the iteration number of the MCMC chain.
#' @param print_progress Whether or not to print the iteration number of the MCMC chain.
#' @return A List of the following sampled values:          
#' Alpha, Beta, A, B, C, D, 
#' Eps1s, Eps2s, Eps1t, Eps2t, 
#' L1t, Sigma1t, L2t, Sigma2t, 
#' Phi_bin, Sigma1s, Phi_nb, Sigma2s, 
#' Sigma_eps1s, Sigma_eps2s, Sigma_eps1t, Sigma_eps2t, 
#' R, R2, 
#' at_risk, Y_pred if save_YPRED = TRUE
#' @export
#' @importFrom MASS glm.nb
#' @importFrom mvtnorm rmvnorm
#' @importFrom mvtnorm dmvnorm
#' @importFrom BayesLogit rpg
#' @importFrom FastGP rcppeigen_invert_matrix
#' @importFrom Matrix bdiag
#' @importFrom Matrix sparseMatrix
#' @importFrom msm rtnorm
#' @importFrom msm dtnorm
#' @importFrom stats dgamma
#' @importFrom stats rnorm
#' @importFrom LaplacesDemon rinvgamma
#' @importFrom stats runif
ZINB_NNGP <- function(X, y, coords, Vs, Vt, Ds, Dt, M = 10, nsim, burn, thin = 1, save_ypred = FALSE, print_iter = 100, print_progress = FALSE) {
    # TODO: Break down the Gibbs sampling and test all steps independently
    # TODO: Remove the need to compute Ds, Dt manually, take in coords for both instead so you can NNGP with large datasets

    # X is the design matrix with dimension N*p
    # x is the vector with length N
    # y is the count response with length N
    n <- nrow(coords) # number of clusters
    N <- nrow(X) # number of observations
    p <- ncol(X) # dimension of alpha and beta
    n_time_points <- ncol(Vt)


    ##########
    # Priors #
    ##########

    ####### priors for alpha and beta ######
    T0a <- diag(1, p)
    T0b <- diag(.01, p)
    s <- 0.0003

    ####### kernel hyperparameters  ######
    a_phi <- 8
    b_phi <- 1 # spatial kernel: phi~Gamma(shape=a,rate=b), mean=a/b
    a_sigmas <- 2
    b_sigmas <- 0.25 # spatial kernel: sigmasq~IG(shape=a, scale=b), mean=b/(a-1)
    a_lt <- 1
    b_lt <- 1 # temporal kernel: l~Gamma(a,b)
    a_sigmat <- 2
    b_sigmat <- 0.04 # temporal kernel: sigmasq~IG(a, b)

    ####### random noise hyperparameters ######
    a_eps1s <- 5
    b_eps1s <- 0.01 # sigmasq~IG(a,b)
    a_eps2s <- 5
    b_eps2s <- 0.01
    a_eps1t <- 5
    b_eps1t <- 0.01 # sigmasq~IG(a,b)
    a_eps2t <- 5
    b_eps2t <- 0.01

    r <- 1
    Acc <- 0
    y1 <- rep(0, N) # At risk indicator (this is W in paper)
    y1[y > 0] <- 1 # If y>0, then at risk w.p. 1
    q <- rep(.5, N) # 1-p=1/(1+exp(X%*%alpha)), used for updating y1


    #########
    # Inits #
    #########

    #################
    # Fixed Effects #
    #################
    y_ind <- rep(0, N) # convert y to a two class indicator and use logistic regression
    y_ind[y != 0] <- 1

    m1 <- glm(y_ind ~ 0 + X, family = "binomial")
    alpha1 <- m1$coefficients # initial for alpha in the binary component
    m2 <- glm.nb(y[y != 0] ~ 0 + X[y != 0, ])
    beta <- m2$coefficients # initial for beta in the binary component

    eta1 <- X %*% alpha1 + 0
    eta2 <- X %*% beta + 0 # Use all n observations
    p_at_risk <- pmax(0.01, pmin(0.99, exp(eta1) / (1 + exp(eta1)))) # at-risk probability

    q <- 1 / (1 + exp(eta2))
    theta <- p_at_risk * (q^r) / (p_at_risk * (q^r) + 1 - p_at_risk) # Conditional prob that y1=1 given y=0 -- i.e. Pr(chance zero|observed zero)
    y1[y == 0] <- rbinom(sum(y == 0), 1, theta[y == 0]) # If y=0, then draw a "chance zero" w.p. theta, otherwise y1=1

    m1 <- glm(y1 ~ 0 + X, family = "binomial")
    alpha2 <- m1$coefficients # initial for alpha in the binary component
    alpha <- alpha2
    alpha <- as.matrix(alpha)

    ##########################
    # Spatial Random Effects #
    ##########################
    l1s <- l2s <- a_phi / b_phi
    sigma1s <- sigma2s <- sqrt(b_sigmas / (a_sigmas - 1))
    Ks_bin <- sigma1s^2 * exp(-l1s * Ds)
    Ks_nb <- sigma2s^2 * exp(-l2s * Ds)
    a <- t(rmvnorm(n = 1, sigma = Ks_bin))
    c <- t(rmvnorm(n = 1, sigma = Ks_nb))

    #################
    # Temporal Random Effects #
    #################
    sigma1t <- sigma2t <- sqrt(b_sigmat / (a_sigmat - 1))
    l1t <- l2t <- a_lt / b_lt
    Kt_bin <- sigma1t^2 * exp(-Dt / (l1t^2))
    Kt_nb <- sigma2t^2 * exp(-Dt / (l2t^2))
    b <- t(rmvnorm(n = 1, sigma = Kt_bin))
    d <- t(rmvnorm(n = 1, sigma = Kt_nb))

    #################
    # Random Noise #
    #################
    sigma_eps1s <- sqrt(b_eps1s / (a_eps1s - 1))
    sigma_eps2s <- sqrt(b_eps2s / (a_eps2s - 1))
    eps1s <- t(rmvnorm(n = 1, sigma = diag(sigma_eps1s^2, nrow = n)))
    eps2s <- t(rmvnorm(n = 1, sigma = diag(sigma_eps2s^2, nrow = n)))
    sigma_eps1t <- sqrt(b_eps1t / (a_eps1t - 1))
    sigma_eps2t <- sqrt(b_eps2t / (a_eps2t - 1))
    eps1t <- t(rmvnorm(n = 1, sigma = diag(sigma_eps1t^2, nrow = n_time_points)))
    eps2t <- t(rmvnorm(n = 1, sigma = diag(sigma_eps2t^2, nrow = n_time_points)))

    # M-H proposals
    sd_l <- 0.2 # l1t~normal(l1t,sd_l)

    ############
    # Num Sims #
    ############
    lastit <- (nsim - burn) / thin # Last stored value

    #########
    # Store #
    #########
    Beta <- Alpha <- matrix(0, lastit, p)
    R <- R2 <- rep(0, lastit)
    A <- C <- matrix(0, lastit, n)
    B <- D <- matrix(0, lastit, n_time_points)
    Eps1s <- Eps2s <- matrix(0, lastit, n)
    Eps1t <- Eps2t <- matrix(0, lastit, n_time_points)
    L1t <- Sigma1t <- rep(0, lastit)
    L2t <- Sigma2t <- rep(0, lastit)
    Phi_bin <- Sigma1s <- rep(0, lastit)
    Phi_nb <- Sigma2s <- rep(0, lastit)
    Sigma_eps1s <- rep(0, lastit)
    Sigma_eps2s <- rep(0, lastit)
    Sigma_eps1t <- rep(0, lastit)
    Sigma_eps2t <- rep(0, lastit)
    if (save_ypred == TRUE) {
        Y_pred <- matrix(NA, lastit, N)
        y1s <- matrix(NA, lastit, N)
    }


    ########
    # MCMC #
    ########
    NN.matrix <- NNMatrix(
        coords = coords, n.neighbors = M,
        n.omp.threads = 1, search.type = "cb"
    )
    XV <- cbind(X, Vs, Vt)
    for (i in 1:nsim)   {
        Ks_bin <- sigma1s^2 * exp(-l1s * Ds)
        Ks_nb <- sigma2s^2 * exp(-l2s * Ds)
        Kt_bin <- sigma1t^2 * exp(-Dt / (l1t^2))
        Kt_nb <- sigma2t^2 * exp(-Dt / (l2t^2)) # TODO: Recomputation of above variables here seems unnecessary
        Sigma0_bin.inv <- as.matrix(bdiag(
            rcppeigen_invert_matrix(Ks_bin),
            rcppeigen_invert_matrix(Kt_bin)
        ))
        Sigma0_nb.inv <- as.matrix(bdiag(
            rcppeigen_invert_matrix(Ks_nb),
            rcppeigen_invert_matrix(Kt_nb)
        ))
        T0_bin <- as.matrix(bdiag(T0a, Sigma0_bin.inv))
        T0_nb <- as.matrix(bdiag(T0b, Sigma0_nb.inv))

        # Update latent variable z
        mu <- X %*% alpha + Vs %*% a + Vt %*% b + Vs %*% eps1s + Vt %*% eps1t
        w <- rpg(N, 1, mu[, 1])
        z <- (y1 - 1 / 2) / w

        # Update alpha, a, b
        v <- rcppeigen_invert_matrix(crossprod(sqrt(w) * XV) + T0_bin)
        m <- v %*% (t(sqrt(w) * XV) %*% (sqrt(w) * (z - Vs %*% eps1s - Vt %*% eps1t)))
        alphaab <- c(rmvnorm(1, m[, 1], v))
        alpha <- alphaab[1:p]
        a <- alphaab[(p + 1):(p + n)]
        b <- alphaab[-(1:(p + n))]

        # Update eps1s
        v <- solve(crossprod(sqrt(w) * Vs) + 1 / (sigma_eps1s^2) * diag(n))
        m <- v %*% (t(sqrt(w) * Vs) %*% (sqrt(w) * (z - (X %*% alpha + Vs %*% a + Vt %*% b + Vt %*% eps1t))))
        eps1s <- c(rmvnorm(1, m[, 1], v))

        # Update eps1t
        v <- solve(crossprod(sqrt(w) * Vt) + 1 / (sigma_eps1t^2) * diag(n_time_points))
        m <- v %*% (t(sqrt(w) * Vt) %*% (sqrt(w) * (z - (X %*% alpha + Vs %*% a + Vt %*% b + Vs %*% eps1s))))
        eps1t <- c(rmvnorm(1, m[, 1], v))

        # Update at-risk indicator y1 (W in paper)
        eta1 <- X %*% alpha + Vs %*% a + Vt %*% b + Vs %*% eps1s + Vt %*% eps1t
        eta2 <- X %*% beta + Vs %*% c + Vt %*% d + Vs %*% eps2s + Vt %*% eps2t # Use all n observations
        pi <- pmax(0.01, pmin(0.99, exp(eta1) / (1 + exp(eta1)))) # at-risk probability
        q <- 1 / (1 + exp(eta2)) # Pr(y=0|y1=1)
        theta <- pi * (q^r) / (pi * (q^r) + 1 - pi) # Conditional prob that y1=1 given y=0 -- i.e. Pr(chance zero|observed zero)
        y1[y == 0] <- rbinom(sum(y == 0), 1, theta[y == 0]) # If y=0, then draw a "chance zero" w.p. theta, otherwise y1=1
        N1 <- sum(y1)

        # Update r
        rnew <- rtnorm(1, r, sqrt(s), lower = 0) # Treat r as continuous
        ratio <- sum(dnbinom(y[y1 == 1], rnew, q[y1 == 1], log = TRUE)) - sum(dnbinom(y[y1 == 1], r, q[y1 == 1], log = TRUE)) +
            dtnorm(r, rnew, sqrt(s), 0, log = TRUE) - dtnorm(rnew, r, sqrt(s), 0, log = TRUE) # Uniform Prior for R
        # Proposal not symmetric
        if (log(runif(1)) < ratio) {
            r <- rnew
            Acc <- Acc + 1
        }

        # Update latent counts, l
        l <- rep(0, N1)
        ytmp <- y[y1 == 1]
        for (j in 1:N1) l[j] <- sum(rbinom(ytmp[j], 1, round(r / (r + 1:ytmp[j] - 1), 6))) # Could try to avoid loop; rounding avoids numerical instability

        # Update r from conjugate gamma distribution given l and psi
        psi <- exp(eta2[y1 == 1]) / (1 + exp(eta2[y1 == 1]))
        r2 <- rgamma(1, 0.01 + sum(l), 0.01 - sum(log(1 - psi)))

        # update l1t, sigma1t
        l1t_star <- rnorm(1, l1t, sd_l)
        if ((l1t_star < 5) && (l1t_star > 0)) {
            Kt_bin_star <- sigma1t^2 * exp(-Dt / (l1t_star^2))
            likelihood_l1t <- dmvnorm(b, mean = rep(0, n_time_points), sigma = Kt_bin_star, log = TRUE) -
                dmvnorm(b, mean = rep(0, n_time_points), sigma = Kt_bin, log = TRUE)
            prior_l1t <- dgamma(x = l1t_star, shape = a_lt, rate = b_lt, log = TRUE) -
                dgamma(x = l1t, shape = a_lt, rate = b_lt, log = TRUE)
            posterior_l1t <- likelihood_l1t + prior_l1t # prior ratio may get too large and dominate

            if (!is.na(posterior_l1t)) {
                if (log(runif(1)) < posterior_l1t) {
                    l1t <- l1t_star
                    Kt_bin <- Kt_bin_star
                }
            }
        }

        ## posterior dist
        K_inv <- rcppeigen_invert_matrix(Kt_bin + diag(exp(-10), nrow = n_time_points))
        a_new <- a_sigmat + 0.5 * n_time_points
        b_new <- b_sigmat + 0.5 * (t(b) %*% K_inv %*% b)
        sigma1t.sq <- rinvgamma(n = 1, shape = a_new, scale = b_new)
        sigma1t <- sqrt(sigma1t.sq)

        # update phi_bin using M-H
        phi_bin_star <- rnorm(1, l1s, 2) # proposal dist

        if ((phi_bin_star < 16) && (phi_bin_star > 0)) {
            Ks_bin_star <- sigma1s^2 * exp(-phi_bin_star * Ds)
            likelihood_phi_bin <- dmvnorm(a, mean = rep(0, n), sigma = Ks_bin_star, log = TRUE) - dmvnorm(a, mean = rep(0, n), sigma = Ks_bin, log = TRUE)
            prior_phi_bin <- dgamma(x = phi_bin_star, shape = a_phi, rate = b_phi, log = TRUE) - dgamma(x = l1s, shape = a_phi, rate = b_phi, log = TRUE)
            posterior_phi_bin <- likelihood_phi_bin + prior_phi_bin # prior ratio may get too large and dominate

            if (!is.na(posterior_phi_bin)) {
                if (log(runif(1)) < posterior_phi_bin) {
                    l1s <- phi_bin_star
                }
            }
        }

        # update sigma1s
        ## obtain A and D using C and N(i)---cluster
        NN.matrix <- NNMatrix(
            coords = coords, n.neighbors = M,
            n.omp.threads = 1, search.type = "cb"
        ) # replace the input NNmatrix
        AD <- getAD(
            neardist = NN.matrix$NN_dist, neardistM = NN.matrix$NN_distM,
            N = n, M = M, phi = l1s
        )

        Dm <- AD[M + 1, ]
        ind_x <- c(c(rep(2:M, times = 1:(M - 1)), rep(((M + 1):n), each = M)), 1:n)
        ind_y <- c(c(t(NN.matrix$NN_ind))[which(c(t(NN.matrix$NN_ind)) > 0)], 1:n)
        DIA <-
            as.matrix(sparseMatrix(
                i = ind_x, j = (ind_y),
                x = c(-na.omit(as.vector(AD[-(M + 1), ])), rep(1, n))
            ) / sqrt(Dm))

        rho_inv <- t(DIA) %*% DIA
        K_inv <- 1 / (sigma1s^2) * rho_inv

        ## posterior dist
        a_new <- a_sigmas + 0.5 * n
        b_new <- b_sigmas + 0.5 * (t(a) %*% K_inv %*% a)
        sigma1s.sq <- rinvgamma(n = 1, shape = a_new, scale = b_new)
        sigma1s <- sqrt(sigma1s.sq)

        # Update sigma_eps1s
        a_eps1s_new <- a_eps1s + 0.5 * n
        b_eps1s_new <- b_eps1s + 0.5 * sum(eps1s^2)
        sigma_eps1s.sq <- rinvgamma(n = 1, shape = a_eps1s_new, scale = b_eps1s_new)
        sigma_eps1s <- sqrt(sigma_eps1s.sq)

        # Update sigma_eps1t
        a_eps1t_new <- a_eps1t + 0.5 * n_time_points
        b_eps1t_new <- b_eps1t + 0.5 * sum(eps1t^2)
        sigma_eps1t.sq <- rinvgamma(n = 1, shape = a_eps1t_new, scale = b_eps1t_new)
        sigma_eps1t <- sqrt(sigma_eps1t.sq)

        # Update beta
        eta <- X[y1 == 1, ] %*% beta + Vs[y1 == 1, ] %*% c + Vt[y1 == 1, ] %*% d + Vs[y1 == 1, ] %*% eps2s + Vt[y1 == 1, ] %*% eps2t
        w <- rpg(N1, y[y1 == 1] + r, eta) # Polya weights
        z <- (y[y1 == 1] - r) / (2 * w)

        # Update beta, c, d
        v <- rcppeigen_invert_matrix(crossprod(sqrt(w) * XV[y1 == 1, ]) + T0_nb)
        m <- v %*% (t(sqrt(w) * XV[y1 == 1, ]) %*% (sqrt(w) * (z - Vs[y1 == 1, ] %*% eps2s - Vt[y1 == 1, ] %*% eps2t)))
        betacd <- c(rmvnorm(1, m[, 1], v))
        beta <- betacd[1:p]
        c <- betacd[(p + 1):(p + n)]
        d <- betacd[-(1:(p + n))]

        # update eps2s
        v <- solve(crossprod(sqrt(w) * Vs[y1 == 1, ]) + 1 / (sigma_eps2s^2) * diag(n))
        # v<-rcppeigen_invert_matrix(crossprod(sqrt(w)*Vs[y1==1,])+Sigma0e_nb.inv)
        m <- v %*% (t(sqrt(w) * Vs[y1 == 1, ]) %*% (sqrt(w) * (z - X[y1 == 1, ] %*% beta - Vs[y1 == 1, ] %*% c - Vt[y1 == 1, ] %*% d - Vt[y1 == 1, ] %*% eps2t)))
        eps2s <- c(rmvnorm(1, m[, 1], v))

        # update eps2t
        v <- solve(crossprod(sqrt(w) * Vt[y1 == 1, ]) + 1 / (sigma_eps2t^2) * diag(n_time_points))
        m <- v %*% (t(sqrt(w) * Vt[y1 == 1, ]) %*% (sqrt(w) * (z - X[y1 == 1, ] %*% beta - Vs[y1 == 1, ] %*% c - Vt[y1 == 1, ] %*% d - Vs[y1 == 1, ] %*% eps2s)))
        eps2t <- c(rmvnorm(1, m[, 1], v))

        # update l2t, sigma2t
        l2t_star <- rnorm(1, l2t, sd_l)

        if ((l2t_star < 5) && (l2t_star > 0)) {
            Kt_nb_star <- sigma2t^2 * exp(-Dt / (l2t_star^2))
            likelihood_l2t <- dmvnorm(d, mean = rep(0, n_time_points), sigma = Kt_nb_star, log = TRUE) -
                dmvnorm(d, mean = rep(0, n_time_points), sigma = Kt_nb, log = TRUE)
            prior_l2t <- dgamma(x = l2t_star, shape = a_lt, rate = b_lt, log = TRUE) -
                dgamma(x = l2t, shape = a_lt, rate = b_lt, log = TRUE)
            posterior_l2t <- likelihood_l2t + prior_l2t # prior ratio may get too large and dominate

            if (!is.na(posterior_l2t)) {
                if (log(runif(1)) < posterior_l2t) {
                    l2t <- l2t_star
                    Kt_nb <- Kt_nb_star
                }
            }
        }

        ## posterior dist
        K_inv <- rcppeigen_invert_matrix(Kt_nb + diag(exp(-10), nrow = n_time_points))
        a_new <- a_sigmat + 0.5 * n_time_points
        b_new <- b_sigmat + 0.5 * (t(d) %*% K_inv %*% d)
        sigma2t.sq <- rinvgamma(n = 1, shape = a_new, scale = b_new)
        sigma2t <- sqrt(sigma2t.sq)

        # update l2s using M-H
        l2s_star <- rnorm(1, l2s, 1)
        if ((l2s_star < 16) && (l2s_star > 0)) {
            Ks_nb_star <- sigma2s^2 * exp(-l2s_star * Ds)
            likelihood_phi_nb <- dmvnorm(c, mean = rep(0, n), sigma = Ks_nb_star, log = TRUE) -
                dmvnorm(c, mean = rep(0, n), sigma = Ks_nb, log = TRUE)
            prior_phi_nb <- dgamma(x = l2s_star, shape = a_phi, rate = b_phi, log = TRUE) -
                dgamma(x = l2s, shape = a_phi, rate = b_phi, log = TRUE)
            posterior_l2s <- likelihood_phi_nb + prior_phi_nb
            if (!is.na(posterior_l2s)) {
                if (log(runif(1)) < posterior_l2s) {
                    l2s <- l2s_star
                }
            }
        }

        # update sigma2s
        ## obtain A and D using C and N(i)---cluster
        NN.matrix <- NNMatrix(
            coords = coords, n.neighbors = M,
            n.omp.threads = 1, search.type = "cb"
        ) # replace the input NNmatrix
        AD <- getAD(
            neardist = NN.matrix$NN_dist, neardistM = NN.matrix$NN_distM,
            N = n, M = M, phi = l2s
        )
        Dm <- AD[M + 1, ]
        ind_x <- c(c(rep(2:M, times = 1:(M - 1)), rep(((M + 1):n), each = M)), 1:n)
        ind_y <- c(c(t(NN.matrix$NN_ind))[which(c(t(NN.matrix$NN_ind)) > 0)], 1:n)
        DIA <-
            as.matrix(sparseMatrix(
                i = ind_x, j = (ind_y),
                x = c(-na.omit(as.vector(AD[-(M + 1), ])), rep(1, n))
            ) / sqrt(Dm))

        rho_inv <- t(DIA) %*% DIA
        K_inv <- 1 / (sigma2s^2) * rho_inv

        ## posterior dist
        a_new <- a_sigmas + 0.5 * n
        b_new <- b_sigmas + 0.5 * (t(c) %*% K_inv %*% c)
        sigma2s.sq <- rinvgamma(n = 1, shape = a_new, scale = b_new)
        sigma2s <- sqrt(sigma2s.sq)

        # Update sigma_eps2s
        a_eps2s_new <- a_eps2s + 0.5 * n
        b_eps2s_new <- b_eps2s + 0.5 * sum(eps2s^2)
        sigma_eps2s.sq <- rinvgamma(n = 1, shape = a_eps2s_new, scale = b_eps2s_new)
        sigma_eps2s <- sqrt(sigma_eps2s.sq)

        # Update sigma_eps2t
        a_eps2t_new <- a_eps2t + 0.5 * n_time_points
        b_eps2t_new <- b_eps2t + 0.5 * sum(eps2t^2)
        sigma_eps2t.sq <- rinvgamma(n = 1, shape = a_eps2t_new, scale = b_eps2t_new)
        sigma_eps2t <- sqrt(sigma_eps2t.sq)

        if (save_ypred) {
            y_pred <- estimate(X, alpha, beta, Vs, Vt, a, b, c, d, eps1s, eps1t, eps2s, eps2t, r) # ?
        }
        # Store
        if ((i > burn) && (i %% thin == 0)) {
            j <- (i - burn) / thin
            Alpha[j, ] <- alpha
            Beta[j, ] <- beta # fixed effects
            A[j, ] <- a
            B[j, ] <- b
            C[j, ] <- c
            D[j, ] <- d # random effects
            Eps1s[j, ] <- eps1s
            Eps2s[j, ] <- eps2s # random noise
            Eps1t[j, ] <- eps1t
            Eps2t[j, ] <- eps2t # random noise
            L1t[j] <- l1t
            Sigma1t[j] <- sigma1t # temporal hyperparameters
            L2t[j] <- l2t
            Sigma2t[j] <- sigma2t # temporal hyperparameters
            Phi_bin[j] <- l1s
            Sigma1s[j] <- sigma1s # spatial hyperparameters
            Phi_nb[j] <- l2s
            Sigma2s[j] <- sigma2s # spatial hyperparameters
            Sigma_eps1s[j] <- sigma_eps1s
            Sigma_eps2s[j] <- sigma_eps2s # random noise parameters
            Sigma_eps1t[j] <- sigma_eps1t
            Sigma_eps2t[j] <- sigma_eps2t # random noise parameters
            R[j] <- r
            R2[j] <- r2
            if (save_ypred == TRUE) {
                Y_pred[j, ] <- y_pred
                y1s[j, ] <- y1
            }
        }
        if ((i %% print_iter == 0) && (print_progress)) print(i)
    }
    # Put the results into a list
    results <- list(
        Alpha = Alpha, Beta = Beta, A = A, B = B, C = C, D = D,
        Eps1s = Eps1s, Eps2s = Eps2s, Eps1t = Eps1t, Eps2t = Eps2t,
        L1t = L1t, Sigma1t = Sigma1t, L2t = L2t, Sigma2t = Sigma2t,
        Phi_bin = Phi_bin, Sigma1s = Sigma1s, Phi_nb = Phi_nb, Sigma2s = Sigma2s,
        Sigma_eps1s = Sigma_eps1s, Sigma_eps2s = Sigma_eps2s,
        Sigma_eps1t = Sigma_eps1t, Sigma_eps2t = Sigma_eps2t,
        R = R, R2 = R2
    )
    if (save_ypred) {
        temp <- list(Y_pred = Y_pred, at_risk = y1s)
        results <- append(results, temp)
    }
    return(results)
}
