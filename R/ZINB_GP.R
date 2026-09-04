#' Mix a Matrix with Identity Noise
#'
#' @description Creates `noise_ratio * A + (1 - noise_ratio) * I`.
#'
#' @param A Numeric square matrix.
#' @param noise_ratio Numeric mixing ratio between zero and one.
#'
#' @return A numeric square matrix with the same dimensions as `A`.
#' @examples
#' correlation <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#' noise_mix(correlation, noise_ratio = 0.8)
#' @export
noise_mix <- function(A, noise_ratio) {
    return(noise_ratio * A + diag(1 - noise_ratio, nrow=nrow(A)))
}

#' Squared Exponential Kernel
#'
#' @description Creates the kernel matrix `exp(-dist / ls^2)`.
#'
#' @param dist Numeric distance matrix.
#' @param ls Positive length scale.
#'
#' @return A numeric kernel matrix with the same dimensions as `dist`.
#' @examples
#' distances <- as.matrix(stats::dist(c(0, 1, 3)))
#' kernel(distances, ls = 2)
#' @export
kernel <- function(dist, ls) {
    return(exp(-dist / (ls^2)))
}

#' nullcheck
#' @description Returns default if value is null
#' 
#' @param value Nullable
#' @param default Default value
nullcheck <- function(value, default) {
    if (is.null(value)) {
        return(default)
    }
    return(value)
}

#' Validate ZINB Model Inputs
#'
#' Checks the fixed-effect design, count response, and MCMC controls shared by
#' the branch-specific ZINB-GP samplers.
#'
#' @param X Fixed-effect design matrix with one row per observation.
#' @param y Non-negative integer count response with one value per row of `X`.
#' @param nsim Positive integer giving the total number of MCMC iterations.
#' @param burn Non-negative integer giving the number of burn-in iterations.
#' @param thin Positive integer giving the post-burn-in thinning interval.
#'
#' @return `TRUE`, invisibly, when all inputs are valid. Otherwise, the function
#'   stops with an input-specific error message.
#' @keywords internal
validate_zinb_inputs <- function(X, y, nsim, burn, thin) {
    if (length(dim(X)) != 2L || nrow(X) < 1L || ncol(X) < 1L) {
        stop("X must be a non-empty two-dimensional design matrix.")
    }
    if (length(y) != nrow(X)) {
        stop("y must have the same number of entries as rows in X.")
    }
    if (any(!is.finite(X))) {
        stop("X must contain only finite values.")
    }
    if (any(!is.finite(y)) || any(y < 0) || any(y != floor(y))) {
        stop("y must contain only finite, non-negative integer counts.")
    }
    if (!any(y > 0)) {
        stop("At least one positive count is required to initialize glm.nb.")
    }

    is_whole_number <- function(value) {
        length(value) == 1L &&
            is.finite(value) &&
            value == floor(value)
    }
    if (!is_whole_number(nsim) || nsim < 1L) {
        stop("nsim must be a positive integer.")
    }
    if (!is_whole_number(burn) || burn < 0L) {
        stop("burn must be a non-negative integer.")
    }
    if (!is_whole_number(thin) || thin < 1L) {
        stop("thin must be a positive integer.")
    }
    if (nsim <= burn) {
        stop("nsim must be greater than burn.")
    }
    if ((nsim - burn) %% thin != 0) {
        stop("(nsim - burn) must be divisible by thin.")
    }

    invisible(TRUE)
}

#' Validate a GP Design and Distance Matrix
#'
#' Checks that a baseline-coded GP design is aligned with the observations and
#' that its distance matrix contains one additional row and column for the
#' omitted baseline level.
#'
#' @param V GP design matrix with one row per observation and one column per
#'   nonbaseline GP level.
#' @param D Pairwise distance matrix for all GP levels, including the baseline.
#' @param n_observations Expected number of rows in `V`.
#' @param design_name Character label used to identify `V` in error messages.
#' @param distance_name Character label used to identify `D` in error messages.
#'
#' @return `TRUE`, invisibly, when the design and distance matrix are valid.
#'   Otherwise, the function stops with an input-specific error message.
#' @keywords internal
validate_gp_design <- function(V, D, n_observations, design_name, distance_name) {
    if (length(dim(V)) != 2L || nrow(V) != n_observations) {
        stop(design_name, " must have the same number of rows as X.")
    }
    if (ncol(V) < 1L) {
        stop(design_name, " must contain at least one nonbaseline GP level.")
    }
    if (any(!is.finite(V))) {
        stop(design_name, " must contain only finite values.")
    }

    D <- as.matrix(D)
    expected_dimension <- ncol(V) + 1L
    if (length(dim(D)) != 2L ||
        nrow(D) != expected_dimension ||
        ncol(D) != expected_dimension) {
        stop(
            distance_name,
            " must be a square matrix with ncol(",
            design_name,
            ") + 1 rows and columns, including the baseline level."
        )
    }
    if (any(!is.finite(D)) || any(D < 0)) {
        stop(distance_name, " must contain finite, non-negative distances.")
    }
    if (any(abs(diag(D)) > sqrt(.Machine$double.eps))) {
        stop("The diagonal of ", distance_name, " must be zero.")
    }
    if (!isTRUE(all.equal(D, t(D), tolerance = sqrt(.Machine$double.eps)))) {
        stop(distance_name, " must be symmetric.")
    }

    invisible(TRUE)
}

#' Validate Spatial Coordinate Alignment
#'
#' Confirms that the full coordinate matrix, the baseline-coded spatial design,
#' and the spatial distance matrix describe the same number of spatial levels.
#'
#' @param coords Spatial coordinate matrix with one row per full spatial level,
#'   including the baseline level.
#' @param Vs Spatial GP design matrix with the baseline column omitted.
#' @param Ds Spatial distance matrix with one row and column per row of `coords`.
#'
#' @return `TRUE`, invisibly, when the spatial inputs are aligned. Otherwise,
#'   the function stops with an input-specific error message.
#' @keywords internal
validate_spatial_coordinates <- function(coords, Vs, Ds) {
    if (length(dim(coords)) != 2L || any(!is.finite(coords))) {
        stop("coords must be a finite two-dimensional coordinate matrix.")
    }

    expected_levels <- ncol(Vs) + 1L
    if (nrow(coords) != expected_levels) {
        stop(
            "nrow(coords) must equal ncol(Vs) + 1. ",
            "Vs must be generated from the supplied coordinates with its ",
            "baseline coordinate column omitted."
        )
    }
    if (nrow(Ds) != nrow(coords) || ncol(Ds) != nrow(coords)) {
        stop(
            "Ds must have one row and column for every row of coords, ",
            "including the baseline coordinate."
        )
    }

    invisible(TRUE)
}

#' Conditional At-Risk Probability for an Observed Zero
#'
#' Computes the conditional probability that an observation belongs to the
#' negative-binomial, or at-risk, component given that its observed count is
#' zero. The calculation is performed on the log-odds scale for numerical
#' stability.
#'
#' @param pi Prior at-risk probabilities.
#' @param q Negative-binomial zero-probability bases, equal to
#'   `1 - sigmoid(eta_count)`.
#' @param r Positive negative-binomial dispersion parameter.
#'
#' @return A numeric vector or matrix with the same shape as `pi` and `q`,
#'   containing probabilities between zero and one.
#' @keywords internal
zero_at_risk_probability <- function(pi, q, r) {
    stats::plogis(stats::qlogis(pi) + r * log(q))
}

#' update_ls_sigma_noise
#' @description Update kernel parameters for a GP
#' 
#' @param ls Current length scale
#' @param sigma Current sigma
#' @param noise_ratio Current noise ratio
#' @param gpdraw Last draw from the gp with these parameters
#' @param K Current kernel matrix
#' @param D Distance matrix
#' @param lsPrior prior information for length scale, needs mh_sd, max, a, b
#' @param sigmaPrior prior information for sigma, needs a, b
#' @param noisePrior prior information for noise_ratio, needs mh_sd, a, b
#' @param kern Kernel function accepting a distance matrix and a length scale.
#' @return A list with the following elements:
#' \describe{
#'   \item{ls}{Updated length scale.}
#'   \item{sigma}{Updated GP scale.}
#'   \item{noise_ratio}{Updated kernel-to-noise mixing ratio.}
#'   \item{K}{Updated kernel matrix.}
#'   \item{K_inv}{Inverse of the updated kernel matrix.}
#' }
#' @importFrom stats dgamma
#' @importFrom stats dbeta
#' @importFrom stats runif
#' @importFrom mvtnorm dmvnorm
#' @importFrom msm rtnorm
#' @importFrom msm dtnorm
update_ls_sigma_noise <- function(ls, sigma, noise_ratio, gpdraw, K, D, lsPrior, sigmaPrior, noisePrior, kern) {
    # update ls
    # Consider using exponential proposals instead
    proposal <- rtnorm(1, mean = ls, sd = lsPrior$mh_sd, lower = 1e-6, upper = lsPrior$max) #max(min(rnorm(1, ls, lsPrior$mh_sd), lsPrior$max - 1), 1e-6)
    if (TRUE) {
        K_star <- sigma^2 * noise_mix(kern(D, proposal), noise_ratio)
        
        # Calculate model likelihood
        likelihood_ls <- dmvnorm(gpdraw, mean = rep(0, length(gpdraw)), sigma = K_star, log = TRUE) -
            dmvnorm(gpdraw, mean = rep(0, length(gpdraw)), sigma = K, log = TRUE)
        
        # Calculate prior likelihood
        prior_ls <- dgamma(x = proposal, shape = lsPrior$a, rate = lsPrior$b, log = TRUE) -
            dgamma(x = ls, shape = lsPrior$a, rate = lsPrior$b, log = TRUE)
        
        # Calculate transition probabilities
        trans_ls <- dtnorm(x = ls, mean = proposal, sd = lsPrior$mh_sd, lower = 1e-6, upper = lsPrior$max, log = 1) - 
            dtnorm(x = proposal, mean = ls, sd = lsPrior$mh_sd, lower = 1e-6, upper = lsPrior$max, log = 1)

        posterior_ls <- likelihood_ls + prior_ls + trans_ls

        if (!is.na(posterior_ls)) {
            if (log(runif(1)) < posterior_ls) {
                ls <- proposal
                K <- K_star
            }
        }
    }
    
    # Update noise ratio
    eps_nr <- 2.06115369216775e-09
    proposal <- rtnorm(1, mean = noise_ratio, sd = noisePrior$mh_sd, lower = eps_nr, upper = 1 - eps_nr) #max(min(rnorm(1, noise_ratio, noisePrior$mh_sd), 1 - 1e-7), 0 + 1e-7)
    if (TRUE) {
        K_star <- sigma^2 * noise_mix(kern(D, ls), proposal)
        
        # Calculate model likelihood
        likelihood_nr <- dmvnorm(gpdraw, mean = rep(0, length(gpdraw)), sigma = K_star, log = TRUE) -
            dmvnorm(gpdraw, mean = rep(0, length(gpdraw)), sigma = K, log = TRUE)
        
        # Calculate prior likelihood
        prior_nr <- dbeta(x = proposal, shape1 = noisePrior$a, shape2 = noisePrior$b, log = TRUE) -
            dbeta(x = noise_ratio, shape1 = noisePrior$a, shape2 = noisePrior$b, log = TRUE)

        # Calculate transition probabilities
        trans_ls <- dtnorm(x = noise_ratio, mean = proposal, sd = noisePrior$mh_sd, lower = eps_nr, upper = 1 - eps_nr, log = 1) - 
            dtnorm(x = proposal, mean = noise_ratio, sd = noisePrior$mh_sd, lower = eps_nr, upper = 1 - eps_nr, log = 1)

        posterior_nr <- likelihood_nr + prior_nr + trans_ls

        if (!is.na(posterior_nr)) {
            if (log(runif(1)) < posterior_nr) {
                noise_ratio <- proposal
            }
        }
    }

    # update sigma1t, kernel matrices
    K_nosigma <- noise_mix(kern(D, ls), noise_ratio)
    K_nosigma_inv <- forceSymmetric(solve(K_nosigma))
    a_new <- sigmaPrior$a + 0.5 * length(gpdraw)
    b_new <- sigmaPrior$b + 0.5 * (t(gpdraw) %*% K_nosigma_inv %*% gpdraw)
    sigma.sq <- rinvgamma(n = 1, shape = a_new, scale = b_new[1,1])
    sigma <- sqrt(sigma.sq)

    K <- K_nosigma * sigma.sq
    K_inv <- K_nosigma_inv * (1 / sigma.sq)

    return(list(ls=ls, sigma=sigma, noise_ratio=noise_ratio, K=K, K_inv=K_inv))
}

#' ZINB_GP_orig
#' @description Fits the spatial-temporal zero-inflated negative-binomial
#' Gaussian-process model described in \doi{10.1016/j.jspi.2023.106098}.
#'
#' @param X Fixed-effect design matrix with one row per observation.
#' @param y Non-negative integer count response.
#' @param coords Spatial coordinate matrix with one row per full spatial level,
#'   including the baseline level omitted from `Vs`.
#' @param Vs Spatial random-effect design matrix with one row per observation.
#' @param Vt Temporal random-effect design matrix with one row per observation.
#' @param Ds Spatial distance matrix; its diagonal must be zero.
#' @param Dt Temporal distance matrix; its diagonal must be zero.
#' @param nsim Total number of MCMC iterations; must exceed `burn`.
#' @param burn Number of burn-in iterations.
#' @param thin Store every `thin`-th iteration after burn-in.
#' @param save_ypred Whether to save posterior predictive draws.
#' @param print_iter Print progress every `print_iter` iterations.
#' @param print_progress Whether to print MCMC progress.
#' @param ltPrior List with `max`, `mh_sd`, `a`, and `b` for temporal length-scale prior and proposal controls.
#' @param lsPrior List with `max`, `mh_sd`, `a`, and `b` for spatial length-scale prior and proposal controls.
#' @param sigmaPrior List with `a` and `b` inverse-gamma prior parameters for GP scales.
#' @param noisePrior List with `a`, `b`, and `mh_sd` for the GP noise-ratio prior and proposal.
#' @param mh_sd_r Proposal standard deviation for the negative-binomial dispersion parameter.
#' @param kern Kernel function accepting a distance matrix and length scale.
#' @return A list containing posterior MCMC draws, including fixed-effect
#'   coefficients, spatial and temporal random effects, their GP parameters,
#'   the negative-binomial dispersion parameter, and optional predictive draws.
#' @importFrom MASS glm.nb
#' @importFrom mvtnorm rmvnorm
#' @importFrom BayesLogit rpg
#' @importFrom Matrix bdiag
#' @importFrom Matrix sparseMatrix
#' @importFrom msm rtnorm
#' @importFrom msm dtnorm
#' @importFrom stats rbinom
#' @importFrom stats dnbinom
#' @importFrom stats rbinom
#' @importFrom stats glm
#' @importFrom stats rnorm
#' @importFrom LaplacesDemon rinvgamma
#' @importFrom stats runif
#' @importFrom Matrix forceSymmetric
ZINB_GP_orig <- function(X, y, coords, Vs, Vt, Ds, Dt, nsim, burn, thin = 1, save_ypred = FALSE, print_iter = 100, print_progress = FALSE, ltPrior = NULL, lsPrior = NULL, sigmaPrior = NULL, noisePrior = NULL, mh_sd_r = NULL, kern = NULL) {
    validate_zinb_inputs(X, y, nsim, burn, thin)
    validate_gp_design(Vs, Ds, nrow(X), "Vs", "Ds")
    validate_gp_design(Vt, Dt, nrow(X), "Vt", "Dt")
    validate_spatial_coordinates(coords, Vs, Ds)

    # X is the design matrix with dimension N*p
    # x is the vector with length N
    # y is the count response with length N
    n <- nrow(coords) - 1 # number of clusters
    N <- nrow(X) # number of observations
    p <- ncol(X) # dimension of alpha and beta
    n_time_points <- ncol(Vt)

    # Sacrifice to the intercept gods
    Ds <- Ds[2:nrow(Ds), 2:ncol(Ds)]
    Dt <- Dt[2:nrow(Dt), 2:ncol(Dt)]

    # Use squared distances
    Ds <- Ds * Ds
    Dt <- Dt * Dt

    # Find reasonable bounds for GP length scales
    kern <- nullcheck(kern, kernel)
    param_bounds <- gp_param_bounds(Ds, Dt, kern)
    lsmax <- param_bounds$lsmax
    ltmax <- param_bounds$ltmax

    ##########
    # Priors #
    ##########

    ####### priors for alpha and beta ######
    T0a <- T0b <- diag(0.01, p)
    sd_r <- nullcheck(mh_sd_r, 0.4)
    
    ####### kernel hyperparameters  ######
    ltPrior <- nullcheck(ltPrior, list(max=ltmax, mh_sd=3, a=1, b=0.001))
    lsPrior <- nullcheck(lsPrior, list(max=lsmax, mh_sd=3, a=1, b=0.001))
    sigmaPrior <- nullcheck(sigmaPrior, list(a=0.01, b=0.1))
    noisePrior <- nullcheck(noisePrior, list(a=1.5, b=1.5, mh_sd=0.2))

    # Model init
    r <- 1
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
    p_at_risk <- sigmoid(eta1) # at-risk probability

    q <- sigmoid(-eta2)
    theta <- zero_at_risk_probability(p_at_risk, q, r)
    y1[y == 0] <- rbinom(sum(y == 0), 1, theta[y == 0]) # If y=0, then draw a "chance zero" w.p. theta, otherwise y1=1

    m1 <- glm(y1 ~ 0 + X, family = "binomial")
    alpha <- as.matrix(m1$coefficients) # initial for alpha in the binary component

    noise_ratio_t1 <- 0.5
    noise_ratio_t2 <- 0.5
    noise_ratio_s1 <- 0.5
    noise_ratio_s2 <- 0.5

    ##########################
    # Spatial Random Effects #
    ##########################
    l1s <- l2s <- 1
    sigma1s <- sigma2s <- 2
    Ks_bin <- sigma1s^2 * noise_mix(kern(Ds, l1s), noise_ratio_s1)
    Ks_bin_inv <- forceSymmetric(solve(Ks_bin))
    Ks_nb <- sigma2s^2 * noise_mix(kern(Ds, l2s), noise_ratio_s2)
    Ks_nb_inv <- forceSymmetric(solve(Ks_nb))
    a <- t(rmvnorm(n = 1, sigma = Ks_bin))
    c <- t(rmvnorm(n = 1, sigma = Ks_nb))

    #################
    # Temporal Random Effects #
    #################
    sigma1t <- sigma2t <- 2
    l1t <- l2t <- 1
    Kt_bin <- sigma1t^2 * noise_mix(kern(Dt, l1t), noise_ratio_t1)
    Kt_bin_inv <- forceSymmetric(solve(Kt_bin))
    Kt_nb <- sigma2t^2 * noise_mix(kern(Dt, l2t), noise_ratio_t2)
    Kt_nb_inv <- forceSymmetric(solve(Kt_nb))
    b <- t(rmvnorm(n = 1, sigma = Kt_bin))
    d <- t(rmvnorm(n = 1, sigma = Kt_nb))

    ############
    # Num Sims #
    ############
    lastit <- (nsim - burn) / thin # Last stored value

    #########
    # Store #
    #########
    Beta <- Alpha <- matrix(0, lastit, p)
    R <- rep(0, lastit)
    A <- C <- matrix(0, lastit, n)
    B <- D <- matrix(0, lastit, n_time_points)
    L1t <- Sigma1t <- Noise1t <- rep(0, lastit)
    L2t <- Sigma2t <- Noise2t <- rep(0, lastit)
    L1s <- Sigma1s <- Noise1s <- rep(0, lastit)
    L2s <- Sigma2s <- Noise2s <- rep(0, lastit)
    if (save_ypred == TRUE) {
        Y_pred <- matrix(NA, lastit, N)
        y1s <- matrix(NA, lastit, N)
    }


    ########
    # MCMC #
    ########
    XV <- cbind(X, Vs, Vt)
    for (i in 1:nsim)   {
        # Ensure these are all updated from last iteration
        Sigma0_bin.inv <- as.matrix(bdiag(
            Ks_bin_inv,
            Kt_bin_inv
        ))
        Sigma0_nb.inv <- as.matrix(bdiag(
            Ks_nb_inv,
            Kt_nb_inv
        ))
        T0_bin <- as.matrix(bdiag(T0a, Sigma0_bin.inv))
        T0_nb <- as.matrix(bdiag(T0b, Sigma0_nb.inv))

        # Update latent variable z
        mu <- X %*% alpha + Vs %*% a + Vt %*% b
        w <- rpg(N, 1, mu[, 1])
        z <- (y1 - 1 / 2) / w

        # Update alpha, a, b
        svd_vinv <- svd(crossprod(sqrt(w) * XV) + T0_bin)
        m <- solve_svd(svd_vinv, (t(sqrt(w) * XV) %*% (sqrt(w) * z)))
        alphaab <- c(mvn_sample_svd(svd_vinv, m))
        alpha <- alphaab[1:p]
        a <- alphaab[(p + 1):(p + n)]
        b <- alphaab[-(1:(p + n))]

        # Update at-risk indicator y1 (W in paper)
        eta1 <- X %*% alpha + Vs %*% a + Vt %*% b
        eta2 <- X %*% beta + Vs %*% c + Vt %*% d # Use all n observations
        pi <- sigmoid(eta1) # at-risk probability
        q <- sigmoid(-eta2) # Pr(y=0|y1=1)
        theta <- zero_at_risk_probability(pi, q, r)
        y1[y == 0] <- rbinom(sum(y == 0), 1, theta[y == 0]) # If y=0, then draw a "chance zero" w.p. theta, otherwise y1=1
        N1 <- sum(y1)

        # Update r, TODO: Check this.
        rnew <- rtnorm(1, r, sd_r, lower = 0) # Treat r as continuous
        ratio <- sum(dnbinom(y[y1 == 1], rnew, q[y1 == 1], log = TRUE)) - sum(dnbinom(y[y1 == 1], r, q[y1 == 1], log = TRUE)) +
            dtnorm(r, rnew, sd_r, 0, log = TRUE) - dtnorm(rnew, r, sd_r, 0, log = TRUE) # Uniform Prior for R
        # Proposal not symmetric
        if (log(runif(1)) < ratio) {
            r <- rnew
        }

        # update l1t, sigma1t, noise_ratio_t1
        out <- update_ls_sigma_noise(l1t, sigma1t, noise_ratio_t1, b, Kt_bin, Dt, ltPrior, sigmaPrior, noisePrior, kern)
        l1t <- out$ls
        sigma1t <- out$sigma
        noise_ratio_t1 <- out$noise_ratio
        Kt_bin <- out$K
        Kt_bin_inv <- out$K_inv

        # update l1s, sigma1s, noise_ratio_s1
        out <- update_ls_sigma_noise(l1s, sigma1s, noise_ratio_s1, a, Ks_bin, Ds, lsPrior, sigmaPrior, noisePrior, kern)
        l1s <- out$ls
        sigma1s <- out$sigma
        noise_ratio_s1 <- out$noise_ratio
        Ks_bin <- out$K
        Ks_bin_inv <- out$K_inv

        # Update beta
        eta <- X[y1 == 1, ] %*% beta + Vs[y1 == 1, ] %*% c + Vt[y1 == 1, ] %*% d
        w <- rpg(N1, y[y1 == 1] + r, eta) # Polya weights
        z <- (y[y1 == 1] - r) / (2 * w)

        # Update beta, c, d
        svd_vinv <- svd(crossprod(sqrt(w) * XV[y1 == 1, ]) + T0_nb)
        m <- solve_svd(svd_vinv, (t(sqrt(w) * XV[y1 == 1, ]) %*% (sqrt(w) * z)))
        betacd <- c(mvn_sample_svd(svd_vinv, m))
        beta <- betacd[1:p]
        c <- betacd[(p + 1):(p + n)]
        d <- betacd[-(1:(p + n))]

        # update l2t, sigma2t, noise_ratio_t2
        out <- update_ls_sigma_noise(l2t, sigma2t, noise_ratio_t2, d, Kt_nb, Dt, ltPrior, sigmaPrior, noisePrior, kern)
        l2t <- out$ls
        sigma2t <- out$sigma
        noise_ratio_t2 <- out$noise_ratio
        Kt_nb <- out$K
        Kt_nb_inv <- out$K_inv

        # update l2s, sigma2s, noise_ratio_s2
        out <- update_ls_sigma_noise(l2s, sigma2s, noise_ratio_s2, c, Ks_nb, Ds, lsPrior, sigmaPrior, noisePrior, kern)
        l2s <- out$ls
        sigma2s <- out$sigma
        noise_ratio_s2 <- out$noise_ratio
        Ks_nb <- out$K
        Ks_nb_inv <- out$K_inv

        # Store
        if ((i > burn) && ((i - burn) %% thin == 0)) {
            j <- (i - burn) / thin
            
            Alpha[j, ] <- alpha
            Beta[j, ] <- beta # fixed effects

            A[j, ] <- a
            B[j, ] <- b
            C[j, ] <- c
            D[j, ] <- d # random effects
            
            L1t[j] <- l1t
            Noise1t[j] <- noise_ratio_t1
            Sigma1t[j] <- sigma1t # temporal hyperparameters
            
            L2t[j] <- l2t
            Noise2t[j] <- noise_ratio_t2
            Sigma2t[j] <- sigma2t # temporal hyperparameters
            
            L1s[j] <- l1s
            Noise1s[j] <- noise_ratio_s1
            Sigma1s[j] <- sigma1s # spatial hyperparameters
            
            L2s[j] <- l2s
            Noise2s[j] <- noise_ratio_s2
            Sigma2s[j] <- sigma2s # spatial hyperparameters
            
            R[j] <- r
            if (save_ypred) {
                eta_at_risk <- X %*% alpha + Vs %*% a + Vt %*% b
                eta_count <- X %*% beta + Vs %*% c + Vt %*% d
                Y_pred[j, ] <- draw_zinb_response(
                    eta_at_risk = eta_at_risk,
                    eta_count = eta_count,
                    r = r
                )
                y1s[j, ] <- y1
            }
        }
        if ((i %% print_iter == 0) && (print_progress)) message("Iteration ", i, " of ", nsim)
    }
    # Put the results into a list
    results <- list(
        Alpha = Alpha, Beta = Beta, A = A, B = B, C = C, D = D,
        L1t = L1t, Sigma1t = Sigma1t, Noise1t=Noise1t, L2t = L2t, Sigma2t = Sigma2t, Noise2t=Noise2t,
        L1s = L1s, Sigma1s = Sigma1s, Noise1s=Noise1s, L2s = L2s, Sigma2s = Sigma2s, Noise2s=Noise2s,
        R = R
    )
    if (save_ypred) {
        temp <- list(Y_pred = Y_pred, at_risk = y1s)
        results <- append(results, temp)
    }
    return(new_zinb_gp_fit(results))
}


#' ZINB_GP_inflation
#' @description Run the ZINB GP model with GP random effects only in the inflation component
#' @inheritParams ZINB_GP_orig
#' @return A list containing posterior MCMC draws for the zero-inflation GP,
#'   fixed effects, and the negative-binomial dispersion parameter.
#' @importFrom MASS glm.nb
#' @importFrom mvtnorm rmvnorm
#' @importFrom BayesLogit rpg
#' @importFrom Matrix bdiag
#' @importFrom Matrix sparseMatrix
#' @importFrom msm rtnorm
#' @importFrom msm dtnorm
#' @importFrom stats rbinom
#' @importFrom stats dnbinom
#' @importFrom stats rbinom
#' @importFrom stats glm
#' @importFrom stats rnorm
#' @importFrom LaplacesDemon rinvgamma
#' @importFrom stats runif
#' @importFrom Matrix forceSymmetric
ZINB_GP_inflation <- function(X, y, coords, Vs, Vt, Ds, Dt, nsim, burn, thin = 1, save_ypred = FALSE, print_iter = 100, print_progress = FALSE, ltPrior = NULL, lsPrior = NULL, sigmaPrior = NULL, noisePrior = NULL, mh_sd_r = NULL, kern = NULL) {
    validate_zinb_inputs(X, y, nsim, burn, thin)
    validate_gp_design(Vs, Ds, nrow(X), "Vs", "Ds")
    validate_gp_design(Vt, Dt, nrow(X), "Vt", "Dt")
    validate_spatial_coordinates(coords, Vs, Ds)

    # X is the design matrix with dimension N*p
    # x is the vector with length N
    # y is the count response with length N
    n <- nrow(coords) - 1 # number of clusters
    N <- nrow(X) # number of observations
    p <- ncol(X) # dimension of alpha and beta
    n_time_points <- ncol(Vt)

    # Sacrifice to the intercept gods
    Ds <- Ds[2:nrow(Ds), 2:ncol(Ds)]
    Dt <- Dt[2:nrow(Dt), 2:ncol(Dt)]

    # Use squared distances
    Ds <- Ds * Ds
    Dt <- Dt * Dt

    # Find reasonable bounds for GP length scales
    kern <- nullcheck(kern, kernel)
    param_bounds <- gp_param_bounds(Ds, Dt, kern)
    lsmax <- param_bounds$lsmax
    ltmax <- param_bounds$ltmax

    ##########
    # Priors #
    ##########

    ####### priors for alpha and beta ######
    T0a <- T0b <- diag(0.01, p)
    sd_r <- nullcheck(mh_sd_r, 0.4)
    
    ####### kernel hyperparameters  ######
    ltPrior <- nullcheck(ltPrior, list(max=ltmax, mh_sd=3, a=1, b=0.001))
    lsPrior <- nullcheck(lsPrior, list(max=lsmax, mh_sd=3, a=1, b=0.001))
    sigmaPrior <- nullcheck(sigmaPrior, list(a=0.01, b=0.1))
    noisePrior <- nullcheck(noisePrior, list(a=1.5, b=1.5, mh_sd=0.2))

    # Model init
    r <- 1
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
    p_at_risk <- sigmoid(eta1) # at-risk probability

    q <- sigmoid(-eta2)
    theta <- zero_at_risk_probability(p_at_risk, q, r)
    y1[y == 0] <- rbinom(sum(y == 0), 1, theta[y == 0]) # If y=0, then draw a "chance zero" w.p. theta, otherwise y1=1

    m1 <- glm(y1 ~ 0 + X, family = "binomial")
    alpha <- as.matrix(m1$coefficients) # initial for alpha in the binary component

    noise_ratio_t1 <- 0.5
    noise_ratio_t2 <- 0.5
    noise_ratio_s1 <- 0.5
    noise_ratio_s2 <- 0.5

    ##########################
    # Spatial Random Effects #
    ##########################
    l1s <- 1
    sigma1s <- 2
    Ks_bin <- sigma1s^2 * noise_mix(kern(Ds, l1s), noise_ratio_s1)
    Ks_bin_inv <- forceSymmetric(solve(Ks_bin))
    a <- t(rmvnorm(n = 1, sigma = Ks_bin))

    #################
    # Temporal Random Effects #
    #################
    sigma1t <- 2
    l1t <- 1
    Kt_bin <- sigma1t^2 * noise_mix(kern(Dt, l1t), noise_ratio_t1)
    Kt_bin_inv <- forceSymmetric(solve(Kt_bin))
    b <- t(rmvnorm(n = 1, sigma = Kt_bin))

    ############
    # Num Sims #
    ############
    lastit <- (nsim - burn) / thin # Last stored value

    #########
    # Store #
    #########
    Beta <- Alpha <- matrix(0, lastit, p)
    R <- rep(0, lastit)
    A <- matrix(0, lastit, n)
    B <- matrix(0, lastit, n_time_points)
    L1t <- Sigma1t <- Noise1t <- rep(0, lastit)
    L1s <- Sigma1s <- Noise1s <- rep(0, lastit)
    if (save_ypred == TRUE) {
        Y_pred <- matrix(NA, lastit, N)
        y1s <- matrix(NA, lastit, N)
    }


    ########
    # MCMC #
    ########
    XV <- cbind(X, Vs, Vt)
    for (i in 1:nsim)   {
        # Ensure these are all updated from last iteration
        Sigma0_bin.inv <- as.matrix(bdiag(
            Ks_bin_inv,
            Kt_bin_inv
        ))
        T0_bin <- as.matrix(bdiag(T0a, Sigma0_bin.inv))
        T0_nb <- as.matrix(T0b)

        # Update latent variable z
        mu <- X %*% alpha + Vs %*% a + Vt %*% b
        w <- rpg(N, 1, mu[, 1])
        z <- (y1 - 1 / 2) / w

        # Update alpha, a, b
        svd_vinv <- svd(crossprod(sqrt(w) * XV) + T0_bin)
        m <- solve_svd(svd_vinv, (t(sqrt(w) * XV) %*% (sqrt(w) * z)))
        alphaab <- c(mvn_sample_svd(svd_vinv, m))
        alpha <- alphaab[1:p]
        a <- alphaab[(p + 1):(p + n)]
        b <- alphaab[-(1:(p + n))]

        # Update at-risk indicator y1 (W in paper)
        eta1 <- X %*% alpha + Vs %*% a + Vt %*% b
        eta2 <- X %*% beta # Use all n observations
        pi <- sigmoid(eta1) # at-risk probability
        q <- sigmoid(-eta2) # Pr(y=0|y1=1)
        theta <- zero_at_risk_probability(pi, q, r)
        y1[y == 0] <- rbinom(sum(y == 0), 1, theta[y == 0]) # If y=0, then draw a "chance zero" w.p. theta, otherwise y1=1
        N1 <- sum(y1)

        # Update r, TODO: Check this.
        rnew <- rtnorm(1, r, sd_r, lower = 0) # Treat r as continuous
        ratio <- sum(dnbinom(y[y1 == 1], rnew, q[y1 == 1], log = TRUE)) - sum(dnbinom(y[y1 == 1], r, q[y1 == 1], log = TRUE)) +
            dtnorm(r, rnew, sd_r, 0, log = TRUE) - dtnorm(rnew, r, sd_r, 0, log = TRUE) # Uniform Prior for R
        # Proposal not symmetric
        if (log(runif(1)) < ratio) {
            r <- rnew
        }

        # update l1t, sigma1t, noise_ratio_t1
        out <- update_ls_sigma_noise(l1t, sigma1t, noise_ratio_t1, b, Kt_bin, Dt, ltPrior, sigmaPrior, noisePrior, kern)
        l1t <- out$ls
        sigma1t <- out$sigma
        noise_ratio_t1 <- out$noise_ratio
        Kt_bin <- out$K
        Kt_bin_inv <- out$K_inv

        # update l1s, sigma1s, noise_ratio_s1
        out <- update_ls_sigma_noise(l1s, sigma1s, noise_ratio_s1, a, Ks_bin, Ds, lsPrior, sigmaPrior, noisePrior, kern)
        l1s <- out$ls
        sigma1s <- out$sigma
        noise_ratio_s1 <- out$noise_ratio
        Ks_bin <- out$K
        Ks_bin_inv <- out$K_inv

        # Update beta
        eta <- X[y1 == 1, ] %*% beta
        w <- rpg(N1, y[y1 == 1] + r, eta) # Polya weights
        z <- (y[y1 == 1] - r) / (2 * w)
        X_nb <- X[y1 == 1, , drop = FALSE]

        # Update beta
        svd_vinv <- svd(crossprod(sqrt(w) * X_nb) + T0_nb)
        m <- solve_svd(
            svd_vinv,
            t(sqrt(w) * X_nb) %*% (sqrt(w) * z)
        )
        betacd <- c(mvn_sample_svd(svd_vinv, m))
        beta <- betacd
        
        # Store
        if ((i > burn) && ((i - burn) %% thin == 0)) {
            j <- (i - burn) / thin
            
            Alpha[j, ] <- alpha
            Beta[j, ] <- beta # fixed effects

            A[j, ] <- a
            B[j, ] <- b
            
            L1t[j] <- l1t
            Noise1t[j] <- noise_ratio_t1
            Sigma1t[j] <- sigma1t # temporal hyperparameters
                
            L1s[j] <- l1s
            Noise1s[j] <- noise_ratio_s1
            Sigma1s[j] <- sigma1s # spatial hyperparameters
            
            R[j] <- r
            if (save_ypred) {
                eta_at_risk <- X %*% alpha + Vs %*% a + Vt %*% b
                eta_count <- X %*% beta
                Y_pred[j, ] <- draw_zinb_response(
                    eta_at_risk = eta_at_risk,
                    eta_count = eta_count,
                    r = r
                )
                y1s[j, ] <- y1
            }
        }
        if ((i %% print_iter == 0) && (print_progress)) message("Iteration ", i, " of ", nsim)
    }
    # Put the results into a list
    results <- list(
        Alpha = Alpha, Beta = Beta, A = A, B = B,
        L1t = L1t, Sigma1t = Sigma1t, Noise1t=Noise1t,
        L1s = L1s, Sigma1s = Sigma1s, Noise1s=Noise1s,
        R = R
    )
    if (save_ypred) {
        temp <- list(Y_pred = Y_pred, at_risk = y1s)
        results <- append(results, temp)
    }
    return(new_zinb_gp_fit(results))
}



#' ZINB_GP_count
#' @description Run the ZINB GP model with GP random effects only in the count part.
#' @inheritParams ZINB_GP_orig
#' @return A list containing posterior MCMC draws for the count GP, fixed
#'   effects, and the negative-binomial dispersion parameter.
#' @importFrom MASS glm.nb
#' @importFrom mvtnorm rmvnorm
#' @importFrom BayesLogit rpg
#' @importFrom Matrix bdiag
#' @importFrom Matrix sparseMatrix
#' @importFrom msm rtnorm
#' @importFrom msm dtnorm
#' @importFrom stats rbinom
#' @importFrom stats dnbinom
#' @importFrom stats rbinom
#' @importFrom stats glm
#' @importFrom stats rnorm
#' @importFrom LaplacesDemon rinvgamma
#' @importFrom stats runif
#' @importFrom Matrix forceSymmetric
ZINB_GP_count <- function(X, y, coords, Vs, Vt, Ds, Dt, nsim, burn, thin = 1, save_ypred = FALSE, print_iter = 100, print_progress = FALSE, ltPrior = NULL, lsPrior = NULL, sigmaPrior = NULL, noisePrior = NULL, mh_sd_r = NULL, kern = NULL) {
    validate_zinb_inputs(X, y, nsim, burn, thin)
    validate_gp_design(Vs, Ds, nrow(X), "Vs", "Ds")
    validate_gp_design(Vt, Dt, nrow(X), "Vt", "Dt")
    validate_spatial_coordinates(coords, Vs, Ds)

    # X is the design matrix with dimension N*p
    # x is the vector with length N
    # y is the count response with length N
    n <- nrow(coords) - 1 # number of clusters
    N <- nrow(X) # number of observations
    p <- ncol(X) # dimension of alpha and beta
    n_time_points <- ncol(Vt)

    # Sacrifice to the intercept gods
    Ds <- Ds[2:nrow(Ds), 2:ncol(Ds)]
    Dt <- Dt[2:nrow(Dt), 2:ncol(Dt)]

    # Use squared distances
    Ds <- Ds * Ds
    Dt <- Dt * Dt

    # Find reasonable bounds for GP length scales
    kern <- nullcheck(kern, kernel)
    param_bounds <- gp_param_bounds(Ds, Dt, kern)
    lsmax <- param_bounds$lsmax
    ltmax <- param_bounds$ltmax

    ##########
    # Priors #
    ##########

    ####### priors for alpha and beta ######
    T0a <- T0b <- diag(0.01, p)
    sd_r <- nullcheck(mh_sd_r, 0.4)
    
    ####### kernel hyperparameters  ######
    ltPrior <- nullcheck(ltPrior, list(max=ltmax, mh_sd=3, a=1, b=0.001))
    lsPrior <- nullcheck(lsPrior, list(max=lsmax, mh_sd=3, a=1, b=0.001))
    sigmaPrior <- nullcheck(sigmaPrior, list(a=0.01, b=0.1))
    noisePrior <- nullcheck(noisePrior, list(a=1.5, b=1.5, mh_sd=0.2))

    # Model init
    r <- 1
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
    p_at_risk <- sigmoid(eta1) # at-risk probability

    q <- sigmoid(-eta2)
    theta <- zero_at_risk_probability(p_at_risk, q, r)
    y1[y == 0] <- rbinom(sum(y == 0), 1, theta[y == 0]) # If y=0, then draw a "chance zero" w.p. theta, otherwise y1=1

    m1 <- glm(y1 ~ 0 + X, family = "binomial")
    alpha <- as.matrix(m1$coefficients) # initial for alpha in the binary component

    noise_ratio_t1 <- 0.5
    noise_ratio_t2 <- 0.5
    noise_ratio_s1 <- 0.5
    noise_ratio_s2 <- 0.5

    ##########################
    # Spatial Random Effects #
    ##########################
    l2s <- 1
    sigma2s <- 2
    Ks_nb <- sigma2s^2 * noise_mix(kern(Ds, l2s), noise_ratio_s2)
    Ks_nb_inv <- forceSymmetric(solve(Ks_nb))
    c <- t(rmvnorm(n = 1, sigma = Ks_nb))

    #################
    # Temporal Random Effects #
    #################
    sigma2t <- 2
    l2t <- 1
    Kt_nb <- sigma2t^2 * noise_mix(kern(Dt, l2t), noise_ratio_t2)
    Kt_nb_inv <- forceSymmetric(solve(Kt_nb))
    d <- t(rmvnorm(n = 1, sigma = Kt_nb))

    ############
    # Num Sims #
    ############
    lastit <- (nsim - burn) / thin # Last stored value

    #########
    # Store #
    #########
    Beta <- Alpha <- matrix(0, lastit, p)
    R <- rep(0, lastit)
    C <- matrix(0, lastit, n)
    D <- matrix(0, lastit, n_time_points)
    L2t <- Sigma2t <- Noise2t <- rep(0, lastit)
    L2s <- Sigma2s <- Noise2s <- rep(0, lastit)
    if (save_ypred == TRUE) {
        Y_pred <- matrix(NA, lastit, N)
        y1s <- matrix(NA, lastit, N)
    }


    ########
    # MCMC #
    ########
    XV <- cbind(X, Vs, Vt)
    for (i in 1:nsim)   {
        # Ensure these are all updated from last iteration
        Sigma0_nb.inv <- as.matrix(bdiag(
            Ks_nb_inv,
            Kt_nb_inv
        ))
        T0_bin <- as.matrix(T0a)
        T0_nb <- as.matrix(bdiag(T0b, Sigma0_nb.inv))

        # Update latent variable z
        mu <- X %*% alpha
        w <- rpg(N, 1, mu[, 1])
        z <- (y1 - 1 / 2) / w

        # Update alpha, a, b
        svd_vinv <- svd(crossprod(sqrt(w) * X) + T0_bin)
        m <- solve_svd(svd_vinv, (t(sqrt(w) * X) %*% (sqrt(w) * z)))
        alphaab <- c(mvn_sample_svd(svd_vinv, m))
        alpha <- alphaab
        
        # Update at-risk indicator y1 (W in paper)
        eta1 <- X %*% alpha
        eta2 <- X %*% beta + Vs %*% c + Vt %*% d # Use all n observations
        pi <- sigmoid(eta1) # at-risk probability
        q <- sigmoid(-eta2) # Pr(y=0|y1=1)
        theta <- zero_at_risk_probability(pi, q, r)
        y1[y == 0] <- rbinom(sum(y == 0), 1, theta[y == 0]) # If y=0, then draw a "chance zero" w.p. theta, otherwise y1=1
        N1 <- sum(y1)

        # Update r
        rnew <- rtnorm(1, r, sd_r, lower = 0) # Treat r as continuous
        ratio <- sum(dnbinom(y[y1 == 1], rnew, q[y1 == 1], log = TRUE)) - sum(dnbinom(y[y1 == 1], r, q[y1 == 1], log = TRUE)) +
            dtnorm(r, rnew, sd_r, 0, log = TRUE) - dtnorm(rnew, r, sd_r, 0, log = TRUE) # Uniform Prior for R
        # Proposal not symmetric
        if (log(runif(1)) < ratio) {
            r <- rnew
        }

        # Update beta
        eta <- X[y1 == 1, ] %*% beta + Vs[y1 == 1, ] %*% c + Vt[y1 == 1, ] %*% d
        w <- rpg(N1, y[y1 == 1] + r, eta) # Polya weights
        z <- (y[y1 == 1] - r) / (2 * w)

        # Update beta, c, d
        svd_vinv <- svd(crossprod(sqrt(w) * XV[y1 == 1, ]) + T0_nb)
        m <- solve_svd(svd_vinv, (t(sqrt(w) * XV[y1 == 1, ]) %*% (sqrt(w) * z)))
        betacd <- c(mvn_sample_svd(svd_vinv, m))
        beta <- betacd[1:p]
        c <- betacd[(p + 1):(p + n)]
        d <- betacd[-(1:(p + n))]

        # update l2t, sigma2t, noise_ratio_t2
        out <- update_ls_sigma_noise(l2t, sigma2t, noise_ratio_t2, d, Kt_nb, Dt, ltPrior, sigmaPrior, noisePrior, kern)
        l2t <- out$ls
        sigma2t <- out$sigma
        noise_ratio_t2 <- out$noise_ratio
        Kt_nb <- out$K
        Kt_nb_inv <- out$K_inv

        # update l2s, sigma2s, noise_ratio_s2
        out <- update_ls_sigma_noise(l2s, sigma2s, noise_ratio_s2, c, Ks_nb, Ds, lsPrior, sigmaPrior, noisePrior, kern)
        l2s <- out$ls
        sigma2s <- out$sigma
        noise_ratio_s2 <- out$noise_ratio
        Ks_nb <- out$K
        Ks_nb_inv <- out$K_inv

        # Store
        if ((i > burn) && ((i - burn) %% thin == 0)) {
            j <- (i - burn) / thin
            
            Alpha[j, ] <- alpha
            Beta[j, ] <- beta # fixed effects

            C[j, ] <- c
            D[j, ] <- d # random effects
            
            L2t[j] <- l2t
            Noise2t[j] <- noise_ratio_t2
            Sigma2t[j] <- sigma2t # temporal hyperparameters
            
            L2s[j] <- l2s
            Noise2s[j] <- noise_ratio_s2
            Sigma2s[j] <- sigma2s # spatial hyperparameters
            
            R[j] <- r
            if (save_ypred) {
                eta_at_risk <- X %*% alpha
                eta_count <- X %*% beta + Vs %*% c + Vt %*% d
                Y_pred[j, ] <- draw_zinb_response(
                    eta_at_risk = eta_at_risk,
                    eta_count = eta_count,
                    r = r
                )
                y1s[j, ] <- y1
            }
        }
        if ((i %% print_iter == 0) && (print_progress)) message("Iteration ", i, " of ", nsim)
    }
    # Put the results into a list
    results <- list(
        Alpha = Alpha, Beta = Beta, C = C, D = D,
        L2t = L2t, Sigma2t = Sigma2t, Noise2t=Noise2t,
        L2s = L2s, Sigma2s = Sigma2s, Noise2s=Noise2s,
        R = R
    )
    if (save_ypred) {
        temp <- list(Y_pred = Y_pred, at_risk = y1s)
        results <- append(results, temp)
    }
    return(new_zinb_gp_fit(results))
}


#' ZINB_GP_spatial
#'
#' Fit a zero-inflated negative-binomial model with spatial Gaussian-process
#' random effects in both the zero-inflation and count components.
#'
#' @param X Fixed-effect design matrix with N rows.
#' @param y Non-negative integer count response of length N.
#' @param Vs Sparse or dense spatial random-effect design matrix. It should
#'   have N rows and one column per spatial location.
#' @param Ds Spatial distance matrix with one row and column per full spatial
#'   level, including the baseline level omitted from `Vs`. Diagonal entries
#'   must be zero.
#' @param nsim Total number of MCMC iterations.
#' @param burn Number of burn-in iterations.
#' @param thin Store every thin-th iteration after burn-in.
#' @param save_ypred Whether to save posterior predictive counts and at-risk draws.
#' @param print_iter Print progress every print_iter iterations.
#' @param print_progress Whether to print MCMC progress.
#' @param lsPrior Prior and proposal controls for spatial GP length scales.
#' @param sigmaPrior Inverse-gamma prior parameters for GP variances.
#' @param noisePrior Beta prior and MH controls for GP noise ratios.
#' @param mh_sd_r Proposal standard deviation for NB dispersion r.
#' @param kern Kernel function accepting a squared-distance matrix and a
#'   length scale.
#'
#' @return A list containing posterior MCMC draws.
#' @importFrom MASS glm.nb
#' @importFrom mvtnorm rmvnorm
#' @importFrom BayesLogit rpg
#' @importFrom Matrix bdiag
#' @importFrom Matrix sparseMatrix
#' @importFrom msm rtnorm
#' @importFrom msm dtnorm
#' @importFrom stats rbinom
#' @importFrom stats dnbinom
#' @importFrom stats rbinom
#' @importFrom stats glm
#' @importFrom stats rnorm
#' @importFrom LaplacesDemon rinvgamma
#' @importFrom stats runif
#' @importFrom Matrix forceSymmetric
ZINB_GP_spatial <- function(
    X,
    y,
    Vs,
    Ds,
    nsim,
    burn,
    thin = 1,
    save_ypred = FALSE,
    print_iter = 100,
    print_progress = FALSE,
    lsPrior = NULL,
    sigmaPrior = NULL,
    noisePrior = NULL,
    mh_sd_r = NULL,
    kern = NULL
) {
    validate_zinb_inputs(X, y, nsim, burn, thin)
    validate_gp_design(Vs, Ds, nrow(X), "Vs", "Ds")

    if (length(y) != nrow(X)) {
        stop("y must have the same number of entries as rows in X.")
    }

    if (nrow(Vs) != nrow(X)) {
        stop("Vs must have the same number of rows as X.")
    }

    if (ncol(Vs) + 1 != nrow(Ds)) {
        stop("Ds must be a square matrix with ncol(Vs) + 1 rows and columns.")
    }

    if (nsim <= burn) {
        stop("nsim must be greater than burn.")
    }

    if ((nsim - burn) %% thin != 0) {
        stop("(nsim - burn) must be divisible by thin.")
    }

    if (!any(y > 0)) {
        stop("At least one positive count is required to initialize glm.nb.")
    }

    N <- nrow(X)
    p <- ncol(X)
    n_space <- ncol(Vs)
    n_saved <- (nsim - burn) / thin

    Ds <- Ds[-1, -1, drop = FALSE]
    # The original model uses squared distances in the kernel.
    Ds <- Ds^2

    kern <- nullcheck(kern, kernel)

    max_distance <- sqrt(max(Ds))
    if (!is.finite(max_distance) || max_distance <= 0) {
        max_distance <- 1
    }

    lsPrior <- nullcheck(
        lsPrior,
        list(max = max_distance, mh_sd = 3, a = 1, b = 0.001)
    )

    sigmaPrior <- nullcheck(
        sigmaPrior,
        list(a = 0.01, b = 0.1)
    )

    noisePrior <- nullcheck(
        noisePrior,
        list(a = 1.5, b = 1.5, mh_sd = 0.2)
    )

    sd_r <- nullcheck(mh_sd_r, 0.4)

    # Fixed-effect prior precisions.
    T0a <- diag(0.01, p)
    T0b <- diag(0.01, p)

    # Initial latent at-risk indicator.
    y_ind <- as.integer(y != 0)

    m1 <- stats::glm(y_ind ~ 0 + X, family = "binomial")
    alpha <- as.numeric(m1$coefficients)

    m2 <- MASS::glm.nb(y[y != 0] ~ 0 + X[y != 0, , drop = FALSE])
    beta <- as.numeric(m2$coefficients)

    eta1 <- as.numeric(X %*% alpha)
    eta2 <- as.numeric(X %*% beta)

    p_at_risk <- sigmoid(eta1)
    q <- sigmoid(-eta2)

    r <- 1

    y1 <- rep(0, N)
    y1[y > 0] <- 1

    theta <- zero_at_risk_probability(p_at_risk, q, r)

    y1[y == 0] <- stats::rbinom(
        sum(y == 0),
        size = 1,
        prob = theta[y == 0]
    )

    # Spatial GP initialization: zero-inflation component.
    l1s <- 1
    sigma1s <- 2
    noise_ratio_s1 <- 0.5

    Ks_bin <- sigma1s^2 * noise_mix(
        kern(Ds, l1s),
        noise_ratio_s1
    )

    Ks_bin_inv <- Matrix::forceSymmetric(solve(Ks_bin))
    a <- as.numeric(mvtnorm::rmvnorm(n = 1, sigma = Ks_bin))

    # Spatial GP initialization: count component.
    l2s <- 1
    sigma2s <- 2
    noise_ratio_s2 <- 0.5

    Ks_nb <- sigma2s^2 * noise_mix(
        kern(Ds, l2s),
        noise_ratio_s2
    )

    Ks_nb_inv <- Matrix::forceSymmetric(solve(Ks_nb))
    c <- as.numeric(mvtnorm::rmvnorm(n = 1, sigma = Ks_nb))

    # Posterior storage.
    Alpha <- matrix(0, n_saved, p)
    Beta <- matrix(0, n_saved, p)

    A <- matrix(0, n_saved, n_space)
    C <- matrix(0, n_saved, n_space)

    L1s <- rep(0, n_saved)
    Sigma1s <- rep(0, n_saved)
    Noise1s <- rep(0, n_saved)

    L2s <- rep(0, n_saved)
    Sigma2s <- rep(0, n_saved)
    Noise2s <- rep(0, n_saved)

    R <- rep(0, n_saved)

    if (save_ypred) {
        Y_pred <- matrix(NA_real_, n_saved, N)
        at_risk <- matrix(NA_integer_, n_saved, N)
    }

    # Fixed and spatial-effect design matrix.
    XV <- cbind(X, Vs)

    for (i in seq_len(nsim)) {
        # ---------------------------------------------------------------
        # Zero-inflation component: update alpha and spatial effect a.
        # ---------------------------------------------------------------
        T0_bin <- as.matrix(Matrix::bdiag(T0a, Ks_bin_inv))

        eta1 <- as.numeric(X %*% alpha + Vs %*% a)
        w_bin <- BayesLogit::rpg(N, 1, eta1)
        z_bin <- (y1 - 0.5) / w_bin

        svd_vinv <- svd(crossprod(sqrt(w_bin) * XV) + T0_bin)
        m <- solve_svd(
            svd_vinv,
            t(sqrt(w_bin) * XV) %*% (sqrt(w_bin) * z_bin)
        )

        draw_bin <- c(mvn_sample_svd(svd_vinv, m))
        alpha <- draw_bin[seq_len(p)]
        a <- draw_bin[p + seq_len(n_space)]

        # ---------------------------------------------------------------
        # Update the latent at-risk indicator.
        # ---------------------------------------------------------------
        eta1 <- as.numeric(X %*% alpha + Vs %*% a)
        eta2 <- as.numeric(X %*% beta + Vs %*% c)

        pi <- sigmoid(eta1)
        q <- sigmoid(-eta2)

        theta <- zero_at_risk_probability(pi, q, r)

        y1[y == 0] <- stats::rbinom(
            sum(y == 0),
            size = 1,
            prob = theta[y == 0]
        )

        N1 <- sum(y1)

        # ---------------------------------------------------------------
        # Update negative-binomial dispersion parameter r.
        # ---------------------------------------------------------------
        rnew <- msm::rtnorm(
            n = 1,
            mean = r,
            sd = sd_r,
            lower = 0
        )

        log_accept_r <-
            sum(stats::dnbinom(
                y[y1 == 1],
                size = rnew,
                prob = q[y1 == 1],
                log = TRUE
            )) -
            sum(stats::dnbinom(
                y[y1 == 1],
                size = r,
                prob = q[y1 == 1],
                log = TRUE
            )) +
            msm::dtnorm(
                x = r,
                mean = rnew,
                sd = sd_r,
                lower = 0,
                log = TRUE
            ) -
            msm::dtnorm(
                x = rnew,
                mean = r,
                sd = sd_r,
                lower = 0,
                log = TRUE
            )

        if (!is.na(log_accept_r) && log(stats::runif(1)) < log_accept_r) {
            r <- rnew
        }

        # ---------------------------------------------------------------
        # Update spatial GP for zero-inflation component.
        # ---------------------------------------------------------------
        out <- update_ls_sigma_noise(
            ls = l1s,
            sigma = sigma1s,
            noise_ratio = noise_ratio_s1,
            gpdraw = a,
            K = Ks_bin,
            D = Ds,
            lsPrior = lsPrior,
            sigmaPrior = sigmaPrior,
            noisePrior = noisePrior,
            kern = kern
        )

        l1s <- out$ls
        sigma1s <- out$sigma
        noise_ratio_s1 <- out$noise_ratio
        Ks_bin <- out$K
        Ks_bin_inv <- out$K_inv

        # ---------------------------------------------------------------
        # Count component: update beta and spatial effect c.
        # ---------------------------------------------------------------
        T0_nb <- as.matrix(Matrix::bdiag(T0b, Ks_nb_inv))

        eta_nb <- as.numeric(
            X[y1 == 1, , drop = FALSE] %*% beta +
                Vs[y1 == 1, , drop = FALSE] %*% c
        )

        w_nb <- BayesLogit::rpg(
            N1,
            y[y1 == 1] + r,
            eta_nb
        )

        z_nb <- (y[y1 == 1] - r) / (2 * w_nb)

        XV_nb <- XV[y1 == 1, , drop = FALSE]

        svd_vinv <- svd(crossprod(sqrt(w_nb) * XV_nb) + T0_nb)
        m <- solve_svd(
            svd_vinv,
            t(sqrt(w_nb) * XV_nb) %*% (sqrt(w_nb) * z_nb)
        )

        draw_nb <- c(mvn_sample_svd(svd_vinv, m))

        beta <- draw_nb[seq_len(p)]
        c <- draw_nb[p + seq_len(n_space)]

        # ---------------------------------------------------------------
        # Update spatial GP for count component.
        # ---------------------------------------------------------------
        out <- update_ls_sigma_noise(
            ls = l2s,
            sigma = sigma2s,
            noise_ratio = noise_ratio_s2,
            gpdraw = c,
            K = Ks_nb,
            D = Ds,
            lsPrior = lsPrior,
            sigmaPrior = sigmaPrior,
            noisePrior = noisePrior,
            kern = kern
        )

        l2s <- out$ls
        sigma2s <- out$sigma
        noise_ratio_s2 <- out$noise_ratio
        Ks_nb <- out$K
        Ks_nb_inv <- out$K_inv

        # ---------------------------------------------------------------
        # Store posterior draw.
        # ---------------------------------------------------------------
        if (i > burn && (i - burn) %% thin == 0) {
            j <- (i - burn) / thin

            Alpha[j, ] <- alpha
            Beta[j, ] <- beta

            A[j, ] <- a
            C[j, ] <- c

            L1s[j] <- l1s
            Sigma1s[j] <- sigma1s
            Noise1s[j] <- noise_ratio_s1

            L2s[j] <- l2s
            Sigma2s[j] <- sigma2s
            Noise2s[j] <- noise_ratio_s2

            R[j] <- r

            if (save_ypred) {
                eta1 <- as.numeric(X %*% alpha + Vs %*% a)
                eta2 <- as.numeric(X %*% beta + Vs %*% c)
                Y_pred[j, ] <- draw_zinb_response(
                    eta_at_risk = eta1,
                    eta_count = eta2,
                    r = r
                )
                at_risk[j, ] <- y1
            }
        }

        if (print_progress && i %% print_iter == 0) {
            message("Iteration ", i, " of ", nsim)
        }
    }

    results <- list(
        Alpha = Alpha,
        Beta = Beta,
        A = A,
        C = C,
        L1s = L1s,
        Sigma1s = Sigma1s,
        Noise1s = Noise1s,
        L2s = L2s,
        Sigma2s = Sigma2s,
        Noise2s = Noise2s,
        R = R
    )

    if (save_ypred) {
        results$Y_pred <- Y_pred
        results$at_risk <- at_risk
    }

    new_zinb_gp_fit(results)
}


#' ZINB_GP_spatial_count
#'
#' Fit a zero-inflated negative-binomial model with spatial Gaussian-process
#' random effects in only the count component.
#'
#' @param X Fixed-effect design matrix with N rows.
#' @param y Non-negative integer count response of length N.
#' @param Vs Sparse or dense spatial random-effect design matrix. It should
#'   have N rows and one column per spatial location.
#' @param Ds Spatial distance matrix with one row and column per full spatial
#'   level, including the baseline level omitted from `Vs`. Diagonal entries
#'   must be zero.
#' @param nsim Total number of MCMC iterations.
#' @param burn Number of burn-in iterations.
#' @param thin Store every thin-th iteration after burn-in.
#' @param save_ypred Whether to save posterior predictive counts and at-risk draws.
#' @param print_iter Print progress every print_iter iterations.
#' @param print_progress Whether to print MCMC progress.
#' @param lsPrior Prior and proposal controls for spatial GP length scales.
#' @param sigmaPrior Inverse-gamma prior parameters for GP variances.
#' @param noisePrior Beta prior and MH controls for GP noise ratios.
#' @param mh_sd_r Proposal standard deviation for NB dispersion r.
#' @param kern Kernel function accepting a squared-distance matrix and a
#'   length scale.
#'
#' @return A list containing posterior MCMC draws.
#' @importFrom MASS glm.nb
#' @importFrom mvtnorm rmvnorm
#' @importFrom BayesLogit rpg
#' @importFrom Matrix bdiag
#' @importFrom Matrix sparseMatrix
#' @importFrom msm rtnorm
#' @importFrom msm dtnorm
#' @importFrom stats rbinom
#' @importFrom stats dnbinom
#' @importFrom stats rbinom
#' @importFrom stats glm
#' @importFrom stats rnorm
#' @importFrom LaplacesDemon rinvgamma
#' @importFrom stats runif
#' @importFrom Matrix forceSymmetric
ZINB_GP_spatial_count <- function(
    X,
    y,
    Vs,
    Ds,
    nsim,
    burn,
    thin = 1,
    save_ypred = FALSE,
    print_iter = 100,
    print_progress = FALSE,
    lsPrior = NULL,
    sigmaPrior = NULL,
    noisePrior = NULL,
    mh_sd_r = NULL,
    kern = NULL
) {
    validate_zinb_inputs(X, y, nsim, burn, thin)
    validate_gp_design(Vs, Ds, nrow(X), "Vs", "Ds")

    if (length(y) != nrow(X)) {
        stop("y must have the same number of entries as rows in X.")
    }

    if (nrow(Vs) != nrow(X)) {
        stop("Vs must have the same number of rows as X.")
    }

    if (ncol(Vs) + 1 != nrow(Ds)) {
        stop("Ds must be a square matrix with ncol(Vs) + 1 rows and columns.")
    }

    if (nsim <= burn) {
        stop("nsim must be greater than burn.")
    }

    if ((nsim - burn) %% thin != 0) {
        stop("(nsim - burn) must be divisible by thin.")
    }

    if (!any(y > 0)) {
        stop("At least one positive count is required to initialize glm.nb.")
    }

    N <- nrow(X)
    p <- ncol(X)
    n_space <- ncol(Vs)
    n_saved <- (nsim - burn) / thin

    Ds <- Ds[-1, -1, drop = FALSE]
    # The original model uses squared distances in the kernel.
    Ds <- Ds^2

    kern <- nullcheck(kern, kernel)

    max_distance <- sqrt(max(Ds))
    if (!is.finite(max_distance) || max_distance <= 0) {
        max_distance <- 1
    }

    lsPrior <- nullcheck(
        lsPrior,
        list(max = max_distance, mh_sd = 3, a = 1, b = 0.001)
    )

    sigmaPrior <- nullcheck(
        sigmaPrior,
        list(a = 0.01, b = 0.1)
    )

    noisePrior <- nullcheck(
        noisePrior,
        list(a = 1.5, b = 1.5, mh_sd = 0.2)
    )

    sd_r <- nullcheck(mh_sd_r, 0.4)

    # Fixed-effect prior precisions.
    T0a <- diag(0.01, p)
    T0b <- diag(0.01, p)

    # Initial latent at-risk indicator.
    y_ind <- as.integer(y != 0)

    m1 <- stats::glm(y_ind ~ 0 + X, family = "binomial")
    alpha <- as.numeric(m1$coefficients)

    m2 <- MASS::glm.nb(y[y != 0] ~ 0 + X[y != 0, , drop = FALSE])
    beta <- as.numeric(m2$coefficients)

    eta1 <- as.numeric(X %*% alpha)
    eta2 <- as.numeric(X %*% beta)

    p_at_risk <- sigmoid(eta1)
    q <- sigmoid(-eta2)

    r <- 1

    y1 <- rep(0, N)
    y1[y > 0] <- 1

    theta <- zero_at_risk_probability(p_at_risk, q, r)

    y1[y == 0] <- stats::rbinom(
        sum(y == 0),
        size = 1,
        prob = theta[y == 0]
    )

    # Spatial GP initialization: count component.
    l2s <- 1
    sigma2s <- 2
    noise_ratio_s2 <- 0.5

    Ks_nb <- sigma2s^2 * noise_mix(
        kern(Ds, l2s),
        noise_ratio_s2
    )

    Ks_nb_inv <- Matrix::forceSymmetric(solve(Ks_nb))
    c <- as.numeric(mvtnorm::rmvnorm(n = 1, sigma = Ks_nb))

    # Posterior storage.
    Alpha <- matrix(0, n_saved, p)
    Beta <- matrix(0, n_saved, p)

    C <- matrix(0, n_saved, n_space)

    L2s <- rep(0, n_saved)
    Sigma2s <- rep(0, n_saved)
    Noise2s <- rep(0, n_saved)

    R <- rep(0, n_saved)

    if (save_ypred) {
        Y_pred <- matrix(NA_real_, n_saved, N)
        at_risk <- matrix(NA_integer_, n_saved, N)
    }

    # Fixed and spatial-effect design matrix.
    XV <- cbind(X, Vs)

    for (i in seq_len(nsim)) {
        # ---------------------------------------------------------------
        # Zero-inflation component: update alpha and spatial effect a.
        # ---------------------------------------------------------------
        T0_bin <- as.matrix(T0a)

        eta1 <- as.numeric(X %*% alpha)
        w_bin <- BayesLogit::rpg(N, 1, eta1)
        z_bin <- (y1 - 0.5) / w_bin

        svd_vinv <- svd(crossprod(sqrt(w_bin) * X) + T0_bin)
        m <- solve_svd(
            svd_vinv,
            t(sqrt(w_bin) * X) %*% (sqrt(w_bin) * z_bin)
        )

        draw_bin <- c(mvn_sample_svd(svd_vinv, m))
        alpha <- draw_bin
        
        # ---------------------------------------------------------------
        # Update the latent at-risk indicator.
        # ---------------------------------------------------------------
        eta1 <- as.numeric(X %*% alpha)
        eta2 <- as.numeric(X %*% beta + Vs %*% c)

        pi <- sigmoid(eta1)
        q <- sigmoid(-eta2)

        theta <- zero_at_risk_probability(pi, q, r)

        y1[y == 0] <- stats::rbinom(
            sum(y == 0),
            size = 1,
            prob = theta[y == 0]
        )

        N1 <- sum(y1)

        # ---------------------------------------------------------------
        # Update negative-binomial dispersion parameter r.
        # ---------------------------------------------------------------
        rnew <- msm::rtnorm(
            n = 1,
            mean = r,
            sd = sd_r,
            lower = 0
        )

        log_accept_r <-
            sum(stats::dnbinom(
                y[y1 == 1],
                size = rnew,
                prob = q[y1 == 1],
                log = TRUE
            )) -
            sum(stats::dnbinom(
                y[y1 == 1],
                size = r,
                prob = q[y1 == 1],
                log = TRUE
            )) +
            msm::dtnorm(
                x = r,
                mean = rnew,
                sd = sd_r,
                lower = 0,
                log = TRUE
            ) -
            msm::dtnorm(
                x = rnew,
                mean = r,
                sd = sd_r,
                lower = 0,
                log = TRUE
            )

        if (!is.na(log_accept_r) && log(stats::runif(1)) < log_accept_r) {
            r <- rnew
        }

        # ---------------------------------------------------------------
        # Count component: update beta and spatial effect c.
        # ---------------------------------------------------------------
        T0_nb <- as.matrix(Matrix::bdiag(T0b, Ks_nb_inv))

        eta_nb <- as.numeric(
            X[y1 == 1, , drop = FALSE] %*% beta +
                Vs[y1 == 1, , drop = FALSE] %*% c
        )

        w_nb <- BayesLogit::rpg(
            N1,
            y[y1 == 1] + r,
            eta_nb
        )

        z_nb <- (y[y1 == 1] - r) / (2 * w_nb)

        XV_nb <- XV[y1 == 1, , drop = FALSE]

        svd_vinv <- svd(crossprod(sqrt(w_nb) * XV_nb) + T0_nb)
        m <- solve_svd(
            svd_vinv,
            t(sqrt(w_nb) * XV_nb) %*% (sqrt(w_nb) * z_nb)
        )

        draw_nb <- c(mvn_sample_svd(svd_vinv, m))

        beta <- draw_nb[seq_len(p)]
        c <- draw_nb[p + seq_len(n_space)]

        # ---------------------------------------------------------------
        # Update spatial GP for count component.
        # ---------------------------------------------------------------
        out <- update_ls_sigma_noise(
            ls = l2s,
            sigma = sigma2s,
            noise_ratio = noise_ratio_s2,
            gpdraw = c,
            K = Ks_nb,
            D = Ds,
            lsPrior = lsPrior,
            sigmaPrior = sigmaPrior,
            noisePrior = noisePrior,
            kern = kern
        )

        l2s <- out$ls
        sigma2s <- out$sigma
        noise_ratio_s2 <- out$noise_ratio
        Ks_nb <- out$K
        Ks_nb_inv <- out$K_inv

        # ---------------------------------------------------------------
        # Store posterior draw.
        # ---------------------------------------------------------------
        if (i > burn && (i - burn) %% thin == 0) {
            j <- (i - burn) / thin

            Alpha[j, ] <- alpha
            Beta[j, ] <- beta

            C[j, ] <- c

            L2s[j] <- l2s
            Sigma2s[j] <- sigma2s
            Noise2s[j] <- noise_ratio_s2

            R[j] <- r

            if (save_ypred) {
                eta1 <- as.numeric(X %*% alpha)
                eta2 <- as.numeric(X %*% beta + Vs %*% c)
                Y_pred[j, ] <- draw_zinb_response(
                    eta_at_risk = eta1,
                    eta_count = eta2,
                    r = r
                )
                at_risk[j, ] <- y1
            }
        }

        if (print_progress && i %% print_iter == 0) {
            message("Iteration ", i, " of ", nsim)
        }
    }

    results <- list(
        Alpha = Alpha,
        Beta = Beta,
        C = C,
        L2s = L2s,
        Sigma2s = Sigma2s,
        Noise2s = Noise2s,
        R = R
    )

    if (save_ypred) {
        results$Y_pred <- Y_pred
        results$at_risk <- at_risk
    }

    new_zinb_gp_fit(results)
}


#' ZINB_GP_spatial_inflation
#'
#' Fit a zero-inflated negative-binomial model with spatial Gaussian-process
#' random effects in only the inflation part.
#'
#' @param X Fixed-effect design matrix with N rows.
#' @param y Non-negative integer count response of length N.
#' @param Vs Sparse or dense spatial random-effect design matrix. It should
#'   have N rows and one column per spatial location.
#' @param Ds Spatial distance matrix with one row and column per full spatial
#'   level, including the baseline level omitted from `Vs`. Diagonal entries
#'   must be zero.
#' @param nsim Total number of MCMC iterations.
#' @param burn Number of burn-in iterations.
#' @param thin Store every thin-th iteration after burn-in.
#' @param save_ypred Whether to save posterior predictive counts and at-risk draws.
#' @param print_iter Print progress every print_iter iterations.
#' @param print_progress Whether to print MCMC progress.
#' @param lsPrior Prior and proposal controls for spatial GP length scales.
#' @param sigmaPrior Inverse-gamma prior parameters for GP variances.
#' @param noisePrior Beta prior and MH controls for GP noise ratios.
#' @param mh_sd_r Proposal standard deviation for NB dispersion r.
#' @param kern Kernel function accepting a squared-distance matrix and a
#'   length scale.
#'
#' @return A list containing posterior MCMC draws.
#' @importFrom MASS glm.nb
#' @importFrom mvtnorm rmvnorm
#' @importFrom BayesLogit rpg
#' @importFrom Matrix bdiag
#' @importFrom Matrix sparseMatrix
#' @importFrom msm rtnorm
#' @importFrom msm dtnorm
#' @importFrom stats rbinom
#' @importFrom stats dnbinom
#' @importFrom stats rbinom
#' @importFrom stats glm
#' @importFrom stats rnorm
#' @importFrom LaplacesDemon rinvgamma
#' @importFrom stats runif
#' @importFrom Matrix forceSymmetric
ZINB_GP_spatial_inflation <- function(
    X,
    y,
    Vs,
    Ds,
    nsim,
    burn,
    thin = 1,
    save_ypred = FALSE,
    print_iter = 100,
    print_progress = FALSE,
    lsPrior = NULL,
    sigmaPrior = NULL,
    noisePrior = NULL,
    mh_sd_r = NULL,
    kern = NULL
) {
    validate_zinb_inputs(X, y, nsim, burn, thin)
    validate_gp_design(Vs, Ds, nrow(X), "Vs", "Ds")

    if (length(y) != nrow(X)) {
        stop("y must have the same number of entries as rows in X.")
    }

    if (nrow(Vs) != nrow(X)) {
        stop("Vs must have the same number of rows as X.")
    }

    if (ncol(Vs) + 1 != nrow(Ds)) {
        stop("Ds must be a square matrix with ncol(Vs) + 1 rows and columns.")
    }

    if (nsim <= burn) {
        stop("nsim must be greater than burn.")
    }

    if ((nsim - burn) %% thin != 0) {
        stop("(nsim - burn) must be divisible by thin.")
    }

    if (!any(y > 0)) {
        stop("At least one positive count is required to initialize glm.nb.")
    }

    N <- nrow(X)
    p <- ncol(X)
    n_space <- ncol(Vs)
    n_saved <- (nsim - burn) / thin

    Ds <- Ds[-1, -1, drop = FALSE]
    # The original model uses squared distances in the kernel.
    Ds <- Ds^2

    kern <- nullcheck(kern, kernel)

    max_distance <- sqrt(max(Ds))
    if (!is.finite(max_distance) || max_distance <= 0) {
        max_distance <- 1
    }

    lsPrior <- nullcheck(
        lsPrior,
        list(max = max_distance, mh_sd = 3, a = 1, b = 0.001)
    )

    sigmaPrior <- nullcheck(
        sigmaPrior,
        list(a = 0.01, b = 0.1)
    )

    noisePrior <- nullcheck(
        noisePrior,
        list(a = 1.5, b = 1.5, mh_sd = 0.2)
    )

    sd_r <- nullcheck(mh_sd_r, 0.4)

    # Fixed-effect prior precisions.
    T0a <- diag(0.01, p)
    T0b <- diag(0.01, p)

    # Initial latent at-risk indicator.
    y_ind <- as.integer(y != 0)

    m1 <- stats::glm(y_ind ~ 0 + X, family = "binomial")
    alpha <- as.numeric(m1$coefficients)

    m2 <- MASS::glm.nb(y[y != 0] ~ 0 + X[y != 0, , drop = FALSE])
    beta <- as.numeric(m2$coefficients)

    eta1 <- as.numeric(X %*% alpha)
    eta2 <- as.numeric(X %*% beta)

    p_at_risk <- sigmoid(eta1)
    q <- sigmoid(-eta2)

    r <- 1

    y1 <- rep(0, N)
    y1[y > 0] <- 1

    theta <- zero_at_risk_probability(p_at_risk, q, r)

    y1[y == 0] <- stats::rbinom(
        sum(y == 0),
        size = 1,
        prob = theta[y == 0]
    )

    # Spatial GP initialization: zero-inflation component.
    l1s <- 1
    sigma1s <- 2
    noise_ratio_s1 <- 0.5

    Ks_bin <- sigma1s^2 * noise_mix(
        kern(Ds, l1s),
        noise_ratio_s1
    )

    Ks_bin_inv <- Matrix::forceSymmetric(solve(Ks_bin))
    a <- as.numeric(mvtnorm::rmvnorm(n = 1, sigma = Ks_bin))

    # Posterior storage.
    Alpha <- matrix(0, n_saved, p)
    Beta <- matrix(0, n_saved, p)

    A <- matrix(0, n_saved, n_space)

    L1s <- rep(0, n_saved)
    Sigma1s <- rep(0, n_saved)
    Noise1s <- rep(0, n_saved)

    R <- rep(0, n_saved)

    if (save_ypred) {
        Y_pred <- matrix(NA_real_, n_saved, N)
        at_risk <- matrix(NA_integer_, n_saved, N)
    }

    # Fixed and spatial-effect design matrix.
    XV <- cbind(X, Vs)

    for (i in seq_len(nsim)) {
        # ---------------------------------------------------------------
        # Zero-inflation component: update alpha and spatial effect a.
        # ---------------------------------------------------------------
        T0_bin <- as.matrix(Matrix::bdiag(T0a, Ks_bin_inv))

        eta1 <- as.numeric(X %*% alpha + Vs %*% a)
        w_bin <- BayesLogit::rpg(N, 1, eta1)
        z_bin <- (y1 - 0.5) / w_bin

        svd_vinv <- svd(crossprod(sqrt(w_bin) * XV) + T0_bin)
        m <- solve_svd(
            svd_vinv,
            t(sqrt(w_bin) * XV) %*% (sqrt(w_bin) * z_bin)
        )

        draw_bin <- c(mvn_sample_svd(svd_vinv, m))
        alpha <- draw_bin[seq_len(p)]
        a <- draw_bin[p + seq_len(n_space)]

        # ---------------------------------------------------------------
        # Update the latent at-risk indicator.
        # ---------------------------------------------------------------
        eta1 <- as.numeric(X %*% alpha + Vs %*% a)
        eta2 <- as.numeric(X %*% beta)

        pi <- sigmoid(eta1)
        q <- sigmoid(-eta2)

        theta <- zero_at_risk_probability(pi, q, r)

        y1[y == 0] <- stats::rbinom(
            sum(y == 0),
            size = 1,
            prob = theta[y == 0]
        )

        N1 <- sum(y1)

        # ---------------------------------------------------------------
        # Update negative-binomial dispersion parameter r.
        # ---------------------------------------------------------------
        rnew <- msm::rtnorm(
            n = 1,
            mean = r,
            sd = sd_r,
            lower = 0
        )

        log_accept_r <-
            sum(stats::dnbinom(
                y[y1 == 1],
                size = rnew,
                prob = q[y1 == 1],
                log = TRUE
            )) -
            sum(stats::dnbinom(
                y[y1 == 1],
                size = r,
                prob = q[y1 == 1],
                log = TRUE
            )) +
            msm::dtnorm(
                x = r,
                mean = rnew,
                sd = sd_r,
                lower = 0,
                log = TRUE
            ) -
            msm::dtnorm(
                x = rnew,
                mean = r,
                sd = sd_r,
                lower = 0,
                log = TRUE
            )

        if (!is.na(log_accept_r) && log(stats::runif(1)) < log_accept_r) {
            r <- rnew
        }

        # ---------------------------------------------------------------
        # Update spatial GP for zero-inflation component.
        # ---------------------------------------------------------------
        out <- update_ls_sigma_noise(
            ls = l1s,
            sigma = sigma1s,
            noise_ratio = noise_ratio_s1,
            gpdraw = a,
            K = Ks_bin,
            D = Ds,
            lsPrior = lsPrior,
            sigmaPrior = sigmaPrior,
            noisePrior = noisePrior,
            kern = kern
        )

        l1s <- out$ls
        sigma1s <- out$sigma
        noise_ratio_s1 <- out$noise_ratio
        Ks_bin <- out$K
        Ks_bin_inv <- out$K_inv

        # ---------------------------------------------------------------
        # Count component: update beta
        # ---------------------------------------------------------------
        T0_nb <- as.matrix(T0b)

        eta_nb <- as.numeric(
            X[y1 == 1, , drop = FALSE] %*% beta)

        w_nb <- BayesLogit::rpg(
            N1,
            y[y1 == 1] + r,
            eta_nb
        )

        z_nb <- (y[y1 == 1] - r) / (2 * w_nb)
        X_nb <- X[y1 == 1, , drop = FALSE]

        svd_vinv <- svd(crossprod(sqrt(w_nb) * X_nb) + T0_nb)
        m <- solve_svd(
            svd_vinv,
            t(sqrt(w_nb) * X_nb) %*% (sqrt(w_nb) * z_nb)
        )

        draw_nb <- c(mvn_sample_svd(svd_vinv, m))

        beta <- draw_nb

        # ---------------------------------------------------------------
        # Store posterior draw.
        # ---------------------------------------------------------------
        if (i > burn && (i - burn) %% thin == 0) {
            j <- (i - burn) / thin

            Alpha[j, ] <- alpha
            Beta[j, ] <- beta

            A[j, ] <- a

            L1s[j] <- l1s
            Sigma1s[j] <- sigma1s
            Noise1s[j] <- noise_ratio_s1

            R[j] <- r

            if (save_ypred) {
                eta1 <- as.numeric(X %*% alpha + Vs %*% a)
                eta2 <- as.numeric(X %*% beta)
                Y_pred[j, ] <- draw_zinb_response(
                    eta_at_risk = eta1,
                    eta_count = eta2,
                    r = r
                )
                at_risk[j, ] <- y1
            }
        }

        if (print_progress && i %% print_iter == 0) {
            message("Iteration ", i, " of ", nsim)
        }
    }

    results <- list(
        Alpha = Alpha,
        Beta = Beta,
        A = A,
        L1s = L1s,
        Sigma1s = Sigma1s,
        Noise1s = Noise1s,
        R = R
    )

    if (save_ypred) {
        results$Y_pred <- Y_pred
        results$at_risk <- at_risk
    }

    new_zinb_gp_fit(results)
}
