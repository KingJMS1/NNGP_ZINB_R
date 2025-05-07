#' i_dist
#' @description distance matrix for location i and its neighbors
i_dist <- function(i, neighbor_index, s) {
    dist(s[c(i, neighbor_index[[i - 1]]), ])
}

#' get_NN_distM
get_NN_distM <- function(ind, ind_distM_d, M) {
    if (ind < M) {
        l <- ind
    } else {
        l <- M
    }
    M_i <- rep(0, M * (M - 1) / 2)
    if (l == 1) {} else {
        M_i[1:(l * (l - 1) / 2)] <-
            c(ind_distM_d[[ind]])[(l + 1):(l * (l + 1) / 2)]
    }
    return(M_i)
}

#' get_NN_dist
get_NN_dist <- function(ind, ind_distM_d, M) {
    if (ind < M) {
        l <- ind
    } else {
        l <- M
    }
    D_i <- rep(0, M)
    D_i[1:l] <- c(ind_distM_d[[ind]])[1:l]
    return(D_i)
}

#' get_NN_ind
get_NN_ind <- function(ind, ind_distM_i, M) {
    if (ind < M) {
        l <- ind
    } else {
        l <- M
    }
    D_i <- rep(0, M)
    D_i[1:l] <- c(ind_distM_i[[ind]])[1:l]
    return(D_i)
}

#' NNMatrix
#' @description Create neighbord matrix with spNNGP
#' @importFrom spNNGP spConjNNGP
NNMatrix <- function(coords, n.neighbors, n.omp.threads = 2,
                     search.type = "brute") {
    N <- nrow(coords)
    m.c <- spConjNNGP(rep(0, N) ~ 1,
        coords = coords, n.neighbors = n.neighbors,
        theta.alpha = c("phi" = 5, "alpha" = 0.5),
        k.fold = NA, n.omp.threads = n.omp.threads,
        search.type = search.type, return.neighbor.info = T,
        sigma.sq.IG = c(2, 1), cov.model = "exponential",
        verbose = F
    )
    NN_ind <- t(sapply(1:(N - 1), get_NN_ind, m.c$neighbor.info$n.indx[-1], n.neighbors))
    ord <- m.c$neighbor.info$ord
    coords.ord <- coords[ord, ]
    neighbor_dist <- sapply(2:N, i_dist, m.c$neighbor.info$n.indx[-1], coords.ord)

    NN_distM <- t(sapply(1:(N - 1), get_NN_distM, neighbor_dist, n.neighbors))
    NN_dist <- t(sapply(1:(N - 1), get_NN_dist, neighbor_dist, n.neighbors))

    return(list(
        ord = ord, coords.ord = coords.ord,
        NN_ind = NN_ind, NN_distM = NN_distM, NN_dist = NN_dist
    ))
}