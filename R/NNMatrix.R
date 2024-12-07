#### distance matrix for location i and its neighbors ####
i_dist <- function(i, neighbor_index, s) {
    dist(s[c(i, neighbor_index[[i - 1]]), ])
}

ho_dist <- function(i, nn.ind, coords) c(dist(coords[nn.ind[i, ], ]))

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

#' Create neighbord matrix with spNNGP
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

#### Function for checking neighbors ####

Check_Neighbors <- function(coords.ord, n.neighbors, NN.matrix, ind) {
    plot(coords.ord)
    points(coords.ord[1:ind, , drop = FALSE], col = "grey", pch = 19)

    # neighbors
    if (ind < n.neighbors) {
        dim <- ind
    } else {
        dim <- n.neighbors
    }
    for (j in 1:dim) {
        points(coords.ord[NN.matrix$NN_ind[ind - 1, j], , drop = FALSE],
            col = "orange", pch = 19
        )
    }
    points(coords.ord[ind, , drop = FALSE], col = "blue", pch = 19)
    legend("topleft",
        inset = .05,
        c(
            "obs", paste0(ind, "th obs"),
            paste0("nb. of ", ind, "th obs"),
            paste0("obs <", ind)
        ), pch = c(1, 19, 19, 19),
        col = c("black", "blue", "orange", "grey"), horiz = F
    )
}
