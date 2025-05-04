getAD <- function(neardist, neardistM, N, M, phi) {
    # coords: n by 2 array,
    # phi, v, alpha: covariance parameter set
    # hold: whether AD is for holded locations or not

    # return: A D
    AD <- matrix(NA, nrow = M + 1, ncol = N)
    # for every row of matrix A and D
    for (i in 2:N) { # i starts from 2 because first row of A = 0
        # the number of neighbors of i
        count_NN <- if (i < M + 1) {
            i - 1
        } else {
            M
        }

        # distance among the neighbors
        temp_neardistM <- matrix(NA, nrow = count_NN, ncol = count_NN)
        u <- rep(NA, count_NN)

        if (count_NN == 1) {
            temp_neardistM[1, 1] <- 1
        } else {
            h <- 0
            for (j in 1:(count_NN - 1)) {
                for (k in (j + 1):count_NN) {
                    h <- h + 1
                    temp_neardistM[j, k] <- exp(-phi * neardistM[(i - 1), h])
                    temp_neardistM[k, j] <- temp_neardistM[j, k]
                }
            }
            for (j in 1:count_NN) {
                temp_neardistM[j, j] <- 1
            }
        }

        L <- chol(temp_neardistM) # m by m when i > m+1

        for (j in 1:count_NN) {
            u[j] <- exp(-phi * neardist[(i - 1), j])
        }

        v <- solve(t(L)) %*% u
        # print(dim(v))
        AD[(M + 1), i] <- (1.0 - (t(v) %*% v))

        v2 <- solve(L) %*% v
        for (j in 1:count_NN) {
            AD[j, i] <- v2[j]
        }
        # A[row.indx, col.indx] = solve(L) %*% v # a vector
    }
    AD[(M + 1), 1] <- 1
    return(AD)
}