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

        temp_neardistM <- diag(count_NN)
        if (count_NN > 1)
        {
            temp_neardistM[lower.tri(temp_neardistM, diag = FALSE)] <- exp(-phi * neardistM[i-1, 1:((count_NN * (count_NN - 1)) %/% 2)])
            temp_neardistM <- matrix(unlist(temp_neardistM), nrow=count_NN, ncol=count_NN)
            temp_neardistM <- t(temp_neardistM) + temp_neardistM - diag(diag(temp_neardistM))
        }
        
        L <- chol(temp_neardistM) # m by m when i > m+1

        u = exp(-phi * unlist(neardist[(i - 1), 1:count_NN]))

        v <- solve(t(L), u)
        AD[(M + 1), i] <- (1.0 - (t(v) %*% v))

        v2 <- solve(L, v)
        AD[1:count_NN, i] = v2
    }
    AD[(M + 1), 1] <- 1
    return(AD)
}