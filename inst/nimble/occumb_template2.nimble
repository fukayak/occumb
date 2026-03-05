            }
            is_pos[j, k] <- step(sum(u[1:I, j, k]) - 0.5)
            sum_ur[j, k] <- inprod(u[1:I, j, k], r[1:I, j, k]) + (1 - is_pos[j, k])
        }
    }

    # Availability of species sequence
    for (i in 1:I) {
        for (j in 1:J) {
            for (k in 1:K) {
