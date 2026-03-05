
        for (m in 1:M_phi) {
            alpha[i, m] <- spec_eff[i, m_phi[m]]
        }
        for (m in 1:M_theta) {
            beta[i, m]  <- spec_eff[i, m_theta[m]]
        }
        for (m in 1:M_psi) {
            gamma[i, m] <- spec_eff[i, m_psi[m]]
        }

        spec_eff[i, 1:M] ~ dmnorm(Mu[1:M], cov = Sigma[1:M, 1:M])
    }

    # Priors
    for (m in 1:M) {
        Mu[m] ~ dnorm(0, prior_prec)
    }

    for (m in 1:M) {
        Sigma[m, m] <- pow(sigma[m], 2)
    }
    for (m1 in 1:(M - 1)) {
        for (m2 in (m1 + 1):M) {
            Sigma[m1, m2] <- rho[rho_index[m1, m2]] * sigma[m1] * sigma[m2]
            Sigma[m2, m1] <- Sigma[m1, m2]
            rho[rho_index[m1, m2]] ~ dunif(-1, 1)
        }
    }

    for (m in 1:M) {
        sigma[m] ~ dunif(0, prior_ulim)
    }

