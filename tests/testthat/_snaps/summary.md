# summary() works for occumbData as expected

    Code
      summary(data)
    Output
      Sequence read counts: 
       Number of species, I = 2 
       Number of sites, J = 2 
       Maximum number of replicates per site, K = 2 
       Number of missing observations = 0 
       Number of replicates per site: 2 (average), 0 (sd) 
       Sequencing depth: 9 (average), 2.9 (sd) 
      
      Species covariates: 
       cov1 (continuous) 
      Site covariates: 
       cov2 (continuous), cov3 (categorical) 
      Replicate covariates: 
       cov4 (continuous) 
      
      Labels for species: 
       (None) 
      Labels for sites: 
       (None) 
      Labels for replicates: 
       (None) 

# summary() works for occumbFit as expected

    Code
      summary(occumb:::internal_fit)
    Output
      Summary for an occumbFit object 
      
      Summary of data:
       Number of species, I = 2 
       Number of sites, J = 2 
       Maximum number of replicates per site, K = 2 
       Number of missing observations = 0 
       Number of replicates per site: 2 (average), 0 (sd) 
       Sequencing depth: 9 (average), 2.9 (sd) 
      
      Model specification:
       formula_phi:          ~ 1 
       formula_theta:        ~ 1 
       formula_psi:          ~ 1 
       formula_phi_shared:   ~ 1 
       formula_theta_shared: ~ 1 
       formula_psi_shared:   ~ 1 
       prior_prec:           1e-04 
       prior_ulim:           10000 
      
      Saved parameters:
       Mu sigma rho alpha beta gamma phi theta psi z pi deviance 
      
      MCMC ran for 0.001 minutes at time 2024-10-21 16:28:26.034164:
       For each of 1 chains:
        Adaptation:            100 iterations (sufficient)
        Burn-in:               10 iterations
        Thin rate:             1 iterations
        Total chain length:    120 iterations
        Posterior sample size: 10 draws
      
      Summary of posterior samples: 
       Mu: 
        Number of parameters: 3 
        Rhat:  (not available) 
        n.eff: (not available) 
       sigma: 
        Number of parameters: 3 
        Rhat:  (not available) 
        n.eff: (not available) 
       rho: 
        Number of parameters: 3 
        Rhat:  (not available) 
        n.eff: (not available) 
       alpha: 
        Number of parameters: 2 
        Rhat:  (not available) 
        n.eff: (not available) 
       beta: 
        Number of parameters: 2 
        Rhat:  (not available) 
        n.eff: (not available) 
       gamma: 
        Number of parameters: 2 
        Rhat:  (not available) 
        n.eff: (not available) 
       phi: 
        Number of parameters: 2 
        Rhat:  (not available) 
        n.eff: (not available) 
       theta: 
        Number of parameters: 2 
        Rhat:  (not available) 
        n.eff: (not available) 
       psi: 
        Number of parameters: 2 
        Rhat:  (not available) 
        n.eff: (not available) 
       z: 
        Number of parameters: 4 
        Rhat:  (not available) 
        n.eff: (not available) 
       pi: 
        Number of parameters: 8 
        Rhat:  (not available) 
        n.eff: (not available) 
       deviance: 
        Rhat:  (not available) 
        n.eff: (not available) 
      
      **WARNING** Set n.chains > 1 to monitor Rhat and n.eff values. 

