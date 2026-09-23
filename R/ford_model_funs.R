# Purpose: Build functions to replicate Ford, M.J. 2002 Conservation Biology
# selection and fitness models
# Author: Ryan N. Kinzer
# Created: 09/15/2026
#
# Reference:
# Ford, M.J. 2002. Selection in captivity during supportive breeding may
# reduce fitness in the wild. Conservation Biology 16:815-825.

# Calculate omega^2
#
# Ford specifies the width of the fitness distribution and stabilizing selection
# as omega^2 while AHA specifies it in a standard deviation. Omega is
# calculated from selection strength and phenotypic variance.
#
# omega^2 is used in Fords phenotypic model equations
#
# selection_sd = width of stabilizing selection in SD units (AHA parameter)
# sigma2 = phenotypic variance (AHA parameter)
#
# omega = selection_sd * sqrt(sigma2)

calc_omega2 <- function(
    selection_sd,
    sigma2
) {
  
  omega <- selection_sd * sqrt(sigma2)
  
  omega^2
}


# Response to selection
#
# Ford 2002 Equation 4
#
# Calculates the expected mean trait value following stabilizing selection
# toward the optimum of the reproductive environment.
#
# P      = current mean phenotypic trait value
# theta  = optimum trait value in the spawning environment (natural or hatchery)
# sigma2 = phenotypic variance
# omega2 = variance describing the strength of stabilizing selection
# h2     = heritability of the trait


selection_response <- function(
    P,
    theta,
    sigma2,
    omega2,
    h2
) {
  
  selected_mean <-
    (P * omega2 + theta * sigma2) /
    (omega2 + sigma2)
  
  P + h2 * (selected_mean - P)
}


# Ford generation
#
# Ford 2002 Equations 5 and 6
#
# Calculates a mean trait value for the next generation for either the natural
# or hatchery population.
#
# pHOS controls the contribution of hatchery-origin spawners to reproduction
# in the natural environment.
#
# pNOB controls the contribution of natural-origin broodstock to reproduction
# in the hatchery environment.
#
# Similar equation in AHA Appendix C:
#   Equation 60 = natural population
#   Equation 61 = hatchery population

ford_generation <- function(
    P_nat, # mean trait value
    P_hat,
    pHOS,
    pNOB,
    pars # list of selection and fitness parameters
) {
  
  omega2 <- calc_omega2(
    selection_sd = pars$selection_sd,
    sigma2 = pars$sigma2
  )
  
  
  # Selection in natural environment
  
  # Natural population reproducing in nature
  nat_to_nat <- selection_response(
    P = P_nat,
    theta = pars$theta_nat,
    sigma2 = pars$sigma2,
    omega2 = omega2,
    h2 = pars$h2
  )
  
  # Hatchery population reproducing in nature
  hat_to_nat <- selection_response(
    P = P_hat,
    theta = pars$theta_nat,
    sigma2 = pars$sigma2,
    omega2 = omega2,
    h2 = pars$h2
  )
  
  
  # Selection in hatchery environment
  
  # Hatchery population reproducing in hatchery
  hat_to_hat <- selection_response(
    P = P_hat,
    theta = pars$theta_hat,
    sigma2 = pars$sigma2,
    omega2 = omega2,
    h2 = pars$h2
  )
  
  # Natural population reproducing in hatchery
  nat_to_hat <- selection_response(
    P = P_nat,
    theta = pars$theta_hat,
    sigma2 = pars$sigma2,
    omega2 = omega2,
    h2 = pars$h2
  )
  
  
  # Natural population in next generation
  # Equation 5
  P_nat_new <-
    (1 - pHOS) * nat_to_nat +
    pHOS * hat_to_nat
  
  
  # Hatchery population in next generation
  #@ Equation 6
  P_hat_new <-
    (1 - pNOB) * hat_to_hat +
    pNOB * nat_to_hat
  
  
  tibble(
    P_nat = P_nat_new,
    P_hat = P_hat_new
  )
}


# Proportionate natural influence
#
# HSRG/AHA management index describing the relative influence of the natural
# and hatchery environments.
#
# Returns NA when both pNOB and pHOS are zero.

pni_calc <- function(
    pNOB,
    pHOS
) {
  
  ifelse(
    pNOB + pHOS > 0,
    pNOB / (pNOB + pHOS),
    NA_real_
  )
}


# Run Ford model
#
# Runs the Ford phenotypic model for a set number of generations while holding
# holding pHOS and pNOB constant.
#
# This function describes the genetic fitness and selection dynamics only.
# Population abundance, productivity, capacity, harvest, and hatchery production
# are handled by the AHA demographic model.

ford_model <- function(
    pHOS,
    pNOB,
    pars
) {
  
  out <- tibble(
    generation = 0:pars$generations,
    P_nat = NA_real_,
    P_hat = NA_real_
  )
  
  
  # Initial trait values
  
  out$P_nat[1] <- pars$P_nat_0
  out$P_hat[1] <- pars$P_hat_0
  
  
  # Iterate generations
  
  for (g in 2:nrow(out)) {
    
    next_generation <- ford_generation(
      P_nat = out$P_nat[g - 1],
      P_hat = out$P_hat[g - 1],
      pHOS = pHOS,
      pNOB = pNOB,
      pars = pars
    )
    
    out$P_nat[g] <- next_generation$P_nat
    out$P_hat[g] <- next_generation$P_hat
  }
  
  
  # Add management metrics
  
  out %>%
    mutate(
      pHOS = pHOS,
      pNOB = pNOB,
      PNI = pni_calc(
        pNOB = pNOB,
        pHOS = pHOS
      )
    )
}
