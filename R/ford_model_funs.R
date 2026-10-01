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

run_historical_ford <- function(
    management_dat,
    pars
) {
  
  # Ford recursion must proceed chronologically.
  management_dat <- management_dat %>%
    arrange(year)
  
  
  # ----------------------------------------------------------
  # Output object
  # ----------------------------------------------------------
  
  out <- management_dat %>%
    mutate(
      P_nat_in = NA_real_,
      P_hat_in = NA_real_,
      P_nat_out = NA_real_,
      P_hat_out = NA_real_
    )
  
  
  # ----------------------------------------------------------
  # Initial phenotypic state
  # ----------------------------------------------------------
  
  P_nat_current <- pars$P_nat_0
  P_hat_current <- pars$P_hat_0
  
  
  # ----------------------------------------------------------
  # Iterate across annual management conditions
  # ----------------------------------------------------------
  
  for (i in seq_len(nrow(out))) {
    
    # Phenotype entering this year's management event.
    out$P_nat_in[i] <- P_nat_current
    out$P_hat_in[i] <- P_hat_current
    
    
    # --------------------------------------------------------
    # Apply this year's realized pHOS and pNOB
    # --------------------------------------------------------
    
    if (
      !is.na(out$pHOS[i]) &&
      !is.na(out$pNOB[i])
    ) {
      
      next_state <- ford_generation(
        P_nat = P_nat_current,
        P_hat = P_hat_current,
        pHOS = out$pHOS[i],
        pNOB = out$pNOB[i],
        pars = pars
      )
      
      P_nat_next <- next_state$P_nat
      P_hat_next <- next_state$P_hat
      
    } else {
      
      # If management metrics are undefined, no Ford update
      # is applied for this model step.
      P_nat_next <- P_nat_current
      P_hat_next <- P_hat_current
    }
    
    
    # --------------------------------------------------------
    # Store resulting phenotype
    # --------------------------------------------------------
    
    out$P_nat_out[i] <- P_nat_next
    out$P_hat_out[i] <- P_hat_next
    
    
    # --------------------------------------------------------
    # Carry state into next model step
    # --------------------------------------------------------
    
    P_nat_current <- P_nat_next
    P_hat_current <- P_hat_next
  }
  
  
  # ----------------------------------------------------------
  # Return results
  # ----------------------------------------------------------
  
  out
}

# with age
# ============================================================
# Age-structured historical Ford model
#
# Applies annual management-specific pHOS and pNOB to returning
# adults and assigns the resulting offspring phenotype to future
# return years using a fixed age-at-return distribution.
#
# Pre-model brood years are initialized at P_nat_0 and P_hat_0.
# ============================================================

run_historical_ford_age <- function(
    management_dat,
    pars,
    age_comp = c(`3` = 0.20, `4` = 0.70, `5` = 0.10)
) {
  
  # ----------------------------------------------------------
  # Check inputs
  # ----------------------------------------------------------
  
  required_cols <- c(
    "year",
    "pHOS",
    "pNOB"
  )
  
  missing_cols <- setdiff(
    required_cols,
    names(management_dat)
  )
  
  if (length(missing_cols) > 0) {
    stop(
      paste(
        "management_dat is missing:",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  if (abs(sum(age_comp) - 1) > 1e-8) {
    stop(
      "age_comp must sum to 1.",
      call. = FALSE
    )
  }
  
  ages <- as.numeric(
    names(age_comp)
  )
  
  
  # ----------------------------------------------------------
  # Arrange management data chronologically
  # ----------------------------------------------------------
  
  management_dat <- management_dat %>%
    arrange(year)
  
  min_year <- min(
    management_dat$year
  )
  
  max_year <- max(
    management_dat$year
  )
  
  
  # ----------------------------------------------------------
  # Create phenotype history
  #
  # We need brood years before the first observed return year
  # because age-3 through age-5 adults returning in the first
  # modeled year were produced before the analysis begins.
  # ----------------------------------------------------------
  
  brood_years <- seq(
    min_year - max(ages),
    max_year
  )
  
  phenotype_history <- tibble(
    brood_year = brood_years,
    P_nat_offspring = NA_real_,
    P_hat_offspring = NA_real_
  )
  
  
  # ----------------------------------------------------------
  # Initialize pre-model brood years
  # ----------------------------------------------------------
  
  phenotype_history <- phenotype_history %>%
    mutate(
      P_nat_offspring = if_else(
        brood_year < min_year,
        pars$P_nat_0,
        P_nat_offspring
      ),
      
      P_hat_offspring = if_else(
        brood_year < min_year,
        pars$P_hat_0,
        P_hat_offspring
      )
    )
  
  
  # ----------------------------------------------------------
  # Output object
  # ----------------------------------------------------------
  
  out <- management_dat %>%
    mutate(
      P_nat_return = NA_real_,
      P_hat_return = NA_real_,
      P_nat_offspring = NA_real_,
      P_hat_offspring = NA_real_
    )
  
  
  # ----------------------------------------------------------
  # Iterate across return years
  # ----------------------------------------------------------
  
  for (i in seq_len(nrow(out))) {
    
    current_year <- out$year[i]
    
    
    # --------------------------------------------------------
    # Identify brood years producing this year's adults
    # --------------------------------------------------------
    
    source_brood_years <-
      current_year - ages
    
    
    # --------------------------------------------------------
    # Get natural-origin parental phenotypes
    # --------------------------------------------------------
    
    nat_source <- phenotype_history %>%
      filter(
        brood_year %in% source_brood_years
      ) %>%
      arrange(
        match(
          brood_year,
          source_brood_years
        )
      )
    
    
    # --------------------------------------------------------
    # Get hatchery-origin parental phenotypes
    # --------------------------------------------------------
    
    hat_source <- phenotype_history %>%
      filter(
        brood_year %in% source_brood_years
      ) %>%
      arrange(
        match(
          brood_year,
          source_brood_years
        )
      )
    
    
    # --------------------------------------------------------
    # Verify that all required brood years exist
    # --------------------------------------------------------
    
    if (
      nrow(nat_source) != length(ages) ||
      nrow(hat_source) != length(ages)
    ) {
      
      stop(
        paste(
          "Missing brood-year phenotype for return year",
          current_year
        ),
        call. = FALSE
      )
    }
    
    
    # --------------------------------------------------------
    # Calculate phenotype of adults returning this year
    #
    # Example:
    #
    # 0.20 * age-3 brood phenotype +
    # 0.70 * age-4 brood phenotype +
    # 0.10 * age-5 brood phenotype
    # --------------------------------------------------------
    
    P_nat_return <- sum(
      age_comp *
        nat_source$P_nat_offspring
    )
    
    P_hat_return <- sum(
      age_comp *
        hat_source$P_hat_offspring
    )
    
    
    # Store returning phenotype
    out$P_nat_return[i] <-
      P_nat_return
    
    out$P_hat_return[i] <-
      P_hat_return
    
    
    # --------------------------------------------------------
    # Apply this year's management
    # --------------------------------------------------------
    
    if (
      !is.na(out$pHOS[i]) &&
      !is.na(out$pNOB[i])
    ) {
      
      next_state <- ford_generation(
        P_nat = P_nat_return,
        P_hat = P_hat_return,
        pHOS = out$pHOS[i],
        pNOB = out$pNOB[i],
        pars = pars
      )
      
      P_nat_offspring <-
        next_state$P_nat
      
      P_hat_offspring <-
        next_state$P_hat
      
    } else {
      
      # If management metrics are unavailable, retain the
      # returning phenotype for the offspring state.
      
      P_nat_offspring <-
        P_nat_return
      
      P_hat_offspring <-
        P_hat_return
    }
    
    
    # --------------------------------------------------------
    # Store offspring phenotype
    # --------------------------------------------------------
    
    out$P_nat_offspring[i] <-
      P_nat_offspring
    
    out$P_hat_offspring[i] <-
      P_hat_offspring
    
    
    # --------------------------------------------------------
    # Add this brood year's offspring phenotype to history
    # --------------------------------------------------------
    
    phenotype_history$P_nat_offspring[
      phenotype_history$brood_year == current_year
    ] <- P_nat_offspring
    
    phenotype_history$P_hat_offspring[
      phenotype_history$brood_year == current_year
    ] <- P_hat_offspring
  }
  
  
  # ----------------------------------------------------------
  # Return results
  # ----------------------------------------------------------
  
  out
}
