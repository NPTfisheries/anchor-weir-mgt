# Purpose: This script iterates across each year of data and runs three weir
# management frameworks with deterministic functions to evaluate Lostine River
# demographic and fitness responses.
#
# Author: Ryan N. Kinzer
# Date: 2026-09-30

# load packages----
library(tidyverse)

source("./R/helpers.R")
source("./R/harvest_funs.R")
source("./R/weir_funs.R")
source("./R/sliding_scale_funs.R")
source("./R/anchor_funs.R")
source("./R/recovery_funs.R")
source("./R/ford_model_funs.R")

# set parameters to calculate anchor----
input_file <- "input_chs.csv"
file_path <- file.path('./data',input_file)

trap <- 'LST'

dat <- read_csv(file_path, show_col_types = FALSE) %>%
  filter(.data$Stock == trap) %>%
  arrange(desc(.data$Year))

# create trap ratio columns safely
dat <- dat %>%
  mutate(
    trap_ratio_no = if_else(Post.NO == 0, NA_real_, Trap.NO / Post.NO), # trap_ratio is the number trapped out of total return
    trap_ratio_ho = if_else(Post.HO == 0, NA_real_, Trap.HO / Post.HO)
  )

m_trap_no <- first5_non_na(dat$trap_ratio_no, dat$Year)
m_trap_ho <- first5_non_na(dat$trap_ratio_ho, dat$Year)
m_weir    <- first5_non_na(dat$Weir.Eff,     dat$Year)
m_abv     <- first5_non_na(dat$PSS.ABV,      dat$Year)
m_blw     <- first5_non_na(dat$PSS.BLW,      dat$Year)


median_params <- list(
  trap_prop_no        = median(m_trap_no$Value), # median proportion of NO captured at the weir; is this equal to weir efficiency?
  trap_prop_ho        = median(m_trap_ho$Value),
  weir_efficiency     = median(m_weir$Value),
  survival_above_weir = median(m_abv$Value),
  survival_below_weir = median(m_blw$Value),
  years_trap_no       = years_used_str(m_trap_no),
  years_trap_ho       = years_used_str(m_trap_ho),
  years_weir          = years_used_str(m_weir),
  years_abv           = years_used_str(m_abv),
  years_blw           = years_used_str(m_blw)
)

# default accounting parameters are the median-derived values.
# These are used for normal planning/review scenarios when current-year estimates
# are unavailable or when use_csv_scenario_values is FALSE.
accounting_params <- list(
  trap_prop_no        = median_params$trap_prop_no,
  trap_prop_ho        = median_params$trap_prop_ho,
  weir_efficiency     = median_params$weir_efficiency,
  survival_above_weir = median_params$survival_above_weir,
  survival_below_weir = median_params$survival_below_weir,
  source_trap_no      = paste0("Median years: ", median_params$years_trap_no),
  source_trap_ho      = paste0("Median years: ", median_params$years_trap_ho),
  source_weir         = paste0("Median years: ", median_params$years_weir),
  source_abv          = paste0("Median years: ", median_params$years_abv),
  source_blw          = paste0("Median years: ", median_params$years_blw)
)

# run scenarios----

scenario_dat <- dat %>%
  mutate(nor = ifelse(!is.na(Post.NO), Post.NO, In.NO),
         hor = ifelse(!is.na(Post.HO), Post.HO, In.HO)) %>%
  select(year = Year, nor, hor)

# ============================================================
# Historical management scenario comparison
# ============================================================

historical_results <- purrr::pmap_dfr(
  scenario_dat,
  function(year, nor, hor) {
    
    # --------------------------------------------------------
    # Create annual scenario inputs
    # --------------------------------------------------------
    
    year_inputs <- scenario_inputs
    
    year_inputs$no_manarea_est <- nor
    year_inputs$ho_manarea_est <- hor
    
    
    # --------------------------------------------------------
    # Anchor
    # --------------------------------------------------------
    
    anchor_result <- run_accounting_scenario(
      scenario_inputs = year_inputs,
      accounting_params = accounting_params,
      anchor_inputs = anchor_inputs
    )
    
    anchor_out <- tidy_accounting_result(
      anchor_result
    ) %>%
      mutate(
        year = year,
        nor = nor,
        hor = hor,
        method = "Anchor"
      )
    
    
    # --------------------------------------------------------
    # Recovery
    # --------------------------------------------------------
    
    recovery_result <- run_recovery_scenario(
      scenario_inputs = year_inputs,
      accounting_params = accounting_params,
      recovery_abundance = nor,
      cbp_low = 1000,
      cbp_medium = 2500,
      cbp_high = 4000,
      recovery_targets = recovery_targets
    )
    
    recovery_out <- tidy_recovery_result(
      recovery_result
    ) %>%
      mutate(
        year = year,
        nor = nor,
        hor = hor,
        method = "Recovery"
      )
    
    
    # --------------------------------------------------------
    # Sliding Scale
    # --------------------------------------------------------
    
    sliding_result <- run_sliding_scale_scenario(
      scenario_inputs = year_inputs,
      accounting_params = accounting_params
    )
    
    sliding_out <- tidy_sliding_scale_result(
      sliding_result
    ) %>%
      mutate(
        year = year,
        nor = nor,
        hor = hor,
        method = "Sliding Scale"
      )
    
    
    # --------------------------------------------------------
    # Combine annual results
    # --------------------------------------------------------
    
    bind_rows(
      anchor_out,
      recovery_out,
      sliding_out
    )
  }
) %>%
  arrange(
    year,
    method
  )

# truncate results
historical_results_20 <- historical_results %>%
  filter(
    year >= max(year, na.rm = TRUE) - 19
  )

#-------------------------------------------------------
# create figures - truncate results to last 20 years
# true abundance
historical_results_20 %>%
  select(
    year,
    nor,
    hor
  ) %>%
  distinct() %>%
  pivot_longer(
    cols = c(
      nor,
      hor
    ),
    names_to = "origin",
    values_to = "adults"
  ) %>%
  mutate(
    origin = recode(
      origin,
      nor = "Natural Origin",
      hor = "Hatchery Origin"
    )
  ) %>%
  ggplot(
    aes(
      x = year,
      y = adults,
      color = origin
    )
  ) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    x = "Year",
    y = "Adult Returns",
    color = "Origin"
  ) +
  theme_bw()

# natural spawners
historical_results_20 %>%
  ggplot(
    aes(
      x = year,
      y = system_spawners,
      color = method
    )
  ) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    x = "Year",
    y = "Total Natural Spawners",
    color = "Management Method"
  ) +
  theme_bw()


# origin specific spawners
historical_results_20 %>%
  select(
    year,
    method,
    NO_spawners,
    HO_spawners
  ) %>%
  pivot_longer(
    cols = c(
      NO_spawners,
      HO_spawners
    ),
    names_to = "origin",
    values_to = "spawners"
  ) %>%
  mutate(
    origin = recode(
      origin,
      NO_spawners = "Natural-Origin Spawners",
      HO_spawners = "Hatchery-Origin Spawners"
    )
  ) %>%
  ggplot(
    aes(
      x = year,
      y = spawners,
      color = method
    )
  ) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(
    ~origin,
    scales = "free_y"
  ) +
  labs(
    x = "Year",
    y = "Natural Spawners",
    color = "Management Method"
  ) +
  theme_bw()

# brood stock
historical_results_20 %>%
  select(
    year,
    method,
    NO_brood,
    HO_brood
  ) %>%
  pivot_longer(
    cols = c(
      NO_brood,
      HO_brood
    ),
    names_to = "origin",
    values_to = "brood"
  ) %>%
  mutate(
    origin = recode(
      origin,
      NO_brood = "Natural-Origin Broodstock",
      HO_brood = "Hatchery-Origin Broodstock"
    )
  ) %>%
  ggplot(
    aes(
      x = year,
      y = brood,
      color = method
    )
  ) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(
    ~origin
  ) +
  labs(
    x = "Year",
    y = "Broodstock",
    color = "Management Method"
  ) +
  theme_bw()

# management metrics
historical_results_20 %>%
  select(
    year,
    method,
    pNOB,
    pHOS,
    PNI
  ) %>%
  pivot_longer(
    cols = c(
      pNOB,
      pHOS,
      PNI
    ),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c(
        "pNOB",
        "pHOS",
        "PNI"
      )
    )
  ) %>%
  ggplot(
    aes(
      x = year,
      y = value,
      color = method
    )
  ) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(
    ~metric,
    ncol = 1
  ) +
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  labs(
    x = "Year",
    y = NULL,
    color = "Management Method"
  ) +
  theme_bw()

# hatchery disposition
historical_results_20 %>%
  select(
    year,
    method,
    total_HO_harvest,
    HO_removed,
    HO_brood,
    HO_spawners
  ) %>%
  pivot_longer(
    cols = c(
      total_HO_harvest,
      HO_removed,
      HO_brood,
      HO_spawners
    ),
    names_to = "disposition",
    values_to = "fish"
  ) %>%
  mutate(
    disposition = recode(
      disposition,
      total_HO_harvest = "Harvest",
      HO_removed = "Removed at Weir",
      HO_brood = "Broodstock",
      HO_spawners = "Natural Spawning"
    )
  ) %>%
  ggplot(
    aes(
      x = year,
      y = fish,
      color = method
    )
  ) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(
    ~disposition,
    scales = "free_y"
  ) +
  labs(
    x = "Year",
    y = "Hatchery-Origin Adults",
    color = "Management Method"
  ) +
  theme_bw()

# Question to ask!!!
#If the same historical sequence of Lostine NOR and HOR returns had occurred under each management framework, how would the theoretical selection and relative-fitness trajectories have differed?

# start P_nat and P_hat from common initial conditions P_nat = P_hat = Phat

# ============================================================
# Historical Ford model
#
# Uses realized annual pNOB and pHOS from each management
# framework to update the Ford phenotypic model.
#
# NOTE:
# Each return year is treated as one Ford model step for this
# initial diagnostic. Model steps should not yet be interpreted
# as biological generations.
# ============================================================


# ------------------------------------------------------------
# Ford parameters
# ------------------------------------------------------------

ford_pars <- list(
  
  # Initial population mean phenotypes
  P_nat_0 = 0,
  P_hat_0 = 0,
  
  # Environmental optima
  theta_nat = 0,
  theta_hat = 1,
  
  # Phenotypic variance
  sigma2 = 1,
  
  # Heritability
  h2 = 0.5,
  
  # Width of stabilizing selection in SD units
  selection_sd = 1
)

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


# ============================================================
# Run Ford model across historical management scenarios
# ============================================================

ford_historical <- historical_results_20 %>%
  group_split(method) %>%
  map_dfr(
    ~run_historical_ford(
      management_dat = .x,
      pars = ford_pars
    )
  ) %>%
  arrange(
    method,
    year
  )


# ============================================================
# Calculate Ford relative fitness
# ============================================================

omega2 <- calc_omega2(
  selection_sd = ford_pars$selection_sd,
  sigma2 = ford_pars$sigma2
)

ford_historical <- ford_historical %>%
  mutate(
    
    # Relative fitness entering annual management step
    fitness_in = exp(
      -0.5 *
        (P_nat_in - ford_pars$theta_nat)^2 /
        (omega2 + ford_pars$sigma2)
    ),
    
    # Relative fitness following annual management step
    fitness_out = exp(
      -0.5 *
        (P_nat_out - ford_pars$theta_nat)^2 /
        (omega2 + ford_pars$sigma2)
    )
  )

ford_historical %>%
  select(
    year,
    method,
    P_nat_out,
    P_hat_out
  ) %>%
  pivot_longer(
    cols = c(
      P_nat_out,
      P_hat_out
    ),
    names_to = "population",
    values_to = "phenotype"
  ) %>%
  mutate(
    population = recode(
      population,
      P_nat_out = "Natural Population",
      P_hat_out = "Hatchery Population"
    )
  ) %>%
  ggplot(
    aes(
      x = year,
      y = phenotype,
      color = method
    )
  ) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~population) +
  geom_hline(
    data = tibble(
      population = c(
        "Natural Population",
        "Hatchery Population"
      ),
      optimum = c(
        ford_pars$theta_nat,
        ford_pars$theta_hat
      )
    ),
    aes(
      yintercept = optimum
    ),
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  labs(
    title = "Ford (2002) Phenotypic Response",
    subtitle = paste(
      "Historical Lostine River returns with management-specific",
      "pNOB and pHOS"
    ),
    x = "Return Year",
    y = "Mean Phenotypic Trait",
    color = "Management Method"
  ) +
  theme_bw()

# relative fitness
ford_historical %>%
  ggplot(
    aes(
      x = year,
      y = fitness_out,
      color = method
    )
  ) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  labs(
    title = "Ford's (2002) Relative Fitness",
    subtitle = paste(
      "Relative fitness reflects the population's mean phenotype",
      "relative to the optimum phenotype in the natural environment."
    ),
    x = "Return Year",
    y = "Natural Relative Fitness",
    color = "Management Method"
  ) +
  theme_bw()
