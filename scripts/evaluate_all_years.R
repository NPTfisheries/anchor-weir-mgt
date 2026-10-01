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
source("./R/aha_model_funs.R")

# set parameters to calculate anchor----
# parameteris follow those produced by Kyle in his original work
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


# historical data for comparisons

scenario_dat <- dat %>%
  mutate(nor = ifelse(!is.na(Post.NO), Post.NO, In.NO),
         hor = ifelse(!is.na(Post.HO), Post.HO, In.HO)) %>%
  select(year = Year, nor, hor)

# run management frameworks against historical Losine River returns

historical_results <- purrr::pmap_dfr(
  scenario_dat,
  function(year, nor, hor) {
    
    year_inputs <- create_scenario_inputs(
      no_manarea_est = nor,  # current/scenario natural-origin abundance estimate
      ho_manarea_est = hor,   # current/scenario hatchery-origin abundance estimate
      brood_need     = 160,   # broodstock need
      spawner_goal   = 800,   # adult spawner goal
      wl_scaling     = 1.4,    # Wild-Lostine scaling factor
      utilization_no  = 1.0,   # proportion of allowed NO impacts actually taken (1.0 = full utilization)
      utilization_ho  = 1.0    # proportion of allowed HO harvest actually taken (1.0 = full utilization)
    )
    
    # year_inputs <- scenario_inputs
    # year_inputs$no_manarea_est <- nor
    # year_inputs$ho_manarea_est <- hor
    
    # anchor method
    
    anchor_inputs <- list(
      brood_need   = year_inputs$brood_need,
      spawner_goal = year_inputs$spawner_goal,
      wl_scaling   = year_inputs$wl_scaling
    )
    
    
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
    
    # recovery method
    
    recovery_targets <- list(
      Preservation = list(
        pNOB = 1.00,
        pHOS = 1.00,
        PNI  = NA_real_
      ),
      Recolonization = list(
        pNOB = 1.00,
        pHOS = 1.00,
        PNI  = NA_real_
      ),
      `Local Adaptation` = list(
        pNOB = 1.00,
        pHOS = 0.50,
        PNI  = 0.67
      ),
      `Full Restoration` = list(
        pNOB = 1.00,
        pHOS = 0.30,
        PNI  = 0.77
      )
    )
    
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
    
    # sliding scale
    
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
    
  # combine results
    
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

# truncate results to last 20 years
historical_results_20 <- historical_results %>%
  filter(
    year >= max(year, na.rm = TRUE) - 19
  )


# Phenotypic and relative fitness calculations----

# Question to ask - 
# If the same historical sequence of Lostine NOR and HOR returns had occurred
# under each management framework, how would the theoretical selection and 
# relative-fitness trajectories have differed?

# start P_nat and P_hat from a common mean phenotypic value
# initial conditions P_nat = P_hat

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
  selection_sd = 1,
  
  fitness_floor = 0
)

# NOTE:
# The original Fork model looks at generational changes with each time step, 
# but historical numbers include overlap across multiple generations. As a result
# We weighted each spawn year phenotypic value by age comp.

ford_historical_age <- historical_results_20 %>%
  group_split(method) %>%
  map_dfr(
    ~run_historical_ford_age(
      management_dat = .x,
      pars = ford_pars,
      age_comp = c(
        `3` = 0.20,
        `4` = 0.70,
        `5` = 0.10
      )
    )
  ) %>%
  arrange(
    method,
    year
  )

# calculate Ford's relative fitness based on mean phenotypic values and the AHA
# relative fitness model which introduces a floor (we set floor to 0 so both
# methods should be equal)

omega2 <- calc_omega2(
  selection_sd = ford_pars$selection_sd,
  sigma2 = ford_pars$sigma2
)

ford_historical_age <- ford_historical_age %>%
  mutate(
    
    # Natural-origin population fitness in natural environment
    fitness_nat = exp(
      -0.5 *
        (P_nat_return - ford_pars$theta_nat)^2 /
        (omega2 + ford_pars$sigma2)
    ),
    
    # Hatchery-origin population fitness in natural environment
    fitness_hat = exp(
      -0.5 *
        (P_hat_return - ford_pars$theta_nat)^2 /
        (omega2 + ford_pars$sigma2)
    ),
    fitness_return = aha_fitness(
      P_nat = P_nat_return,
      pars = ford_pars
    ),
    
    fitness_offspring = aha_fitness(
      P_nat = P_nat_offspring,
      pars = ford_pars
    )
  )


saveRDS(ford_historical_age, file = './data/output/ford_historical_age.rds')


# create figures
fig_path <- './figures/historical_evaluation/'

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
  scale_colour_viridis_d(option = "H", begin = .1, end = .9)+
  theme_bw()

ggsave(paste0(fig_path,'historical_abundance.png'))

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
  scale_colour_viridis_d(option = "A", begin = .1, end = .9)+
  labs(
    x = "Year",
    y = "Total Natural Spawners",
    color = "Management Method"
  ) +
  theme_bw()

ggsave(paste0(fig_path,'total_spawners.png'))

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
  scale_colour_viridis_d(option = "A", begin = .1, end = .9)+
  facet_wrap(
    ~origin#,
#    scales = "free_y"
  ) +
  labs(
    x = "Year",
    y = "Natural Spawners",
    color = "Management Method"
  ) +
  theme_bw()

ggsave(paste0(fig_path,'origin_spawners.png'))

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
  scale_colour_viridis_d(option = "A", begin = .1, end = .9)+
  facet_wrap(
    ~origin
  ) +
  labs(
    x = "Year",
    y = "Broodstock",
    color = "Management Method"
  ) +
  theme_bw()

ggsave(paste0(fig_path,'broodstock_take.png'))

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
  scale_colour_viridis_d(option = "A", begin = .1, end = .9)+
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

ggsave(paste0(fig_path,'hatchery_disposition.png'))

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
  scale_colour_viridis_d(option = "A", begin = .1, end = .9)+
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

ggsave(paste0(fig_path,'management_metrics.png'))


# AHA relative fitness estimate for natural population from the optimum, this
# estimate could include a fitness floor, but ours is set to 0 so the output
# matches Ford's estimates below

# ford_historical_age %>%
#   ggplot(
#     aes(
#       x = year,
#       y = fitness_return,
#       color = method
#     )
#   ) +
#   geom_line(linewidth = 1) +
#   geom_point(size = 2) +
#   scale_y_continuous(
#     limits = c(0, 1)
#   ) +
#   labs(
#     x = "Return Year",
#     y = "Relative Natural-Population Fitness",
#     color = "Management Method"
#   ) +
#   theme_bw()

# prep data for Ford's phenotypic and fitness estimates


# phenotypic response
ford_historical_age %>%
  select(
    year,
    method,
    P_nat_return,
    P_hat_return
  ) %>%
  pivot_longer(
    cols = c(
      P_nat_return,
      P_hat_return
    ),
    names_to = "origin",
    values_to = "phenotype"
  ) %>%
  mutate(
    origin = recode(
      origin,
      P_nat_return = "Natural Origin",
      P_hat_return = "Hatchery Origin"
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
  scale_colour_viridis_d(option = "A", begin = .1, end = .9)+
  facet_wrap(~origin) +
  labs(
    x = "Return Year",
    y = "Mean Phenotypic Trait",
    color = "Management Method"
  ) +
  theme_bw()

ggsave(paste0(fig_path,'phenotypic_response.png'))

fitness_plot_dat <- ford_historical_age %>%
  select(
    year,
    method,
    fitness_nat,
    fitness_hat
  ) %>%
  pivot_longer(
    cols = c(
      fitness_nat,
      fitness_hat
    ),
    names_to = "origin",
    values_to = "relative_fitness"
  ) %>%
  mutate(
    origin = recode(
      origin,
      fitness_nat = "Natural Origin",
      fitness_hat = "Hatchery Origin"
    )
  )


# relative fitness
fitness_plot_dat %>%
  ggplot(
    aes(
      x = year,
      y = relative_fitness,
      color = method
    )
  ) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_colour_viridis_d(option = "A", begin = .1, end = .9)+
  facet_wrap(~origin) +
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  labs(
    x = "Return Year",
    y = "Relative Fitness in Natural Environment",
    color = "Management Method"
  ) +
  theme_bw()

ggsave(paste0(fig_path,'relative_fitness.png'))
