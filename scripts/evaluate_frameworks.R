# Purpose: This script uses the three weir management frameworks with
# deterministic functions to evaluate Lostine River demographic and fitness
# responses. Functions were originally developed by Ford and AHA.
#
# Author: Ryan N. Kinzer
# Date: 2026-09-24

# load packages----
library(tidyverse)

source("./R/helpers.R")
source("./R/harvest_funs.R")
source("./R/weir_funs.R")
source("./R/sliding_scale_funs.R")
source("./R/anchor_funs.R")
source("./R/recovery_funs.R")

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


# validate management frameworks at varying abundance levels

# function to run one return year abundance

run_management_test <- function(
    no_abundance,
    ho_abundance,
    abundance_type
) {
  
  # Update adult abundance.
  
  test_inputs <- create_scenario_inputs(
    no_manarea_est = no_abundance,  # current/scenario natural-origin abundance estimate
    ho_manarea_est = ho_abundance,   # current/scenario hatchery-origin abundance estimate
    brood_need     = 160,   # broodstock need
    spawner_goal   = 800,   # adult spawner goal
    wl_scaling     = 1.4,    # Wild-Lostine scaling factor
    utilization_no  = 1.0,   # proportion of allowed NO impacts actually taken (1.0 = full utilization)
    utilization_ho  = 1.0    # proportion of allowed HO harvest actually taken (1.0 = full utilization)
  )
  
# anchor method
  
  anchor_inputs <- list(
    brood_need   = test_inputs$brood_need,
    spawner_goal = test_inputs$spawner_goal,
    wl_scaling   = test_inputs$wl_scaling
  )
  
  anchor_test <- run_accounting_scenario(
    scenario_inputs = test_inputs,
    accounting_params = accounting_params,
    anchor_inputs = anchor_inputs
  )
  
  anchor_out <- tidy_accounting_result(
    anchor_test
  ) %>%
    mutate(
      method = "Anchor"
    )
  
 # recovery framework
  
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
  
  recovery_test <- run_recovery_scenario(
    scenario_inputs = test_inputs,
    accounting_params = accounting_params,
    recovery_abundance = no_abundance,
    cbp_low = 1000,
    cbp_medium = 2500,
    cbp_high = 4000,
    recovery_targets = recovery_targets
  )
  
  recovery_out <- tidy_recovery_result(
    recovery_test
  ) %>%
    mutate(
      method = "Recovery"
    )
  
  
 # sliding scale
  
  sliding_test <- run_sliding_scale_scenario(
    scenario_inputs = test_inputs,
    accounting_params = accounting_params
  )
  
  sliding_out <- tidy_sliding_scale_result(
    sliding_test
  ) %>%
    mutate(
      method = "Sliding Scale"
    )
  
  
# combine results
  
  bind_rows(
    anchor_out,
    recovery_out,
    sliding_out
  ) %>%
    mutate(
      abundance_type = abundance_type,
      abundance = if_else(
        abundance_type == "NOR varies",
        no_abundance,
        ho_abundance
      ),
      NO_abundance = no_abundance,
      HO_abundance = ho_abundance
    )
}

# validation for each origin
# 1. NOR abundance varies 0-2000; HOR held at the average return value
# 2. HOR abundance varies 0-2000; NOR held at the average return value

abundance_test <- seq(0, 2000, by = 10)

# NOR validate
#
# NOR = 0-2000
# HOR = mean(dat$Post.HO, na.rm = TRUE) # 800

nor_test_results <- purrr::map_dfr(
  abundance_test,
  function(no_abundance) {
    
    run_management_test(
      no_abundance = no_abundance,
      ho_abundance = mean(dat$Post.HO, na.rm = TRUE),
      abundance_type = "NOR varies"
    )
  }
)

# HOR validation
#
# HOR = 0-2000
# NOR = mean(dat$Post.NO, na.rm = TRUE) # 369

hor_test_results <- purrr::map_dfr(
  abundance_test,
  function(ho_abundance) {
    
    run_management_test(
      no_abundance = mean(dat$Post.NO, na.rm = TRUE),
      ho_abundance = ho_abundance,
      abundance_type = "HOR varies"
    )
  }
)


 # combine results across both origins

test_results <- bind_rows(
  nor_test_results,
  hor_test_results
) %>%
  mutate(
    abundance_type = factor(
      abundance_type,
      levels = c(
        "NOR varies",
        "HOR varies"
      )
    )
  )


# figures

# create figures
fig_path <- './figures/validation/'

plot_results <- test_results %>%
  mutate(
    abundance_type = recode(
      as.character(abundance_type),
      `NOR varies` = "NOR varies; HOR = 800",
      `HOR varies` = "HOR varies; NOR = 369"
    ),
    abundance_type = factor(
      abundance_type,
      levels = c(
        "NOR varies; HOR = 800",
        "HOR varies; NOR = 369"
      )
    )
  )

# total spawners

plot_results %>%
  ggplot(
    aes(
      x = abundance,
      y = system_spawners,
      color = method
    )
  ) +
  geom_line(linewidth = 1) +
  scale_colour_viridis_d(option = "A", begin = .1, end = .9)+
  facet_wrap(
    ~abundance_type,
    scales = "free_y"
  ) +
  labs(
    x = "Adult Abundance",
    y = "Total Natural Spawners",
    color = "Management Method"
  ) +
  theme_bw()


ggsave(paste0(fig_path,'val_total_spawners.png'))


# origin spawners

plot_results %>%
  select(
    abundance,
    abundance_type,
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
      x = abundance,
      y = spawners,
      color = method
    )
  ) +
  geom_line(linewidth = 1) +
  scale_colour_viridis_d(option = "A", begin = .1, end = .9)+
  facet_grid(
    origin ~ abundance_type,
    scales = "free_y"
  ) +
  labs(
    x = "Adult Abundance",
    y = "Natural Spawners",
    color = "Management Method"
  ) +
  theme_bw()

ggsave(paste0(fig_path,'val_origin_spawners.png'))

# broodstock

plot_results %>%
  select(
    abundance,
    abundance_type,
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
      x = abundance,
      y = brood,
      color = method
    )
  ) +
  geom_line(linewidth = 1) +
  scale_colour_viridis_d(option = "A", begin = .1, end = .9)+
  facet_grid(
    origin ~ abundance_type
  ) +
  labs(
    x = "Adult Abundance",
    y = "Broodstock",
    color = "Management Method"
  ) +
  theme_bw()

ggsave(paste0(fig_path,'val_broodstock.png'))

# management metrics

plot_results %>%
  select(
    abundance,
    abundance_type,
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
      x = abundance,
      y = value,
      color = method
    )
  ) +
  geom_line(linewidth = 1) +
  scale_colour_viridis_d(option = "A", begin = .1, end = .9)+
  facet_grid(
    metric ~ abundance_type
  ) +
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  labs(
    x = "Adult Abundance",
    y = NULL,
    color = "Management Method"
  ) +
  theme_bw()

ggsave(paste0(fig_path,'val_management_metrics.png'))

# harvest

plot_results %>%
  select(
    abundance,
    abundance_type,
    total_NO_impact,
    total_HO_harvest
  ) %>%
  distinct() %>%
  pivot_longer(
    cols = c(
      total_NO_impact,
      total_HO_harvest
    ),
    names_to = "origin",
    values_to = "harvest"
  ) %>%
  mutate(
    origin = recode(
      origin,
      total_NO_impact = "Natural-Origin Impact",
      total_HO_harvest = "Hatchery-Origin Harvest"
    )
  ) %>%
  ggplot(
    aes(
      x = abundance,
      y = harvest
    )
  ) +
  geom_line(linewidth = 1) +
  facet_grid(
    origin ~ abundance_type,
    scales = "free_y"
  ) +
  labs(
    x = "Adult Abundance",
    y = "Fishery Harvest"
  ) +
  theme_bw()

ggsave(paste0(fig_path,'val_harvest.png'))

# ============================================================
# Sport versus treaty fishery
# ============================================================

harvest_plot_dat <- plot_results %>%
  select(
    abundance,
    abundance_type,
    sport_NO_impact,
    treaty_NO_impact,
    HO_sport_harvest,
    HO_treaty_harvest
  ) %>%
  distinct() %>%
  transmute(
    abundance,
    abundance_type,
    `NO - Sport` = sport_NO_impact,
    `NO - Treaty` = treaty_NO_impact,
    `HO - Sport` = HO_sport_harvest,
    `HO - Treaty` = HO_treaty_harvest
  ) %>%
  pivot_longer(
    cols = -c(
      abundance,
      abundance_type
    ),
    names_to = "group",
    values_to = "harvest"
  ) %>%
  separate(
    group,
    into = c(
      "origin",
      "fishery"
    ),
    sep = " - "
  ) %>%
  mutate(
    origin = recode(
      origin,
      NO = "Natural-Origin Impact",
      HO = "Hatchery-Origin Harvest"
    )
  )


harvest_plot_dat %>%
  ggplot(
    aes(
      x = abundance,
      y = harvest,
      color = fishery
    )
  ) +
  geom_line(linewidth = 1) +
  facet_grid(
    origin ~ abundance_type,
    scales = "free_y"
  ) +
  labs(
    x = "Adult Abundance",
    y = "Fishery Mortality / Harvest",
    color = "Fishery"
  ) +
  theme_bw()
