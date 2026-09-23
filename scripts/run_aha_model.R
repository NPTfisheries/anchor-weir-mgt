# Purpose: Run and evaluate Ford / All-H Analyzer population model
# Author: Ryan N. Kinzer
# Created: 9/21/2026
# This script:
#
#   1. Loads Ford and AHA functions
#   2. Defines model parameters
#   3. Validates individual model components
#   4. Runs Ford genetic scenarios
#   5. Evaluates pHOS / pNOB / PNI relationships
#   6. Evaluates AHA fitness
#   7. Evaluates AHA natural population dynamics
#

# load packages and functions
library(tidyverse)
source("./R/ford_model_funs.R")
source("./R/aha_model_funs.R")

# model parameters

aha_pars <- list(

  # Ford fitness parameters
  # optimal trait value in natural environment
  theta_nat = 100,
  
  # optimal trait value in hatchery environment
  theta_hat = 80,
  
  # phenotypic variance
  sigma2 = 10,
  
  # S\strength of stabilizing selection in standard deviation units
  selection_sd = 3,
  
  # heritability
  h2 = 0.5,
  
  # Initial natural-population trait value
  P_nat_0 = 100,
  
  # Initial hatchery-population trait value
  P_hat_0 = 80,
  
  # Minimum natural relative fitness allowed by AHA
  fitness_floor = 0.50,
  
  # Population parameters
  
  # average adult-to-adult productivity
  prod_adult = 8.0,
  
  # adult recruitment capacity
  cap_adult = 800,  # set at Kyle's approximation
  
  # natural-origin fecundity
  fecundity_nor = 4000,
  
  # proportion female
  prop_female_nor = 0.50,
  
  # natural-origin smolt-to-adult survival
  sar_nor = 0.04,
  
  # unlimited capacity
  cap_unlimited = 1e6,
  
  # AHA parameters for allocating fitness across lifestages
  
  # spawner -> egg
  rel_loss_sp_egg = 0.50,
  
  # egg -> smolt
  rel_loss_egg_sm = 0.40,
  
  # Smolt -> adult
  rel_loss_sm_ad = 0.10,
  
  # hatchery parameters
  
  # relative reproductive success of HOS compared with NOS
  rrs_hos = 1.0,
  
  prespawn_survival = 0.95, # this if for hatchery environment-- should I add prespawn for natural spawnerd
  
  fecundity_hor = 4000,
  
  prop_female_hor = 0.50,
  
  surv_egg_yearling = 0.85,
  
  # ============================================================================
  # SIMULATION
  # ============================================================================
  
  generations = 10
)

# derive AHA life-stage parameters

life_pars <- derive_aha_lifestages(
  prod_adult = aha_pars$prod_adult,
  cap_adult = aha_pars$cap_adult,
  fecundity = aha_pars$fecundity_nor,
  prop_female = aha_pars$prop_female_nor,
  sar = aha_pars$sar_nor,
  cap_unlimited = aha_pars$cap_unlimited
)

life_pars


# ==============================================================================
# FORD MODEL
# ==============================================================================

# pni = .6667 at phos .3 and pnob .6 and phos = .5 and pnob = 1.0

ford_test <- ford_model(
  pHOS = 0.3,
  pNOB = 0.6,
  pars = aha_pars
) %>%
  
  mutate(
    fitness = aha_fitness(
      P_nat = P_nat,
      pars = aha_pars
    )
  )

ford_test %>%
  
  select(
    generation,
    Natural = P_nat,
    Hatchery = P_hat
  ) %>%
  
  pivot_longer(
    cols = c(Natural, Hatchery),
    names_to = "population",
    values_to = "trait"
  ) %>%
  
  ggplot(
    aes(
      x = generation,
      y = trait,
      color = population
    )
  ) +
  
  geom_line(
    linewidth = 1
  ) +
  
  geom_hline(
    yintercept = aha_pars$theta_nat,
    linetype = "dashed"
  ) +
  
  geom_hline(
    yintercept = aha_pars$theta_hat,
    linetype = "dashed"
  ) +
  
  labs(
    x = "Generation",
    y = "Mean Trait Value",
    color = "Population"
  ) +
  
  theme_bw()


# Plot natural relative fitness -------------------------------------------------

ggplot(
  ford_test,
  aes(
    x = generation,
    y = fitness
  )
) +
  
  geom_line(
    linewidth = 1
  ) +
  
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  
  labs(
    x = "Generation",
    y = "Natural Relative Fitness"
  ) +
  
  theme_bw()


# ==============================================================================
# MANAGEMENT SCENARIOS
# ==============================================================================

phos_range <- seq(.1, .5, by = .1)
pnob_range <- seq(.25,1,.25)

management_scenarios <- expand.grid('pHOS' = phos_range,
                                    'pNOB' = pnob_range) %>%
  mutate(
    scenario = 1:n(),
    PNI = round(pni_calc(pNOB = pNOB, pHOS = pHOS), 2)
  )

management_scenarios

scenario_results <- management_scenarios %>%
  
  mutate(
    model = map2(
      pHOS,
      pNOB,
      ~ ford_model(
        pHOS = .x,
        pNOB = .y,
        pars = aha_pars
      )
    )
  ) %>%
  
  select(
    scenario,
    model
  ) %>%
  unnest(model) %>%
  mutate(
    fitness = aha_fitness(
      P_nat = P_nat,
      pars = aha_pars
    ),
    PNI = round(PNI, 2),
    label = paste0('PNI = ', PNI, ', pNOB = ',pNOB,', pHOS = ', pHOS),
    label = forcats::fct_reorder(
      label,
      PNI,
      .fun = first,
      .desc = TRUE
    )
  )


# Trait trajectories among management scenarios --------------------------------

ggplot(
  scenario_results,
  aes(
    x = generation,
    y = P_nat,
    color = label #as.factor(scenario)
  )
) +
  
  geom_line(
    linewidth = 1
  ) +
  
  geom_hline(
    yintercept = aha_pars$theta_nat,
    linetype = "dashed"
  ) +
  geom_label(data = scenario_results[scenario_results$generation==9,],
             aes(x = generation,
                 y = P_nat,
                 label = PNI)) +
  labs(
    x = "Generation",
    y = "Natural Population Trait Value",
    color = "Scenario"
  ) +
  
  theme_bw()


# Fitness trajectories among management scenarios ------------------------------

ggplot(
  scenario_results,
  aes(
    x = generation,
    y = fitness,
    color = label #as.factor(scenario)
  )
) +
  
  geom_line(
    linewidth = 1
  ) +
  geom_label(data = scenario_results[scenario_results$generation == 9,],
             aes(x = generation,
                 y = fitness,
                 label = PNI)) +
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  
  labs(
    title = "Ford's (2002) Relative Fitness",
    subtitle =  "Relative fitness reflects the reduction in fitness as the population's,
    mean phenotype moves away from the natural optimum.",
    x = "Generation",
    y = "Natural Relative Fitness",
    color = "Scenario"
  ) +
  
  theme_bw()


# ==============================================================================
# pHOS x pNOB EQUILIBRIUM SURFACE
# ==============================================================================


scenario_grid <- expand_grid(
  pHOS = seq(
    0.01,
    1.00,
    by = 0.01
  ),
  pNOB = seq(
    0.01,
    1.00,
    by = 0.01
  )
)


ford_grid <- scenario_grid %>%
  
  mutate(
    model = map2(
      pHOS,
      pNOB,
      ~ ford_model(
        pHOS = .x,
        pNOB = .y,
        pars = aha_pars
      )
    )
  )


ford_equilibrium <- ford_grid %>%
  
  mutate(
    equilibrium = map(
      model,
      ~ slice_tail(
        .x,
        n = 1
      )
    )
  ) %>%
  
  select(
    pHOS,
    pNOB,
    equilibrium
  ) %>%
  
  unnest(
    equilibrium,
    names_sep = "_"
  ) %>%
  
  mutate(
    PNI = pni_calc(
      pNOB = pNOB,
      pHOS = pHOS
    ),
    
    fitness = aha_fitness(
      P_nat = equilibrium_P_nat,
      pars = aha_pars
    )
  )


# PNI surface ------------------------------------------------------------------

ggplot(
  ford_equilibrium, aes(x = pNOB, y = pHOS)) +
  geom_raster(aes(fill = PNI)) +
  geom_contour(aes(z = PNI),
    breaks = c(
      0.50,
      0.67,
      0.80
    )
  ) +
  coord_equal() +
  scale_fill_viridis_c(limits = c(0, 1)) +
  labs(
    x = "pNOB",
    y = "pHOS",
    fill = "PNI"
  ) +
  theme_bw()


# Ford/AHA fitness surface ------------------------------------------------------

ggplot(
  ford_equilibrium,
  aes(
    x = pNOB,
    y = pHOS)) +
  geom_raster(aes(fill = fitness)) +
  geom_contour(aes(z = PNI),
    breaks = c(
      0.50,
      0.67,
      0.80
    )
  ) +
  coord_equal() +
  scale_fill_viridis_c(
    limits = c(aha_pars$fitness_floor, 1),
    name = "Relative\nFitness"
  ) +
  labs(
    x = "pNOB",
    y = "pHOS"
  ) +
  
  theme_bw()


# ==============================================================================
# SAME PNI, DIFFERENT pHOS / pNOB
# ==============================================================================


same_pni <- tribble(
  ~scenario, ~pHOS, ~pNOB,
  "A",         0.05,  0.10,
  "B",         0.10,  0.20,
  "C",         0.20,  0.40,
  "D",         0.30,  0.60,
  "E",         0.40,  0.80
) %>%
  mutate(
    PNI = pni_calc(pNOB = pNOB, pHOS = pHOS)
  )


same_pni_results <- same_pni %>%
  mutate(
    model = map2(
      pHOS,
      pNOB,
      ~ ford_model(
        pHOS = .x,
        pNOB = .y,
        pars = aha_pars
      )
    ),
    equilibrium = map(
      model,
      ~ slice_tail(
        .x,
        n = 1
      )
    )
  ) %>%
  select(
    scenario,
    pHOS,
    pNOB,
    PNI,
    equilibrium
  ) %>%
  unnest(
    equilibrium,
    names_sep = "_"
  ) %>%
  mutate(
    fitness = aha_fitness(
      P_nat = equilibrium_P_nat,
      pars = aha_pars
    )
  )

same_pni_results %>%
  select(
    scenario,
    pHOS,
    pNOB,
    PNI,
    equilibrium_P_nat,
    equilibrium_P_hat,
    fitness
  )


# ==============================================================================
# VALIDATE AHA NATURAL POPULATION DYNAMICS
# ==============================================================================


# Natural production with no HOS and no fitness reduction ----------------------

natural_test <- aha_natural_production(
  NOS = 400,
  HOS = 0,
  fitness = 1,
  pars = aha_pars
)

natural_test

# loop over the multiple generations using the predicted return
n_gen <- 50

NOS <- 400

natural_results <- vector("list", n_gen)

for (g in seq_len(n_gen)) {
  
  # Run natural production for this generation
  result <- aha_natural_production(
    NOS = NOS,
    HOS = 0,
    fitness = 1,
    pars = aha_pars
  )
  
  # Save results
  natural_results[[g]] <- tibble(
    generation = g,
    spawners = NOS,
    eggs = result$eggs[1],
    smolts = result$smolts[1],
    adults = result$adults[1]
  )
  
  # Adults become next generation's spawners
  NOS <- result$adults[1]
}

natural_results <- bind_rows(natural_results)

ggplot(
  natural_results,
  aes(
    x = generation,
    y = adults
  )
) +
  geom_line(linewidth = 1) +
  geom_point() +
  labs(
    x = "Generation",
    y = "Natural-Origin Adults",
    title = "Natural Population Dynamics"
  ) +
  theme_bw()


# Adult stock-recruit relationship ---------------------------------------------

sr_test <- tibble(
  spawners = seq(
    1,
    2000,
    by = 10
  )
) %>%
  
  mutate(
    adults = map_dbl(
      spawners,
      ~ aha_natural_production(
        NOS = .x,
        HOS = 0,
        fitness = 1,
        pars = aha_pars
      )$adults["NOS"]
    )
  )


ggplot(
  sr_test,
  aes(
    x = spawners,
    y = adults
  )
) +
  
  geom_line(
    linewidth = 1
  ) +
  
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  
  labs(
    x = "Natural Spawners",
    y = "Adult Recruits"
  ) +
  
  theme_bw()

# ==============================================================================
# EFFECT OF FORD FITNESS ON NATURAL PRODUCTION
# ==============================================================================

fitness_sr <- expand_grid(
  spawners = seq(1, 2000, by = 10),
  fitness = c(1.00, 0.90, 0.75, 0.50)
) %>%
  mutate(
    adults = map2_dbl(
      spawners,
      fitness,
      ~ aha_natural_production(
        NOS = .x,
        HOS = 0,
        fitness = .y,
        pars = aha_pars
      )$adults["NOS"]
    ),
    fitness = factor(
      fitness,
      levels = c(1.00, 0.90, 0.75, 0.50)
    )
  )


ggplot(
  fitness_sr,
  aes(
    x = spawners,
    y = adults,
    color = fitness
  )
) +
  geom_line(linewidth = 1) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title = "Fitness-Adjusted Beverton-Holt Recruitment",
    subtitle = paste(
      "Relative fitness modifies the productivity and capacity",
      "of the natural population."
    ),
    x = "Natural Spawners",
    y = "Adult Recruits",
    color = "Relative Fitness"
  ) +
  theme_bw()


# ==============================================================================
# pHOS / EFFECTIVE pHOS EXAMPLE
# ==============================================================================


aha_calc_phos(
  NOS = 600,
  HOS = 200,
  rrs_hos = aha_pars$rrs_hos
)


# ==============================================================================
# HATCHERY PRODUCTION EXAMPLE
# ==============================================================================


hatchery_test <- aha_hatchery_releases(
  NOB = 100,
  HOB = 60,
  pars = aha_pars
)

hatchery_test


# ==============================================================================
# HARVEST EXAMPLE
# ==============================================================================


NOR_harvest <- aha_harvest(
  adults = 1000,
  
  harvest_rates = c(
    0.10, # below Bonneville Dam
    0.05, # Zone 6 fishery
    0.05, # Snake River
    0.10 # terminal
  ),
  
  adult_passage_survival = 1 # passage survival between fisheries 3 and 4
)

NOR_harvest

HOR_harvest <- aha_harvest(
  adults = 1000,
  
  harvest_rates = c(
    0.20,
    0.10,
    0.10,
    0.30
  ),
  
  adult_passage_survival = 1
)

HOR_harvest

