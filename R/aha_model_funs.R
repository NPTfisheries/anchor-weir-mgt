# Purpose: Population-dynamics functions for the All-H Analyzer (AHA).
# Author: Ryan N. Kinzer
# Created: 9/16/2026
#
# The AHA demographic model links the Ford 2002 selection models with demographic
# responses using life-state productivity and capacity parameters. Other
# functions exists to calculate pHOS, pNOB, hatchery production, and harvest.
# 
# Ford trait recursion is implemented in AHA Appendix C, Equations 60-61.

source('./R/ford_model_funs.R')

# AHA relative fitness
#
# Transforms the mean trait value from Ford into a relative
# fitness for the natural population spawning in the natural environment.
# How well adapted is the natural population to the natural environment?

# P_nat = mean trait value of natural population
# AHA includes a minimum fitness value as a fitness_floor.

aha_fitness <- function(
    P_nat,
    pars
) {
  
  omega2 <- calc_omega2(
    selection_sd = pars$selection_sd,
    sigma2 = pars$sigma2
  )
  
  fitness <- exp(
    -0.5 *
      (P_nat - pars$theta_nat)^2 /
      (omega2 + pars$sigma2)
  )
  
  pmax(
    fitness,
    pars$fitness_floor
  )
}


# AHA Beverton-Holt relationship:
#
# N_next = (N * productivity) /
#          (1 + N * productivity / capacity)

aha_bh <- function(
    N,
    productivity,
    capacity
) {
  
  (N * productivity) /
    (1 + (N * productivity / capacity))
}


# Derive AHA life-stage parameters
#
# Converts adult-scale productivity and capacity into three life stages:
#
#   spawner -> egg
#   egg -> smolt
#   smolt -> adult
#
# productivity is the maximum recruitment at low levels of abundance
#
# Spawner-to-egg productivity is determined by fecundity and sex ratio. 
# The max number of eggs that can contribute to future generations.
#
# Smolt-to-adult productivity is SAR.
#
# egg-to-smolt productivity and capacity are derived so the combined
# life cycle corresponds to the specified adult-scale productivity/capacity.

derive_aha_lifestages <- function(
    prod_adult,
    cap_adult,
    fecundity,
    prop_female,
    sar,
    cap_unlimited = 1e6 # capacity of emergent egg and smolt
) {
  
  # Spawner -> emergence
  
  prod_sp_egg <-
    fecundity * prop_female
  
  cap_sp_egg <-
    cap_unlimited
  
  
  # Smolt -> adult
  
  prod_sm_ad <-
    sar
  
  cap_sm_ad <-
    cap_unlimited
  
  
  # egg -> smolt
  
  prod_egg_sm <-
    prod_adult /
    (prod_sp_egg * prod_sm_ad)
  
  
  # Derive emergence-to-smolt capacity
  
  cap_egg_sm <-
    1 / (
      prod_sm_ad *
        (
          (1 / cap_adult) -
            (1 / cap_sm_ad)
        )
    )
  
  
  tibble(
    stage = c(
      "spawner_egg",
      "egg_smolt",
      "smolt_adult"
    ),
    
    productivity = c(
      prod_sp_egg,
      prod_egg_sm,
      prod_sm_ad
    ),
    
    capacity = c(
      cap_sp_egg,
      cap_egg_sm,
      cap_sm_ad
    )
  )
}


# Allocate fitness among life stages
#
# AHA allocates the total fitness effect among life stages:
#
# stage fitness = total fitness ^ relative loss
#
# The loss parameters should sum to 1.

allocate_fitness <- function(
    fitness,
    pars
) {
  
  tibble(
    stage = c(
      "spawner_egg",
      "egg_smolt",
      "smolt_adult"
    ),
    
    rel_loss = c(
      pars$rel_loss_sp_egg,
      pars$rel_loss_egg_sm,
      pars$rel_loss_sm_ad
    )
  ) %>%
    
    mutate(
      fitness_multiplier =
        fitness^rel_loss
    )
}


# Apply fitness to life-stage parameters
#
# AHA modifies both productivity and capacity using the stage-specific
# fitness multiplier.

apply_fitness_to_lifestages <- function(
    life_pars,
    fitness,
    pars
) {
  
  fitness_pars <-
    allocate_fitness(
      fitness = fitness,
      pars = pars
    )
  
  life_pars %>%
    
    left_join(
      fitness_pars,
      by = "stage"
    ) %>%
    
    mutate(
      fit_productivity =
        productivity * fitness_multiplier,
      
      fit_capacity =
        capacity * fitness_multiplier
    )
}


# Calculate pHOS
#
# Calculates both census pHOS and a effective pHOS
#
# Effective pHOS weights HOS by relative reproductive success in the same way
# as our proposal.

aha_calc_phos <- function(
    NOS,
    HOS,
    rrs_hos = 1
) {
  
  HOS_effective <-
    HOS * rrs_hos
  
  
  pHOS_census <-
    ifelse(
      NOS + HOS > 0,
      HOS / (NOS + HOS),
      NA_real_
    )
  
  
  pHOS_effective <-
    ifelse(
      NOS + HOS_effective > 0,
      HOS_effective /
        (NOS + HOS_effective),
      NA_real_
    )
  
  
  tibble(
    HOS_effective = HOS_effective,
    pHOS_census = pHOS_census,
    pHOS_effective = pHOS_effective
  )
}


# Calculate pNOB

aha_calc_pnob <- function(
    NOB,
    HOB
) {
  
  ifelse(
    NOB + HOB > 0,
    NOB / (NOB + HOB),
    NA_real_
  )
}


# Single AHA life stage
#
# Applies a shared density-dependent Beverton-Holt denominator to multiple
# origin components.
#
# Example:
#
# N_components = c(NOS = 500, HOS = 100)
#
# The components retain their identities while experiencing the same
# density-dependent environment.

aha_stage <- function(
    N_components,
    productivity,
    capacity
) {
  
  N_total <-
    sum(N_components)
  
  denominator <-
    1 +
    productivity * N_total /
    capacity
  
  N_next <-
    N_components *
    productivity /
    denominator
  
  N_next
}


# Natural production
#
# Moves NOS and effective HOS through the three AHA life stages:
#
#   spawner -> egg -> smolt -> adult
#
# Fitness modifies productivity and capacity at each stage.

aha_natural_production <- function(
    NOS,
    HOS,
    fitness,
    pars
) {
  
  # Baseline life-stage parameters
  
  life_pars <-
    derive_aha_lifestages(
      prod_adult = pars$prod_adult,
      cap_adult = pars$cap_adult,
      fecundity = pars$fecundity_nor,
      prop_female = pars$prop_female_nor,
      sar = pars$sar_nor,
      cap_unlimited = pars$cap_unlimited
    )
  
  # Apply genetic fitness effect
  
  life_pars <-
    apply_fitness_to_lifestages(
      life_pars = life_pars,
      fitness = fitness,
      pars = pars
    )
  
  # Effective HOS
  
  HOS_effective <-
    HOS * pars$rrs_hos
  
  spawners <- c(
    NOS = NOS,
    HOS = HOS_effective
  )
  
  # spawner -> egg
  
  stage_1 <-
    life_pars %>%
    filter(stage == "spawner_egg")
  
  egg <-
    aha_stage(
      N_components = spawners,
      productivity = stage_1$fit_productivity,
      capacity = stage_1$fit_capacity
    )
  
  
  # egg -> smolt
  
  stage_2 <-
    life_pars %>%
    filter(stage == "egg_smolt")
  
  smolts <-
    aha_stage(
      N_components = egg,
      productivity = stage_2$fit_productivity,
      capacity = stage_2$fit_capacity
    )
  
  
  # smolt -> adult
  
  stage_3 <-
    life_pars %>%
    filter(stage == "smolt_adult")
  
  adults <-
    aha_stage(
      N_components = smolts,
      productivity = stage_3$fit_productivity,
      capacity = stage_3$fit_capacity
    )
  
  
  list(
    life_pars = life_pars,
    spawners = spawners,
    egg = egg,
    smolts = smolts,
    adults = adults
  )
}


# Hatchery production
#
# Converts broodstock abundance to hatchery juvenile releases.
#
# NOB      = natural-origin broodstock
# HOB      = hatchery-origin broodstock
# imported = broodstock imported from another population/program

aha_hatchery_releases <- function(
    NOB,
    HOB,
    pars
) {
  
  eggs_per_spawner <-
    pars$prespawn_survival *
    pars$fecundity_hor *
    pars$prop_female_hor
  
  brood_total <- NOB + HOB
  
  eggs <-brood_total * eggs_per_spawner
  
  releases <- eggs * pars$surv_egg_yearling
  
  tibble(
    brood_total = brood_total,
    eggs_per_spawner = eggs_per_spawner,
    eggs = eggs,
    releases = releases
  )
}

# AHA permits four sequential fisheries/harvest as fish move upstream
#
# Fisheries 1-3 occur before adult passage mortality.
# Fishery 4 occurs after adult passage mortality.
#
# NOR and HOR can be passed through this function separately using
# origin-specific harvest rates.

aha_harvest <- function(
    adults,
    harvest_rates,
    adult_passage_survival = 1
) {
  
  stopifnot(
    length(harvest_rates) == 4
  )
  
  remaining <-
    adults
  
  harvest <-
    numeric(4)
  
  
  # Fisheries 1-3
  
  for (i in 1:3) {
    
    harvest[i] <-
      remaining *
      harvest_rates[i]
    
    remaining <-
      remaining -
      harvest[i]
  }
  
  
  # Adult passage survival
  
  remaining <-
    remaining *
    adult_passage_survival
  
  
  # Terminal fishery
  
  harvest[4] <-
    remaining *
    harvest_rates[4]
  
  remaining <-
    remaining -
    harvest[4]
  
  
  tibble(
    fishery_1 = harvest[1],
    fishery_2 = harvest[2],
    fishery_3 = harvest[3],
    fishery_4 = harvest[4],
    total_harvest = sum(harvest),
    escapement = remaining
  )
}

# ==============================================================================
# AHA ADULT ALLOCATION
#
# Purpose:
#   Allocate post-harvest adult NOR and HOR escapement among broodstock,
#   natural spawning, export, and hatchery surplus.
#
# Inputs:
#   NOR_esc              = NOR escapement after harvest and adult passage
#   HOR_esc              = HOR escapement after harvest and adult passage
#   local_brood          = local broodstock requirement
#   export_brood_goal    = goal for exported hatchery-origin broodstock
#   pNOB_goal            = management target for proportion natural-origin brood
#   pHOS_goal            = management target for proportion hatchery-origin
#                          natural spawners
#   max_prop_nor_brood   = maximum fraction of NOR escapement available
#                          for broodstock
#   prop_hatchery        = proportion of HOR escapement returning to hatchery
#   prop_spawn           = proportion of HOR escapement going to natural spawn
#   total_strays         = additional hatchery-origin strays entering the
#                          natural spawning population
#   rrs_hos              = relative reproductive success of HOS
#
# Outputs:
#   NOB                  = natural-origin broodstock
#   HOB                  = hatchery-origin broodstock
#   NOS                  = natural-origin natural spawners
#   HOS_local            = locally returning HOR natural spawners
#   HOS                  = total hatchery-origin natural spawners
#   HOS_effective        = RRS-adjusted HOS
#   pNOB                 = realized proportion natural-origin broodstock
#   pHOS_census          = census pHOS
#   pHOS_effective       = RRS-adjusted pHOS used by Ford model
#   PNI                  = pNOB / (pNOB + effective pHOS)
#
# Notes:
#   AHA uses census HOS to calculate the management pHOS goal, but uses
#   RRS-adjusted HOS when calculating the pHOS value entering the Ford
#   fitness model.
#
# ==============================================================================

aha_allocate_adults <- function(
    NOR_esc,
    HOR_esc,
    local_brood,
    export_brood_goal = 0,
    pNOB_goal = 0,
    pHOS_goal = 0,
    max_prop_nor_brood = 1,
    prop_hatchery = 1,
    prop_spawn = 0,
    total_strays = 0,
    rrs_hos = 1
) {
  
  # ============================================================================
  # 1. INPUT CHECKS
  # ============================================================================
  
  NOR_esc <-
    max(0, NOR_esc)
  
  HOR_esc <-
    max(0, HOR_esc)
  
  local_brood <-
    max(0, local_brood)
  
  export_brood_goal <-
    max(0, export_brood_goal)
  
  total_strays <-
    max(0, total_strays)
  
  
  pNOB_goal <-
    min(
      max(pNOB_goal, 0),
      1
    )
  
  pHOS_goal <-
    min(
      max(pHOS_goal, 0),
      1
    )
  
  max_prop_nor_brood <-
    min(
      max(max_prop_nor_brood, 0),
      1
    )
  
  prop_hatchery <-
    min(
      max(prop_hatchery, 0),
      1
    )
  
  prop_spawn <-
    min(
      max(prop_spawn, 0),
      1
    )
  
  rrs_hos <-
    max(0, rrs_hos)
  
  
  # ============================================================================
  # 2. NATURAL-ORIGIN BROODSTOCK
  #
  # NOB is constrained by:
  #
  #   1. the pNOB management goal,
  #   2. the maximum allowable proportion of NOR escapement,
  #   3. available NOR escapement, and
  #   4. the total local brood requirement.
  #
  # Scenario 2 generation 2:
  #
  #   pNOB target:
  #       627 * 0.30 = 188.1
  #
  #   maximum NOR brood:
  #       301 * 0.30 = 90.3
  #
  # Therefore:
  #
  #       NOB = 90.3
  #
  # ============================================================================
  
  NOB_goal <-
    local_brood *
    pNOB_goal
  
  
  NOB_max <-
    NOR_esc *
    max_prop_nor_brood
  
  
  NOB <-
    min(
      NOB_goal,
      NOB_max,
      NOR_esc,
      local_brood
    )
  
  
  # ============================================================================
  # 3. NATURAL-ORIGIN NATURAL SPAWNERS
  # ============================================================================
  
  NOS <-
    max(
      0,
      NOR_esc - NOB
    )
  
  
  # ============================================================================
  # 4. HOR ESCAPEMENT DISPOSITION
  #
  # Post-harvest HOR escapement is divided between:
  #
  #   - adults returning to the hatchery
  #   - adults entering natural spawning
  #
  # External strays are treated separately below.
  #
  # Scenario 2:
  #
  #   90% -> hatchery
  #   10% -> natural spawning
  #
  # ============================================================================
  
  HOR_hatchery <-
    HOR_esc *
    prop_hatchery
  
  
  HOS_local <-
    HOR_esc *
    prop_spawn
  
  
  # ============================================================================
  # 5. TOTAL HATCHERY-ORIGIN NATURAL SPAWNERS
  #
  # Total HOS consists of locally produced HORs that naturally spawn plus
  # external hatchery-origin strays.
  # ============================================================================
  
  HOS <-
    HOS_local +
    total_strays
  
  
  # ============================================================================
  # 6. EFFECTIVE HOS
  #
  # RRS-adjusted HOS is used to represent the effective contribution of
  # hatchery-origin natural spawners in the Ford model.
  # ============================================================================
  
  HOS_effective <-
    HOS *
    rrs_hos
  
  
  # ============================================================================
  # 7. HOS MANAGEMENT GOAL
  #
  # The AHA management pHOS target appears to use census HOS:
  #
  #                 HOS_goal
  #   pHOS_goal = ----------------
  #                NOS + HOS_goal
  #
  # Solving for HOS_goal:
  #
  #                  pHOS_goal * NOS
  #   HOS_goal = -------------------------
  #                   1 - pHOS_goal
  #
  # RRS is NOT applied in this calculation.
  #
  # ============================================================================
  
  if (
    pHOS_goal > 0 &&
    pHOS_goal < 1
  ) {
    
    HOS_goal <-
      (
        pHOS_goal *
          NOS
      ) /
      (
        1 - pHOS_goal
      )
    
  } else if (
    pHOS_goal >= 1
  ) {
    
    HOS_goal <- Inf
    
  } else {
    
    HOS_goal <- 0
  }
  
  
  # ============================================================================
  # 8. HOS ABOVE MANAGEMENT GOAL
  # ============================================================================
  
  HOS_above_goal <-
    max(
      0,
      HOS - HOS_goal
    )
  
  
  # ============================================================================
  # 9. HATCHERY-ORIGIN BROODSTOCK
  #
  # HORs returning to the hatchery fill the remainder of the local
  # broodstock requirement after NOB are selected.
  # ============================================================================
  
  HOB_goal <-
    max(
      0,
      local_brood - NOB
    )
  
  
  HOB <-
    min(
      HOB_goal,
      HOR_hatchery
    )
  
  
  # ============================================================================
  # 10. TOTAL LOCAL BROODSTOCK
  # ============================================================================
  
  brood_total <-
    NOB +
    HOB
  
  
  # ============================================================================
  # 11. HORs REMAINING AT HATCHERY AFTER LOCAL BROODSTOCK
  # ============================================================================
  
  HOR_after_brood <-
    max(
      0,
      HOR_hatchery - HOB
    )
  
  
  # ============================================================================
  # 12. EXPORTED HATCHERY FISH
  #
  # Exported broodstock is removed from hatchery returns remaining after
  # local broodstock collection.
  # ============================================================================
  
  HOR_export <-
    min(
      export_brood_goal,
      HOR_after_brood
    )
  
  
  # ============================================================================
  # 13. HATCHERY SURPLUS AFTER EXPORT
  # ============================================================================
  
  HOR_surplus <-
    max(
      0,
      HOR_after_brood - HOR_export
    )
  
  
  # ============================================================================
  # 14. REALIZED pNOB
  #
  #                NOB
  #   pNOB = ----------------
  #            NOB + HOB
  #
  # ============================================================================
  
  if (brood_total > 0) {
    
    pNOB <-
      NOB /
      brood_total
    
  } else {
    
    pNOB <-
      NA_real_
  }
  
  
  # ============================================================================
  # 15. CENSUS pHOS
  #
  #                   HOS
  #   pHOS = ---------------------
  #              NOS + HOS
  #
  # ============================================================================
  
  if (
    (NOS + HOS) > 0
  ) {
    
    pHOS_census <-
      HOS /
      (
        NOS + HOS
      )
    
  } else {
    
    pHOS_census <-
      NA_real_
  }
  
  
  # ============================================================================
  # 16. EFFECTIVE pHOS
  #
  # AHA/Ford pHOS uses RRS-adjusted HOS:
  #
  #                   RRS * HOS
  #   pHOS = -----------------------------
  #              NOS + RRS * HOS
  #
  # ============================================================================
  
  if (
    (NOS + HOS_effective) > 0
  ) {
    
    pHOS_effective <-
      HOS_effective /
      (
        NOS + HOS_effective
      )
    
  } else {
    
    pHOS_effective <-
      NA_real_
  }
  
  
  # ============================================================================
  # 17. PROPORTIONATE NATURAL INFLUENCE
  #
  #                     pNOB
  #   PNI = ---------------------------
  #              pNOB + pHOS
  #
  # Use effective pHOS because this is the value entering the Ford
  # fitness calculation.
  #
  # ============================================================================
  
  if (
    !is.na(pNOB) &&
    !is.na(pHOS_effective) &&
    (pNOB + pHOS_effective) > 0
  ) {
    
    PNI <-
      pNOB /
      (
        pNOB +
          pHOS_effective
      )
    
  } else {
    
    PNI <-
      NA_real_
  }
  
  
  # ============================================================================
  # 18. TOTAL SPAWNERS
  # ============================================================================
  
  natural_spawners <-
    NOS +
    HOS
  
  
  effective_natural_spawners <-
    NOS +
    HOS_effective
  
  
  # ============================================================================
  # 19. RETURN RESULTS
  # ============================================================================
  
  list(
    
    # Adult escapement
    NOR_esc = NOR_esc,
    HOR_esc = HOR_esc,
    
    # NOR broodstock
    NOB_goal = NOB_goal,
    NOB_max = NOB_max,
    NOB = NOB,
    
    # NOR natural spawning
    NOS = NOS,
    
    # HOR disposition
    HOR_hatchery = HOR_hatchery,
    HOS_local = HOS_local,
    total_strays = total_strays,
    
    # HOS
    HOS_goal = HOS_goal,
    HOS = HOS,
    HOS_effective = HOS_effective,
    HOS_above_goal = HOS_above_goal,
    
    # HOR broodstock
    HOB_goal = HOB_goal,
    HOB = HOB,
    
    # Broodstock totals
    brood_total = brood_total,
    
    # Hatchery disposition
    HOR_after_brood = HOR_after_brood,
    HOR_export = HOR_export,
    HOR_surplus = HOR_surplus,
    
    # Genetic metrics
    pNOB = pNOB,
    pHOS_census = pHOS_census,
    pHOS_effective = pHOS_effective,
    PNI = PNI,
    
    # Spawner totals
    natural_spawners = natural_spawners,
    effective_natural_spawners =
      effective_natural_spawners
  )
}


# ==============================================================================
# AHA ADULT ALLOCATION + TWO-POPULATION FORD MODEL
#
# Purpose:
#   Connect adult demographic allocation to the two-population Ford genetic
#   model.
#
#   Adult allocation determines realized pNOB and effective pHOS. These
#   realized values, rather than the management goals, are passed to Ford.
#
# Returns:
#   Adult allocation results
#   Realized pNOB and pHOS
#   Updated natural and hatchery phenotypes
#   Natural-population fitness
#
# ==============================================================================

aha_ford_generation <- function(
    P_nat,
    P_hat,
    NOR_esc,
    HOR_esc,
    pars
) {
  
  # ============================================================================
  # 1. ADULT ALLOCATION
  # ============================================================================
  
  allocation <- aha_allocate_adults(
    
    NOR_esc = NOR_esc,
    HOR_esc = HOR_esc,
    
    local_brood =
      pars$local_brood,
    
    export_brood_goal =
      pars$export_brood_goal,
    
    pNOB_goal =
      pars$pNOB_goal,
    
    pHOS_goal =
      pars$pHOS_goal,
    
    max_prop_nor_brood =
      pars$max_prop_nor_brood,
    
    prop_hatchery =
      pars$prop_hatchery,
    
    prop_spawn =
      pars$prop_natural_spawn,
    
    total_strays =
      pars$total_strays,
    
    rrs_hos =
      pars$rrs_hos
  )
  
  
  # ============================================================================
  # 2. REALIZED GENETIC INPUTS
  #
  # These come from the demographic allocation, NOT the management goals.
  # ============================================================================
  
  pNOB_realized <-
    allocation$pNOB
  
  pHOS_realized <-
    allocation$pHOS_effective
  
  
  # ============================================================================
  # 3. HANDLE MISSING pNOB / pHOS
  #
  # If there are no broodstock or natural spawners, the corresponding
  # proportion may be undefined.
  #
  # For now, set undefined influence to zero so that the model can continue.
  # We can revisit this behavior if an AHA test scenario indicates otherwise.
  # ============================================================================
  
  if (is.na(pNOB_realized)) {
    pNOB_realized <- 0
  }
  
  if (is.na(pHOS_realized)) {
    pHOS_realized <- 0
  }
  
  
  # ============================================================================
  # 4. TWO-POPULATION FORD UPDATE
  # ============================================================================
  
  ford <-
    ford_generation(
      
      P_nat = P_nat,
      P_hat = P_hat,
      
      pHOS = pHOS_realized,
      pNOB = pNOB_realized,
      
      pars = pars
    )
  
  
  # ============================================================================
  # 5. NATURAL-POPULATION FITNESS
  #
  # Fitness is calculated from the updated natural phenotype.
  # ============================================================================
  
  fitness <-
    aha_fitness(
      P_nat = ford$P_nat,
      pars = pars
    )
  
  
  # ============================================================================
  # 6. RETURN RESULTS
  # ============================================================================
  
  list(
    
    # Adult allocation
    allocation = allocation,
    
    # Realized genetic inputs
    pNOB = pNOB_realized,
    pHOS = pHOS_realized,
    
    # Previous phenotypes
    P_nat_previous = P_nat,
    P_hat_previous = P_hat,
    
    # Updated phenotypes
    P_nat = ford$P_nat,
    P_hat = ford$P_hat,
    
    # Fitness
    fitness = fitness
  )
}

# ==============================================================================
# RUN TWO-POPULATION AHA / FORD MODEL
#
# Iterates adult allocation and Ford genetic change through generations.
#
# NOR_esc and HOR_esc can either be:
#   - single values, held constant across generations, or
#   - vectors containing generation-specific escapements.
#
# ==============================================================================

aha_ford_sim <- function(
    generations,
    P_nat_start,
    P_hat_start,
    NOR_esc,
    HOR_esc,
    pars
) {
  
  # ---------------------------------------------------------------------------
  # Expand constant escapements if necessary
  # ---------------------------------------------------------------------------
  
  if (length(NOR_esc) == 1) {
    NOR_esc <- rep(NOR_esc, generations)
  }
  
  if (length(HOR_esc) == 1) {
    HOR_esc <- rep(HOR_esc, generations)
  }
  
  
  if (length(NOR_esc) < generations) {
    stop("NOR_esc must contain at least 'generations' values.")
  }
  
  if (length(HOR_esc) < generations) {
    stop("HOR_esc must contain at least 'generations' values.")
  }
  
  
  # ---------------------------------------------------------------------------
  # Storage
  # ---------------------------------------------------------------------------
  
  out <- data.frame(
    
    generation = seq_len(generations),
    
    NOR_esc = NA_real_,
    HOR_esc = NA_real_,
    
    NOB = NA_real_,
    HOB = NA_real_,
    
    NOS = NA_real_,
    HOS = NA_real_,
    HOS_effective = NA_real_,
    
    pNOB = NA_real_,
    pHOS_census = NA_real_,
    pHOS_effective = NA_real_,
    PNI = NA_real_,
    
    P_nat_start = NA_real_,
    P_hat_start = NA_real_,
    
    P_nat = NA_real_,
    P_hat = NA_real_,
    
    fitness = NA_real_
  )
  
  
  # ---------------------------------------------------------------------------
  # Initial phenotype values
  # ---------------------------------------------------------------------------
  
  P_nat <- P_nat_start
  P_hat <- P_hat_start
  
  
  # ===========================================================================
  # GENERATION LOOP
  # ===========================================================================
  
  for (g in seq_len(generations)) {
    
    res <- aha_ford_generation(
      
      P_nat = P_nat,
      P_hat = P_hat,
      
      NOR_esc = NOR_esc[g],
      HOR_esc = HOR_esc[g],
      
      pars = pars
    )
    
    
    # -------------------------------------------------------------------------
    # Save results
    # -------------------------------------------------------------------------
    
    out$NOR_esc[g] <-
      NOR_esc[g]
    
    out$HOR_esc[g] <-
      HOR_esc[g]
    
    
    out$NOB[g] <-
      res$allocation$NOB
    
    out$HOB[g] <-
      res$allocation$HOB
    
    
    out$NOS[g] <-
      res$allocation$NOS
    
    out$HOS[g] <-
      res$allocation$HOS
    
    out$HOS_effective[g] <-
      res$allocation$HOS_effective
    
    
    out$pNOB[g] <-
      res$allocation$pNOB
    
    out$pHOS_census[g] <-
      res$allocation$pHOS_census
    
    out$pHOS_effective[g] <-
      res$allocation$pHOS_effective
    
    out$PNI[g] <-
      res$allocation$PNI
    
    
    out$P_nat_start[g] <-
      res$P_nat_previous
    
    out$P_hat_start[g] <-
      res$P_hat_previous
    
    
    out$P_nat[g] <-
      res$P_nat
    
    out$P_hat[g] <-
      res$P_hat
    
    
    out$fitness[g] <-
      res$fitness
    
    
    # -------------------------------------------------------------------------
    # Update phenotype state for next generation
    # -------------------------------------------------------------------------
    
    P_nat <-
      res$P_nat
    
    P_hat <-
      res$P_hat
  }
  
  
  out
}