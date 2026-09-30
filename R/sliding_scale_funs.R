# Purpose: Functions implementing the historical Lostine River
# spring Chinook Sliding Scale management framework.
#
# Author: Ryan N. Kinzer
# Created: 09/15/2026
#
# Reference:
# HGMP document.
#
# The Sliding Scale uses natural-origin adult abundance to determine:
#   1. target proportion of natural-origin broodstock (pNOB)
#   2. proportion of hatchery-origin adults released above the weir
#
# Fisheries are calculated independently using the common
# calculate_fisheries() function. Sliding Scale management begins
# after fisheries and determines broodstock composition and
# hatchery-origin disposition at the weir.


#------------------------------------------------------------
# Sliding Scale targets
#------------------------------------------------------------

get_sliding_scale_targets <- function(
    no_abundance
) {
  
  if (is.na(no_abundance)) {
    stop(
      "NO abundance cannot be NA.",
      call. = FALSE
    )
  }
  
  
  # ------------------------------------------------------------
  # Historical Sliding Scale
  # ------------------------------------------------------------
  
  if (no_abundance < 8) {
    
    pNOB_target <- NA_real_
    ho_release_prop <- NA_real_
    
  } else if (no_abundance < 75) {
    
    pNOB_target <- NA_real_
    ho_release_prop <- NA_real_
    
  } else if (no_abundance < 150) {
    
    pNOB_target <- 0.20
    ho_release_prop <- 0.70
    
  } else if (no_abundance < 250) {
    
    pNOB_target <- 0.25
    ho_release_prop <- 0.60
    
  } else if (no_abundance < 500) {
    
    pNOB_target <- 0.30
    ho_release_prop <- 0.50
    
  } else if (no_abundance < 750) {
    
    pNOB_target <- 0.40
    ho_release_prop <- 0.40
    
  } else if (no_abundance < 1000) {
    
    pNOB_target <- 0.50
    ho_release_prop <- 0.25
    
  } else {
    
    pNOB_target <- 1.00
    
    # HGMP specifies <10% HO release for the highest abundance
    # category. A value of 0.10 is used here as the operational
    # approximation for deterministic simulation.
    #
    # The HGMP identifies the preceding abundance category as
    # 750-999 and the highest category as >1000. Because the
    # table does not explicitly assign exactly 1000 adults,
    # the model treats 1000 as the start of the highest category.
    ho_release_prop <- 0.10
  }
  
  
  list(
    pNOB_target =
      pNOB_target,
    
    ho_release_prop =
      ho_release_prop
  )
}


#------------------------------------------------------------
# Broodstock allocation
#------------------------------------------------------------

calculate_sliding_scale_brood <- function(
    weir_result,
    scenario_inputs,
    sliding_targets
) {
  
  pNOB_target <-
    sliding_targets$pNOB_target
  
  
  # ------------------------------------------------------------
  # Desired NO broodstock
  # ------------------------------------------------------------
  
  # Desired NO brood based on Sliding Scale pNOB.
  #
  # NOTE:
  # The first two historical Sliding Scale rows do not specify
  # a pNOB target. The existing model convention is retained here
  # and assigns a target of zero NO broodstock when pNOB is NA.
  
  no_brood_target <- if (
    is.na(pNOB_target)
  ) {
    
    0
    
  } else {
    
    round_fish(
      scenario_inputs$brood_need *
        pNOB_target
    )
  }
  
  
  # ------------------------------------------------------------
  # Actual broodstock allocation
  # ------------------------------------------------------------
  
  # Actual NO brood cannot exceed NO captured at the weir.
  no_brood <- min(
    no_brood_target,
    weir_result$no_captured_weir
  )
  
  # Fill remaining brood need with captured HO adults.
  ho_brood_needed <- max(
    0,
    scenario_inputs$brood_need -
      no_brood
  )
  
  ho_brood <- min(
    ho_brood_needed,
    weir_result$ho_captured_weir
  )
  
  total_brood <-
    no_brood +
    ho_brood
  
  brood_deficit <- max(
    0,
    scenario_inputs$brood_need -
      total_brood
  )
  
  
  # ------------------------------------------------------------
  # Adults remaining after broodstock collection
  # ------------------------------------------------------------
  
  no_captured_available <- max(
    0,
    weir_result$no_captured_weir -
      no_brood
  )
  
  ho_captured_available <- max(
    0,
    weir_result$ho_captured_weir -
      ho_brood
  )
  
  
  # ------------------------------------------------------------
  # Return broodstock accounting
  # ------------------------------------------------------------
  
  list(
    no_brood_target =
      no_brood_target,
    
    no_brood =
      no_brood,
    
    ho_brood =
      ho_brood,
    
    total_brood =
      total_brood,
    
    brood_deficit =
      brood_deficit,
    
    no_captured_available =
      no_captured_available,
    
    ho_captured_available =
      ho_captured_available
  )
}


#------------------------------------------------------------
# Spawner allocation
#------------------------------------------------------------

calculate_sliding_scale_spawners <- function(
    weir_result,
    brood_result,
    accounting_params,
    sliding_targets
) {
  
  # ------------------------------------------------------------
  # Natural-origin adults
  # ------------------------------------------------------------
  
  # All remaining captured NO adults are released above the weir.
  no_captured_released <-
    brood_result$no_captured_available
  
  # NO spawning above the weir includes captured/released fish plus
  # fish estimated to have passed upstream without being captured.
  no_spawners_above <- round_fish(
    (
      no_captured_released +
        weir_result$no_above_weir_uncaptured
    ) *
      accounting_params$survival_above_weir
  )
  
  # NO adults remaining below the weir.
  no_spawners_below <- round_fish(
    weir_result$no_below_weir *
      accounting_params$survival_below_weir
  )
  
  no_spawners <-
    no_spawners_above +
    no_spawners_below
  
  
  # ------------------------------------------------------------
  # Uncontrolled hatchery-origin adults
  # ------------------------------------------------------------
  
  # HOR that pass upstream without being captured cannot be
  # selectively removed and therefore contribute to natural spawning.
  ho_spawners_uncontrolled_above <- round_fish(
    weir_result$ho_above_weir_uncaptured *
      accounting_params$survival_above_weir
  )
  
  # HOR remaining below the weir also contribute to natural spawning.
  ho_spawners_below <- round_fish(
    weir_result$ho_below_weir *
      accounting_params$survival_below_weir
  )
  
  
  # ------------------------------------------------------------
  # Controlled hatchery-origin adults
  # ------------------------------------------------------------
  
  # Historical Sliding Scale specifies the proportion of captured
  # hatchery adults remaining after broodstock collection that are
  # released above the weir.
  
  ho_release_prop <-
    sliding_targets$ho_release_prop
  
  if (is.na(ho_release_prop)) {
    
    ho_captured_released <-
      NA_real_
    
    ho_removed <-
      NA_real_
    
  } else {
    
    # Convert the historical release proportion into a whole-fish
    # management disposition.
    ho_captured_released <- round_fish(
      brood_result$ho_captured_available *
        ho_release_prop
    )
    
    # Removal is the remainder so that released + removed always
    # equals the number of captured HO adults available.
    ho_removed <-
      brood_result$ho_captured_available -
      ho_captured_released
  }
  
  
  # ------------------------------------------------------------
  # Hatchery-origin natural spawners
  # ------------------------------------------------------------
  
  # Apply pre-spawn survival to captured HOR released above the weir.
  ho_spawners_released <- if (
    is.na(ho_captured_released)
  ) {
    
    NA_real_
    
  } else {
    
    round_fish(
      ho_captured_released *
        accounting_params$survival_above_weir
    )
  }
  
  # Total HOR natural spawners.
  ho_spawners <- if (
    is.na(ho_spawners_released)
  ) {
    
    NA_real_
    
  } else {
    
    ho_spawners_uncontrolled_above +
      ho_spawners_below +
      ho_spawners_released
  }
  
  
  # ------------------------------------------------------------
  # Total natural spawning population
  # ------------------------------------------------------------
  
  system_spawners <-
    no_spawners +
    ho_spawners
  
  # Realized proportion of hatchery-origin natural spawners.
  pHOS <- if (
    length(system_spawners) == 1 &&
    !is.na(system_spawners) &&
    system_spawners > 0
  ) {
    
    ho_spawners /
      system_spawners
    
  } else {
    
    NA_real_
  }
  
  
  # ------------------------------------------------------------
  # Return spawner accounting
  # ------------------------------------------------------------
  
  list(
    no_captured_released =
      no_captured_released,
    
    ho_captured_released =
      ho_captured_released,
    
    ho_removed =
      ho_removed,
    
    no_spawners_above =
      no_spawners_above,
    
    no_spawners_below =
      no_spawners_below,
    
    no_spawners =
      no_spawners,
    
    ho_spawners_uncontrolled_above =
      ho_spawners_uncontrolled_above,
    
    ho_spawners_below =
      ho_spawners_below,
    
    ho_spawners_released =
      ho_spawners_released,
    
    ho_spawners =
      ho_spawners,
    
    system_spawners =
      system_spawners,
    
    pHOS =
      pHOS
  )
}


#------------------------------------------------------------
# Performance metrics
#------------------------------------------------------------

calculate_sliding_scale_metrics <- function(
    brood_result,
    spawner_result
) {
  
  pNOB <- if (
    brood_result$total_brood > 0
  ) {
    
    brood_result$no_brood /
      brood_result$total_brood
    
  } else {
    
    NA_real_
  }
  
  pHOS <-
    spawner_result$pHOS
  
  PNI <- if (
    !is.na(pNOB) &&
    !is.na(pHOS) &&
    (pNOB + pHOS) > 0
  ) {
    
    pNOB /
      (pNOB + pHOS)
    
  } else {
    
    NA_real_
  }
  
  list(
    pNOB =
      pNOB,
    
    pHOS =
      pHOS,
    
    PNI =
      PNI
  )
}


#------------------------------------------------------------
# Complete Sliding Scale scenario
#------------------------------------------------------------

run_sliding_scale_scenario <- function(
    scenario_inputs,
    accounting_params
) {
  
  # ------------------------------------------------------------
  # 1. Determine Sliding Scale targets
  # ------------------------------------------------------------
  
  sliding_targets <- get_sliding_scale_targets(
    scenario_inputs$no_manarea_est
  )
  
  
  # ------------------------------------------------------------
  # 2. Calculate common fisheries
  # ------------------------------------------------------------
  
  # Fisheries are independent of the management framework.
  # The same function and fishery rules are used for Anchor,
  # Recovery Phase, and Sliding Scale scenarios.
  
  fishery_result <- calculate_fisheries(
    scenario_inputs = scenario_inputs
  )
  
  
  # ------------------------------------------------------------
  # 3. Calculate post-fishery weir accounting
  # ------------------------------------------------------------
  
  weir_result <- calculate_scenario_weir_accounting(
    scenario_inputs = scenario_inputs,
    accounting_params = accounting_params,
    fishery_result = fishery_result
  )
  
  
  # ------------------------------------------------------------
  # 4. Allocate broodstock
  # ------------------------------------------------------------
  
  # Sliding Scale pNOB is applied here, after fisheries.
  
  brood_result <- calculate_sliding_scale_brood(
    weir_result = weir_result,
    scenario_inputs = scenario_inputs,
    sliding_targets = sliding_targets
  )
  
  
  # ------------------------------------------------------------
  # 5. Allocate natural spawners and HOR disposition
  # ------------------------------------------------------------
  
  # Historical HOR release proportion is applied after broodstock
  # collection.
  
  spawner_result <- calculate_sliding_scale_spawners(
    weir_result = weir_result,
    brood_result = brood_result,
    accounting_params = accounting_params,
    sliding_targets = sliding_targets
  )
  
  
  # ------------------------------------------------------------
  # 6. Calculate realized pNOB, pHOS, and PNI
  # ------------------------------------------------------------
  
  metrics <- calculate_sliding_scale_metrics(
    brood_result = brood_result,
    spawner_result = spawner_result
  )
  
  
  # ------------------------------------------------------------
  # Return complete accounting
  # ------------------------------------------------------------
  
  list(
    sliding_targets =
      sliding_targets,
    
    fishery_result =
      fishery_result,
    
    weir_result =
      weir_result,
    
    brood_result =
      brood_result,
    
    spawner_result =
      spawner_result,
    
    metrics =
      metrics
  )
}


#------------------------------------------------------------
# Tidy output
#------------------------------------------------------------

tidy_sliding_scale_result <- function(result) {
  
  tibble(
    
    # Sliding Scale targets
    pNOB_target =
      result$sliding_targets$pNOB_target,
    
    ho_release_prop =
      result$sliding_targets$ho_release_prop,
    
    
    # Fisheries
    sport_NO_impact =
      result$fishery_result$sport_no_impacts,
    
    treaty_NO_impact =
      result$fishery_result$treaty_no_impacts,
    
    total_NO_impact =
      result$fishery_result$total_no_impacts,
    
    HO_sport_harvest =
      result$fishery_result$ho_sport_harvest,
    
    HO_treaty_harvest =
      result$fishery_result$ho_treaty_harvest,
    
    total_HO_harvest =
      result$fishery_result$total_ho_harvest,
    
    
    # Broodstock
    NO_brood =
      result$brood_result$no_brood,
    
    HO_brood =
      result$brood_result$ho_brood,
    
    total_brood =
      result$brood_result$total_brood,
    
    brood_deficit =
      result$brood_result$brood_deficit,
    
    
    # Natural spawning
    HO_removed =
      result$spawner_result$ho_removed,
    
    NO_spawners =
      result$spawner_result$no_spawners,
    
    HO_spawners =
      result$spawner_result$ho_spawners,
    
    system_spawners =
      result$spawner_result$system_spawners,
    
    
    # Realized hatchery influence
    pHOS =
      result$metrics$pHOS,
    
    pNOB =
      result$metrics$pNOB,
    
    PNI =
      result$metrics$PNI
  )
}