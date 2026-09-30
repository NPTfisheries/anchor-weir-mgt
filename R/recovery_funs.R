# Purpose: Functions for evaluating the Recovery Phase management
# framework for Lostine River spring Chinook.
#
# Author: Ryan N. Kinzer
#
# Recovery-phase management targets are supplied externally through
# a named list. No default management targets are specified here.
#
# Fisheries are calculated independently using the common
# calculate_fisheries() function. Recovery Phase management begins
# after fisheries and determines broodstock composition and
# hatchery-origin disposition at the weir.


#------------------------------------------------------------
# Determine recovery phase
#------------------------------------------------------------

calculate_recovery_phase <- function(
    abundance,
    cbp_low,
    cbp_medium,
    cbp_high
) {
  
  if (is.na(abundance)) {
    return(NA_character_)
  }
  
  if (abundance < cbp_low) {
    return("Preservation")
  }
  
  if (abundance < cbp_medium) {
    return("Recolonization")
  }
  
  if (abundance < cbp_high) {
    return("Local Adaptation")
  }
  
  "Full Restoration"
}


#------------------------------------------------------------
# Get recovery-phase management targets
#------------------------------------------------------------

get_recovery_targets <- function(
    recovery_phase,
    recovery_targets
) {
  
  if (missing(recovery_targets) || is.null(recovery_targets)) {
    stop(
      "`recovery_targets` must be supplied.",
      call. = FALSE
    )
  }
  
  required_phases <- c(
    "Preservation",
    "Recolonization",
    "Local Adaptation",
    "Full Restoration"
  )
  
  missing_phases <- setdiff(
    required_phases,
    names(recovery_targets)
  )
  
  if (length(missing_phases) > 0) {
    stop(
      paste0(
        "`recovery_targets` is missing phase(s): ",
        paste(missing_phases, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  if (
    is.na(recovery_phase) ||
    !recovery_phase %in% required_phases
  ) {
    stop(
      paste0(
        "Unknown recovery phase: ",
        recovery_phase
      ),
      call. = FALSE
    )
  }
  
  phase_targets <-
    recovery_targets[[recovery_phase]]
  
  required_targets <- c(
    "pNOB",
    "pHOS",
    "PNI"
  )
  
  missing_targets <- setdiff(
    required_targets,
    names(phase_targets)
  )
  
  if (length(missing_targets) > 0) {
    stop(
      paste0(
        "Recovery targets for '",
        recovery_phase,
        "' are missing: ",
        paste(missing_targets, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  list(
    pNOB_target = phase_targets$pNOB,
    pHOS_max = phase_targets$pHOS,
    PNI_min = phase_targets$PNI
  )
}


#------------------------------------------------------------
# Allocate recovery-phase broodstock
#------------------------------------------------------------

calculate_recovery_brood <- function(
    brood_need,
    pNOB_target,
    no_captured_weir,
    ho_captured_weir
) {
  
  # Desired number of NO broodstock based on recovery-phase target.
  # pNOB remains continuous, but broodstock must be whole fish.
  no_brood_target <- round_fish(
    brood_need * pNOB_target
  )
  
  # Actual NO broodstock cannot exceed available captured NO fish.
  no_brood_actual <- min(
    no_brood_target,
    no_captured_weir
  )
  
  # Any remaining broodstock need is backfilled with HO fish.
  ho_brood_needed <- max(
    0,
    brood_need - no_brood_actual
  )
  
  # Actual HO broodstock cannot exceed available captured HO fish.
  ho_brood_actual <- min(
    ho_brood_needed,
    ho_captured_weir
  )
  
  # Total broodstock collected.
  brood_total <-
    no_brood_actual +
    ho_brood_actual
  
  # Remaining broodstock shortfall.
  brood_deficit <- max(
    0,
    brood_need - brood_total
  )
  
  # Realized pNOB based on actual broodstock composition.
  pNOB <- ifelse(
    brood_total > 0,
    no_brood_actual / brood_total,
    NA_real_
  )
  
  list(
    no_brood_target =
      no_brood_target,
    
    no_brood_actual =
      no_brood_actual,
    
    ho_brood_needed =
      ho_brood_needed,
    
    ho_brood_actual =
      ho_brood_actual,
    
    brood_total =
      brood_total,
    
    brood_deficit =
      brood_deficit,
    
    pNOB_target =
      pNOB_target,
    
    pNOB =
      pNOB,
    
    no_captured_available_above =
      no_captured_weir -
      no_brood_actual,
    
    ho_captured_available_above =
      ho_captured_weir -
      ho_brood_actual
  )
}


#------------------------------------------------------------
# Allocate recovery-phase natural spawners
#------------------------------------------------------------

calculate_recovery_spawners <- function(
    accounting_params,
    weir_result,
    brood_result,
    pHOS_max
) {
  
  # ------------------------------------------------------------
  # Natural-origin spawners
  # ------------------------------------------------------------
  
  # Captured NO fish not used for broodstock are released above
  # the weir. These fish and uncaptured NO fish above the weir
  # are subject to survival_above_weir.
  no_spawners_above <- round_fish(
    (
      weir_result$no_above_weir_uncaptured +
        brood_result$no_captured_available_above
    ) *
      accounting_params$survival_above_weir
  )
  
  # NO fish remaining below the weir.
  no_spawners_below <- round_fish(
    weir_result$no_below_weir *
      accounting_params$survival_below_weir
  )
  
  NOS <-
    no_spawners_above +
    no_spawners_below
  
  
  # ------------------------------------------------------------
  # Hatchery-origin spawners that cannot be controlled
  # ------------------------------------------------------------
  
  # HO fish that pass above the weir without being captured.
  ho_spawners_above_uncaptured <- round_fish(
    weir_result$ho_above_weir_uncaptured *
      accounting_params$survival_above_weir
  )
  
  # HO fish remaining below the weir.
  ho_spawners_below <- round_fish(
    weir_result$ho_below_weir *
      accounting_params$survival_below_weir
  )
  
  HOS_uncontrolled <-
    ho_spawners_above_uncaptured +
    ho_spawners_below
  
  
  # ------------------------------------------------------------
  # Maximum HOS allowed by pHOS target
  # ------------------------------------------------------------
  
  # pHOS = HOS / (NOS + HOS)
  #
  # Solving for HOS:
  #
  # HOS = pHOS * NOS / (1 - pHOS)
  #
  # This value remains continuous because it defines the
  # theoretical pHOS boundary rather than an actual fish count.
  HOS_target_max <- if (
    is.na(pHOS_max)
  ) {
    
    Inf
    
  } else if (
    pHOS_max <= 0
  ) {
    
    0
    
  } else if (
    pHOS_max >= 1
  ) {
    
    Inf
    
  } else {
    
    pHOS_max *
      NOS /
      (1 - pHOS_max)
  }
  
  
  # ------------------------------------------------------------
  # Captured HO fish that may be released
  # ------------------------------------------------------------
  
  # Remaining room under the pHOS limit after accounting for
  # uncontrolled HO spawners.
  HOS_available <- max(
    0,
    HOS_target_max -
      HOS_uncontrolled
  )
  
  # Convert allowable successful spawners back to the number of
  # captured HO adults that can be released above the weir.
  #
  # This is a management disposition and therefore must be a
  # whole number of fish.
  ho_release_needed <- if (
    accounting_params$survival_above_weir > 0
  ) {
    
    if (is.infinite(HOS_available)) {
      
      Inf
      
    } else {
      
      round_fish(
        HOS_available /
          accounting_params$survival_above_weir
      )
    }
    
  } else {
    
    0
  }
  
  # Cannot release more HO fish than remain after broodstock.
  ho_captured_released <- min(
    brood_result$ho_captured_available_above,
    ho_release_needed
  )
  
  # Remaining captured HO fish are removed.
  ho_captured_removed <- max(
    0,
    brood_result$ho_captured_available_above -
      ho_captured_released
  )
  
  
  # ------------------------------------------------------------
  # Final HO spawners
  # ------------------------------------------------------------
  
  # Released HO adults are subject to above-weir survival.
  ho_spawners_above_released <- round_fish(
    ho_captured_released *
      accounting_params$survival_above_weir
  )
  
  HOS <-
    HOS_uncontrolled +
    ho_spawners_above_released
  
  system_spawners <-
    NOS +
    HOS
  
  pHOS <- ifelse(
    system_spawners > 0,
    HOS / system_spawners,
    NA_real_
  )
  
  
  # ------------------------------------------------------------
  # pHOS diagnostics
  # ------------------------------------------------------------
  
  pHOS_exceeded <- ifelse(
    is.na(pHOS_max),
    FALSE,
    pHOS > pHOS_max
  )
  
  pHOS_excess <- ifelse(
    is.na(pHOS_max),
    NA_real_,
    max(
      0,
      pHOS - pHOS_max
    )
  )
  
  
  # ------------------------------------------------------------
  # Return accounting
  # ------------------------------------------------------------
  
  list(
    no_spawners_above =
      no_spawners_above,
    
    no_spawners_below =
      no_spawners_below,
    
    ho_spawners_above_uncaptured =
      ho_spawners_above_uncaptured,
    
    ho_spawners_below =
      ho_spawners_below,
    
    HOS_uncontrolled =
      HOS_uncontrolled,
    
    HOS_target_max =
      HOS_target_max,
    
    HOS_available =
      HOS_available,
    
    ho_captured_released =
      ho_captured_released,
    
    ho_captured_removed =
      ho_captured_removed,
    
    ho_spawners_above_released =
      ho_spawners_above_released,
    
    no_spawners_total =
      NOS,
    
    ho_spawners_total =
      HOS,
    
    system_spawners_total =
      system_spawners,
    
    pHOS_max =
      pHOS_max,
    
    pHOS =
      pHOS,
    
    pHOS_exceeded =
      pHOS_exceeded,
    
    pHOS_excess =
      pHOS_excess
  )
}


#------------------------------------------------------------
# Calculate realized management metrics
#------------------------------------------------------------

calculate_recovery_metrics <- function(
    brood_result,
    spawner_result,
    PNI_min = NA_real_
) {
  
  pNOB <-
    brood_result$pNOB
  
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
  
  PNI_met <- if (
    is.na(PNI_min) ||
    is.na(PNI)
  ) {
    
    NA
    
  } else {
    
    PNI >= PNI_min
  }
  
  PNI_shortfall <- if (
    is.na(PNI_min) ||
    is.na(PNI)
  ) {
    
    NA_real_
    
  } else {
    
    max(
      0,
      PNI_min - PNI
    )
  }
  
  list(
    pNOB =
      pNOB,
    
    pHOS =
      pHOS,
    
    PNI =
      PNI,
    
    PNI_min =
      PNI_min,
    
    PNI_met =
      PNI_met,
    
    PNI_shortfall =
      PNI_shortfall
  )
}


#------------------------------------------------------------
# Run complete Recovery Phase scenario
#------------------------------------------------------------

run_recovery_scenario <- function(
    scenario_inputs,
    accounting_params,
    recovery_abundance,
    cbp_low,
    cbp_medium,
    cbp_high,
    recovery_targets
) {
  
  # ------------------------------------------------------------
  # 1. Validate required recovery targets
  # ------------------------------------------------------------
  
  if (missing(recovery_targets) || is.null(recovery_targets)) {
    stop(
      "`recovery_targets` must be supplied.",
      call. = FALSE
    )
  }
  
  
  # ------------------------------------------------------------
  # 2. Determine recovery phase
  # ------------------------------------------------------------
  
  recovery_phase <- calculate_recovery_phase(
    abundance = recovery_abundance,
    cbp_low = cbp_low,
    cbp_medium = cbp_medium,
    cbp_high = cbp_high
  )
  
  
  # ------------------------------------------------------------
  # 3. Get management targets for recovery phase
  # ------------------------------------------------------------
  
  recovery_targets_phase <- get_recovery_targets(
    recovery_phase = recovery_phase,
    recovery_targets = recovery_targets
  )
  
  
  # ------------------------------------------------------------
  # 4. Calculate common fisheries
  # ------------------------------------------------------------
  
  # Fisheries are independent of the management framework.
  # The same function and fishery rules are used for Anchor,
  # Recovery Phase, and Sliding Scale scenarios.
  
  fishery_result <- calculate_fisheries(
    scenario_inputs = scenario_inputs
  )
  
  
  # ------------------------------------------------------------
  # 5. Calculate post-fishery weir accounting
  # ------------------------------------------------------------
  
  weir_result <- calculate_scenario_weir_accounting(
    scenario_inputs = scenario_inputs,
    accounting_params = accounting_params,
    fishery_result = fishery_result
  )
  
  
  # ------------------------------------------------------------
  # 6. Allocate broodstock
  # ------------------------------------------------------------
  
  # Recovery-phase pNOB is applied here, after fisheries.
  # Fisheries are not reduced to achieve the pNOB target.
  
  brood_result <- calculate_recovery_brood(
    brood_need = scenario_inputs$brood_need,
    pNOB_target = recovery_targets_phase$pNOB_target,
    no_captured_weir = weir_result$no_captured_weir,
    ho_captured_weir = weir_result$ho_captured_weir
  )
  
  
  # ------------------------------------------------------------
  # 7. Allocate natural spawners
  # ------------------------------------------------------------
  
  # Recovery-phase pHOS is applied to disposition of captured
  # hatchery-origin adults remaining after broodstock collection.
  
  spawner_result <- calculate_recovery_spawners(
    accounting_params = accounting_params,
    weir_result = weir_result,
    brood_result = brood_result,
    pHOS_max = recovery_targets_phase$pHOS_max
  )
  
  
  # ------------------------------------------------------------
  # 8. Calculate realized management metrics
  # ------------------------------------------------------------
  
  metrics <- calculate_recovery_metrics(
    brood_result = brood_result,
    spawner_result = spawner_result,
    PNI_min = recovery_targets_phase$PNI_min
  )
  
  
  # ------------------------------------------------------------
  # Return complete accounting
  # ------------------------------------------------------------
  
  list(
    recovery_abundance =
      recovery_abundance,
    
    recovery_phase =
      recovery_phase,
    
    recovery_targets =
      recovery_targets_phase,
    
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
# Tidy Recovery Phase output
#------------------------------------------------------------

tidy_recovery_result <- function(result) {
  
  tibble::tibble(
    
    # Recovery framework
    recovery_abundance =
      result$recovery_abundance,
    
    recovery_phase =
      result$recovery_phase,
    
    pNOB_target =
      result$recovery_targets$pNOB_target,
    
    pHOS_max =
      result$recovery_targets$pHOS_max,
    
    PNI_min =
      result$recovery_targets$PNI_min,
    
    
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
      result$brood_result$no_brood_actual,
    
    HO_brood =
      result$brood_result$ho_brood_actual,
    
    total_brood =
      result$brood_result$brood_total,
    
    brood_deficit =
      result$brood_result$brood_deficit,
    
    
    # Natural spawning
    HO_removed =
      result$spawner_result$ho_captured_removed,
    
    NO_spawners =
      result$spawner_result$no_spawners_total,
    
    HO_spawners =
      result$spawner_result$ho_spawners_total,
    
    system_spawners =
      result$spawner_result$system_spawners_total,
    
    
    # Realized hatchery influence
    pHOS =
      result$metrics$pHOS,
    
    pNOB =
      result$metrics$pNOB,
    
    PNI =
      result$metrics$PNI,
    
    
    # Target diagnostics
    pHOS_exceeded =
      result$spawner_result$pHOS_exceeded,
    
    pHOS_excess =
      result$spawner_result$pHOS_excess,
    
    PNI_met =
      result$metrics$PNI_met,
    
    PNI_shortfall =
      result$metrics$PNI_shortfall
  )
}