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

get_recovery_targets <- function(
    recovery_phase,
    
    # Preservation
    preservation_pNOB = NA_real_,
    preservation_pHOS = NA_real_,
    preservation_PNI  = NA_real_,
    
    # Recolonization
    recolonization_pNOB = 0.20,
    recolonization_pHOS = NA_real_,
    recolonization_PNI  = NA_real_,
    
    # Local Adaptation
    local_adaptation_pNOB = 0.50,
    local_adaptation_pHOS = 0.30,
    local_adaptation_PNI  = 0.62,
    
    # Full Restoration
    full_restoration_pNOB = 0.60,
    full_restoration_pHOS = 0.30,
    full_restoration_PNI  = 0.67
) {
  
  switch(
    recovery_phase,
    
    "Preservation" = list(
      pNOB_target = preservation_pNOB,
      pHOS_max = preservation_pHOS,
      PNI_min = preservation_PNI
    ),
    
    "Recolonization" = list(
      pNOB_target = recolonization_pNOB,
      pHOS_max = recolonization_pHOS,
      PNI_min = recolonization_PNI
    ),
    
    "Local Adaptation" = list(
      pNOB_target = local_adaptation_pNOB,
      pHOS_max = local_adaptation_pHOS,
      PNI_min = local_adaptation_PNI
    ),
    
    "Full Restoration" = list(
      pNOB_target = full_restoration_pNOB,
      pHOS_max = full_restoration_pHOS,
      PNI_min = full_restoration_PNI
    ),
    
    stop(
      paste("Unknown recovery phase:", recovery_phase),
      call. = FALSE
    )
  )
}

calculate_recovery_brood <- function(
    brood_need,
    pNOB_target,
    no_captured_weir,
    ho_captured_weir
) {
  
  # Desired number of NO broodstock based on recovery-phase target.
  no_brood_target <- brood_need * pNOB_target
  
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
  brood_total <- no_brood_actual + ho_brood_actual
  
  # Remaining broodstock shortfall, if insufficient NO + HO are available.
  brood_deficit <- max(
    0,
    brood_need - brood_total
  )
  
  # Realized pNOB.
  pNOB <- ifelse(
    brood_total > 0,
    no_brood_actual / brood_total,
    NA_real_
  )
  
  list(
    no_brood_target = no_brood_target,
    no_brood_actual = no_brood_actual,
    ho_brood_needed = ho_brood_needed,
    ho_brood_actual = ho_brood_actual,
    brood_total = brood_total,
    brood_deficit = brood_deficit,
    pNOB_target = pNOB_target,
    pNOB = pNOB,
    
    no_captured_available_above =
      no_captured_weir - no_brood_actual,
    
    ho_captured_available_above =
      ho_captured_weir - ho_brood_actual
  )
}

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
  no_spawners_above <-
    (
      weir_result$no_above_weir_uncaptured +
        brood_result$no_captured_available_above
    ) *
    accounting_params$survival_above_weir
  
  # NO fish remaining below the weir.
  no_spawners_below <-
    weir_result$no_below_weir *
    accounting_params$survival_below_weir
  
  NOS <-
    no_spawners_above +
    no_spawners_below
  
  
  # ------------------------------------------------------------
  # Hatchery-origin spawners that cannot be controlled
  # ------------------------------------------------------------
  
  # HO fish that pass above the weir without being captured.
  ho_spawners_above_uncaptured <-
    weir_result$ho_above_weir_uncaptured *
    accounting_params$survival_above_weir
  
  # HO fish remaining below the weir.
  ho_spawners_below <-
    weir_result$ho_below_weir *
    accounting_params$survival_below_weir
  
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
    pHOS_max * NOS / (1 - pHOS_max)
  }
  
  
  # ------------------------------------------------------------
  # Captured HO fish that may be released
  # ------------------------------------------------------------
  
  # Remaining room under the pHOS limit after accounting for
  # uncontrolled HO spawners.
  HOS_available <-
    max(
      0,
      HOS_target_max - HOS_uncontrolled
    )
  
  # Convert allowable successful spawners back to the number of
  # captured HO fish that can be released above the weir.
  ho_release_needed <- if (
    accounting_params$survival_above_weir > 0
  ) {
    HOS_available /
      accounting_params$survival_above_weir
  } else {
    0
  }
  
  # Cannot release more HO fish than remain after broodstock.
  ho_captured_released <-
    min(
      brood_result$ho_captured_available_above,
      ho_release_needed
    )
  
  # Remaining captured HO fish are removed.
  ho_captured_removed <-
    max(
      0,
      brood_result$ho_captured_available_above -
        ho_captured_released
    )
  
  
  # ------------------------------------------------------------
  # Final HO spawners
  # ------------------------------------------------------------
  
  ho_spawners_above_released <-
    ho_captured_released *
    accounting_params$survival_above_weir
  
  HOS <-
    HOS_uncontrolled +
    ho_spawners_above_released
  
  system_spawners <-
    NOS + HOS
  
  pHOS <- ifelse(
    system_spawners > 0,
    HOS / system_spawners,
    NA_real_
  )
  
  pHOS_exceeded <- ifelse(
    is.na(pHOS_max),
    FALSE,
    pHOS > pHOS_max
  )
  
  pHOS_excess <- ifelse(
    is.na(pHOS_max),
    NA_real_,
    max(0, pHOS - pHOS_max)
  )
  
  
  # ------------------------------------------------------------
  # Return accounting
  # ------------------------------------------------------------
  
  list(
    no_spawners_above = no_spawners_above,
    no_spawners_below = no_spawners_below,
    
    ho_spawners_above_uncaptured =
      ho_spawners_above_uncaptured,
    
    ho_spawners_below =
      ho_spawners_below,
    
    HOS_uncontrolled =
      HOS_uncontrolled,
    
    HOS_target_max = HOS_target_max,
    HOS_available = HOS_available,
    
    ho_captured_released =
      ho_captured_released,
    
    ho_captured_removed =
      ho_captured_removed,
    
    ho_spawners_above_released =
      ho_spawners_above_released,
    
    no_spawners_total = NOS,
    ho_spawners_total = HOS,
    system_spawners_total = system_spawners,
    
    pHOS_max = pHOS_max,
    pHOS = pHOS,
    pHOS_exceeded = pHOS_exceeded,
    pHOS_excess = pHOS_excess
  )
}


calculate_recovery_fisheries <- function(
    scenario_inputs,
    accounting_params,
    pNOB_target
) {
  
  # ------------------------------------------------------------
  # Project NO fishery impacts
  # ------------------------------------------------------------
  
  no_wl <-
    scenario_inputs$wl_scaling *
    scenario_inputs$no_manarea_est
  
  sport_impact_wl_projected <-
    max(0, impact_sport_wl(no_wl))
  
  sport_no_impacts_projected <-
    max(
      0,
      sport_impact_wl_projected /
        scenario_inputs$wl_scaling
    )
  
  treaty_no_impacts_projected <-
    max(
      0,
      impact_treaty_wl_direct(no_wl)
    )
  
  total_no_impacts_projected <-
    sport_no_impacts_projected +
    treaty_no_impacts_projected
  
  
  # ------------------------------------------------------------
  # Protect NO broodstock target
  # ------------------------------------------------------------
  
  # Recovery-phase pNOB determines the desired NO contribution
  # to broodstock.
  no_brood_target <-
    pNOB_target *
    scenario_inputs$brood_need
  
  # Number of NO fish that must remain after fisheries so the
  # desired NO broodstock can be captured at the weir.
  required_no_after_fisheries <- if (
    !is.na(accounting_params$trap_prop_no) &&
    accounting_params$trap_prop_no > 0
  ) {
    no_brood_target /
      accounting_params$trap_prop_no
  } else {
    Inf
  }
  
  # Maximum allowable NO fishery impacts while protecting the
  # desired NO broodstock contribution.
  max_no_impacts_allowed <-
    max(
      0,
      scenario_inputs$no_manarea_est -
        required_no_after_fisheries
    )
  
  no_impacts_allowed <-
    scale_to_cap(
      values = c(
        sport_no_impacts_projected,
        treaty_no_impacts_projected
      ),
      cap = max_no_impacts_allowed
    )
  
  sport_no_impacts <-
    no_impacts_allowed[1]
  
  treaty_no_impacts <-
    no_impacts_allowed[2]
  
  total_no_impacts_allowed <-
    sport_no_impacts +
    treaty_no_impacts
  
  
  # ------------------------------------------------------------
  # Project HO sport harvest
  # ------------------------------------------------------------
  
  # Convert allowed Lostine sport NO impact back to the
  # Wallowa/Lostine scale.
  sport_impact_wl <-
    max(
      0,
      sport_no_impacts *
        scenario_inputs$wl_scaling
    )
  
  # Sport NO impacts represent 10% hooking mortality.
  # Therefore, total NO encounters are 10x the NO impact.
  # Expected HO harvest is proportional to the HO:NO abundance
  # ratio.
  if (
    is.na(scenario_inputs$ho_manarea_est) ||
    scenario_inputs$ho_manarea_est < 20 ||
    is.na(sport_no_impacts) ||
    sport_no_impacts <= 0 ||
    scenario_inputs$no_manarea_est <= 0
  ) {
    
    ho_sport_harvest_projected <- 0
    sport_closed <- TRUE
    
  } else {
    
    allowed_no_handle <-
      10 * sport_no_impacts
    
    ho_sport_harvest_projected <-
      max(
        0,
        round(
          allowed_no_handle *
            (
              scenario_inputs$ho_manarea_est /
                scenario_inputs$no_manarea_est
            )
        )
      )
    
    sport_closed <- FALSE
  }
  
  
  # ------------------------------------------------------------
  # Project HO treaty harvest
  # ------------------------------------------------------------
  
  ho_treaty_harvest_projected <-
    calculate_ho_treaty_harvest_balance(
      sport_impact_wl = sport_impact_wl,
      ho_sport_harvest =
        ho_sport_harvest_projected,
      treaty_no_impacts =
        treaty_no_impacts
    )
  
  total_ho_harvest_projected <-
    ho_sport_harvest_projected +
    ho_treaty_harvest_projected
  
  
  # ------------------------------------------------------------
  # Determine actual NO brood availability after fisheries
  # ------------------------------------------------------------
  
  no_after_fisheries_for_brood_check <-
    scenario_inputs$no_manarea_est -
    sport_no_impacts -
    treaty_no_impacts
  
  no_captured_for_brood_check <-
    accounting_params$trap_prop_no *
    no_after_fisheries_for_brood_check
  
  no_brood_actual_for_brood_check <-
    min(
      no_brood_target,
      no_captured_for_brood_check
    )
  
  
  # ------------------------------------------------------------
  # Protect HO needed to backfill broodstock
  # ------------------------------------------------------------
  
  ho_brood_needed <-
    max(
      0,
      scenario_inputs$brood_need -
        no_brood_actual_for_brood_check
    )
  
  required_ho_after_fisheries <- if (
    !is.na(accounting_params$trap_prop_ho) &&
    accounting_params$trap_prop_ho > 0
  ) {
    ho_brood_needed /
      accounting_params$trap_prop_ho
  } else {
    Inf
  }
  
  max_ho_harvest_allowed <-
    max(
      0,
      scenario_inputs$ho_manarea_est -
        required_ho_after_fisheries
    )
  
  ho_harvest_allowed <-
    scale_to_cap(
      values = c(
        ho_sport_harvest_projected,
        ho_treaty_harvest_projected
      ),
      cap = max_ho_harvest_allowed
    )
  
  ho_sport_harvest <-
    ho_harvest_allowed[1]
  
  ho_treaty_harvest <-
    ho_harvest_allowed[2]
  
  total_ho_harvest_allowed <-
    ho_sport_harvest +
    ho_treaty_harvest
  
  
  # ------------------------------------------------------------
  # Return fishery accounting
  # ------------------------------------------------------------
  
  list(
    no_wl = no_wl,
    
    sport_no_impacts_projected =
      sport_no_impacts_projected,
    
    treaty_no_impacts_projected =
      treaty_no_impacts_projected,
    
    total_no_impacts_projected =
      total_no_impacts_projected,
    
    no_brood_target =
      no_brood_target,
    
    required_no_after_fisheries =
      required_no_after_fisheries,
    
    max_no_impacts_allowed =
      max_no_impacts_allowed,
    
    sport_no_impacts =
      sport_no_impacts,
    
    treaty_no_impacts =
      treaty_no_impacts,
    
    total_no_impacts_allowed =
      total_no_impacts_allowed,
    
    ho_sport_harvest_projected =
      ho_sport_harvest_projected,
    
    ho_treaty_harvest_projected =
      ho_treaty_harvest_projected,
    
    total_ho_harvest_projected =
      total_ho_harvest_projected,
    
    ho_brood_needed =
      ho_brood_needed,
    
    required_ho_after_fisheries =
      required_ho_after_fisheries,
    
    max_ho_harvest_allowed =
      max_ho_harvest_allowed,
    
    ho_sport_harvest =
      ho_sport_harvest,
    
    ho_treaty_harvest =
      ho_treaty_harvest,
    
    total_ho_harvest_allowed =
      total_ho_harvest_allowed,
    
    sport_closed =
      sport_closed
  )
}


calculate_recovery_metrics <- function(
    brood_result,
    spawner_result,
    PNI_min = NA_real_
) {
  
  pNOB <- brood_result$pNOB
  pHOS <- spawner_result$pHOS
  
  PNI <- if (
    !is.na(pNOB) &&
    !is.na(pHOS) &&
    (pNOB + pHOS) > 0
  ) {
    pNOB / (pNOB + pHOS)
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
    pNOB = pNOB,
    pHOS = pHOS,
    PNI = PNI,
    
    PNI_min = PNI_min,
    PNI_met = PNI_met,
    PNI_shortfall = PNI_shortfall
  )
}

run_recovery_scenario <- function(
    scenario_inputs,
    accounting_params,
    recovery_abundance,
    cbp_low,
    cbp_medium,
    cbp_high,
    
    # Preservation targets
    preservation_pNOB = NA_real_,
    preservation_pHOS = NA_real_,
    preservation_PNI  = NA_real_,
    
    # Recolonization targets
    recolonization_pNOB = 0.20,
    recolonization_pHOS = NA_real_,
    recolonization_PNI  = NA_real_,
    
    # Local Adaptation targets
    local_adaptation_pNOB = 0.50,
    local_adaptation_pHOS = 0.30,
    local_adaptation_PNI  = 0.62,
    
    # Full Restoration targets
    full_restoration_pNOB = 0.60,
    full_restoration_pHOS = 0.30,
    full_restoration_PNI  = 0.67
) {
  
  # ------------------------------------------------------------
  # 1. Determine recovery phase
  # ------------------------------------------------------------
  
  recovery_phase <- calculate_recovery_phase(
    abundance = recovery_abundance,
    cbp_low = cbp_low,
    cbp_medium = cbp_medium,
    cbp_high = cbp_high
  )
  
  
  # ------------------------------------------------------------
  # 2. Get management targets for recovery phase
  # ------------------------------------------------------------
  
  recovery_targets <- get_recovery_targets(
    recovery_phase = recovery_phase,
    
    preservation_pNOB = preservation_pNOB,
    preservation_pHOS = preservation_pHOS,
    preservation_PNI = preservation_PNI,
    
    recolonization_pNOB = recolonization_pNOB,
    recolonization_pHOS = recolonization_pHOS,
    recolonization_PNI = recolonization_PNI,
    
    local_adaptation_pNOB = local_adaptation_pNOB,
    local_adaptation_pHOS = local_adaptation_pHOS,
    local_adaptation_PNI = local_adaptation_PNI,
    
    full_restoration_pNOB = full_restoration_pNOB,
    full_restoration_pHOS = full_restoration_pHOS,
    full_restoration_PNI = full_restoration_PNI
  )
  
  
  # ------------------------------------------------------------
  # 3. Calculate fisheries
  # ------------------------------------------------------------
  
  fishery_result <- calculate_recovery_fisheries(
    scenario_inputs = scenario_inputs,
    accounting_params = accounting_params,
    pNOB_target = recovery_targets$pNOB_target
  )
  
  
  # ------------------------------------------------------------
  # 4. Calculate post-fishery weir accounting
  # ------------------------------------------------------------
  
  weir_result <- calculate_scenario_weir_accounting(
    scenario_inputs = scenario_inputs,
    accounting_params = accounting_params,
    fishery_result = fishery_result
  )
  
  
  # ------------------------------------------------------------
  # 5. Allocate broodstock
  # ------------------------------------------------------------
  
  brood_result <- calculate_recovery_brood(
    brood_need = scenario_inputs$brood_need,
    pNOB_target = recovery_targets$pNOB_target,
    no_captured_weir = weir_result$no_captured_weir,
    ho_captured_weir = weir_result$ho_captured_weir
  )
  
  
  # ------------------------------------------------------------
  # 6. Allocate natural spawners
  # ------------------------------------------------------------
  
  spawner_result <- calculate_recovery_spawners(
    accounting_params = accounting_params,
    weir_result = weir_result,
    brood_result = brood_result,
    pHOS_max = recovery_targets$pHOS_max
  )
  
  
  # ------------------------------------------------------------
  # 7. Calculate realized management metrics
  # ------------------------------------------------------------
  
  metrics <- calculate_recovery_metrics(
    brood_result = brood_result,
    spawner_result = spawner_result,
    PNI_min = recovery_targets$PNI_min
  )
  
  
  # ------------------------------------------------------------
  # Return complete accounting
  # ------------------------------------------------------------
  
  list(
    recovery_abundance = recovery_abundance,
    recovery_phase = recovery_phase,
    recovery_targets = recovery_targets,
    fishery_result = fishery_result,
    weir_result = weir_result,
    brood_result = brood_result,
    spawner_result = spawner_result,
    metrics = metrics
  )
}

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
      result$fishery_result$total_no_impacts_allowed,
    
    HO_sport_harvest =
      result$fishery_result$ho_sport_harvest,
    
    HO_treaty_harvest =
      result$fishery_result$ho_treaty_harvest,
    
    total_HO_harvest =
      result$fishery_result$total_ho_harvest_allowed,
    
    
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


#-------------------------------------------
# validation
#--------------------------------------------

source("./R/helpers.R")
source("./R/harvest_funs.R")
source("./R/weir_funs.R")
source("./R/recovery_funs.R")

# Existing Lostine parameter block

test_inputs <- create_scenario_inputs(
  no_manarea_est = 1000,
  ho_manarea_est = 1000,
  brood_need = 160,
  spawner_goal = 800,
  wl_scaling = 1.4,
  utilization_no = 1,
  utilization_ho = 1
)

test_result <- run_recovery_scenario(
  scenario_inputs = test_inputs,
  accounting_params = accounting_params,
  
  # Force the test into Local Adaptation
  recovery_abundance = 600,
  cbp_low = 200,
  cbp_medium = 500,
  cbp_high = 1000
)

tidy_recovery_result(test_result)

