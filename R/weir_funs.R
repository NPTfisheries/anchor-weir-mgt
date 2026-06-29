# Purpose: The functions within this script develop weir management numbers.
# Original Author: Kyle Bratcher
# Modified: Ryan Kinzer


# Weir efficiency identity:
#   weir_efficiency = captured / (captured + uncaptured_above)
# Rearranged:
#   uncaptured_above = (captured / weir_efficiency) - captured
calculate_uncaptured_above <- function(captured, weir_efficiency) {
  if (is.na(captured) || is.na(weir_efficiency) || weir_efficiency <= 0) return(NA_real_)
  (captured / weir_efficiency) - captured
}

calculate_scenario_weir_accounting <- function(scenario_inputs, accounting_params, fishery_result) {
  # Utilization rates allow scenarios where fisheries do not fully take the
  # allowed impact or harvest. At 1.0 (default), full utilization is assumed.
  # Reducing these values leaves more fish in the system post-fishery, which
  # affects weir partitioning, broodstock availability, and spawner estimates.
  # Utilization is applied after broodstock-protection caps, so brood
  # availability is always calculated at full allowed impact regardless of
  # whether the fishery is actually expected to be fully utilized.
  no_actual_impacts <-
    (fishery_result$sport_no_impacts + fishery_result$treaty_no_impacts) *
    scenario_inputs$utilization_no
  
  ho_actual_harvest <-
    (fishery_result$ho_sport_harvest + fishery_result$ho_treaty_harvest) *
    scenario_inputs$utilization_ho
  
  no_after_fisheries <- scenario_inputs$no_manarea_est - no_actual_impacts
  ho_after_fisheries <- scenario_inputs$ho_manarea_est - ho_actual_harvest
  
  no_captured_weir <- accounting_params$trap_prop_no * no_after_fisheries
  ho_captured_weir <- accounting_params$trap_prop_ho * ho_after_fisheries
  
  no_above_weir_uncaptured <- calculate_uncaptured_above(no_captured_weir, accounting_params$weir_efficiency)
  ho_above_weir_uncaptured <- calculate_uncaptured_above(ho_captured_weir, accounting_params$weir_efficiency)
  
  no_below_weir <- no_after_fisheries - (no_captured_weir + no_above_weir_uncaptured)
  ho_below_weir <- ho_after_fisheries - (ho_captured_weir + ho_above_weir_uncaptured)
  
  list(
    no_after_fisheries = no_after_fisheries,
    ho_after_fisheries = ho_after_fisheries,
    no_captured_weir = no_captured_weir,
    ho_captured_weir = ho_captured_weir,
    no_above_weir_uncaptured = no_above_weir_uncaptured,
    ho_above_weir_uncaptured = ho_above_weir_uncaptured,
    no_below_weir = no_below_weir,
    ho_below_weir = ho_below_weir
  )
}



calculate_pnob_slope <- function(scenario_inputs, anchor_result) {
  min(1, scenario_inputs$no_manarea_est / anchor_result$anchor_required_no)
}

allocate_final_brood <- function(scenario_inputs, weir_result, pnob_slope) {
  no_brood_desired <- pnob_slope * scenario_inputs$brood_need
  no_brood_actual <- min(no_brood_desired, weir_result$no_captured_weir)
  ho_brood_actual <- min(scenario_inputs$brood_need - no_brood_actual, weir_result$ho_captured_weir)
  
  list(
    pnob_slope = pnob_slope,
    no_brood_desired = no_brood_desired,
    no_brood_actual = no_brood_actual,
    ho_brood_actual = ho_brood_actual,
    no_captured_available_above = weir_result$no_captured_weir - no_brood_actual,
    ho_captured_available_above = weir_result$ho_captured_weir - ho_brood_actual
  )
}

calculate_final_spawners <- function(scenario_inputs, accounting_params, weir_result, brood_result) {
  no_spawners_above <-
    (weir_result$no_above_weir_uncaptured + brood_result$no_captured_available_above) *
    accounting_params$survival_above_weir
  
  ho_spawners_above_uncapt <-
    weir_result$ho_above_weir_uncaptured * accounting_params$survival_above_weir
  
  no_spawners_below <- weir_result$no_below_weir * accounting_params$survival_below_weir
  ho_spawners_below <- weir_result$ho_below_weir * accounting_params$survival_below_weir
  
  current_spawners_without_HO_captured_above <-
    no_spawners_above +
    no_spawners_below +
    ho_spawners_above_uncapt +
    ho_spawners_below
  
  spawner_deficit <- max(0, scenario_inputs$spawner_goal - current_spawners_without_HO_captured_above)
  ho_needed_pre_survival <- ifelse(
    accounting_params$survival_above_weir > 0,
    spawner_deficit / accounting_params$survival_above_weir,
    0
  )
  
  ho_captured_above <- min(
    brood_result$ho_captured_available_above,
    ho_needed_pre_survival
  )
  
  ho_captured_removed <- max(
    0,
    brood_result$ho_captured_available_above - ho_captured_above
  )
  
  ho_spawners_above <-
    ho_spawners_above_uncapt +
    ho_captured_above * accounting_params$survival_above_weir
  
  no_spawners_total <- no_spawners_above + no_spawners_below
  ho_spawners_total <- ho_spawners_above + ho_spawners_below
  system_spawners_total <- no_spawners_total + ho_spawners_total
  
  phos <- ifelse(system_spawners_total > 0, ho_spawners_total / system_spawners_total, NA_real_)
  pnob <- ifelse(scenario_inputs$brood_need > 0, brood_result$no_brood_actual / scenario_inputs$brood_need, NA_real_)
  pni <- ifelse(!is.na(phos) && !is.na(pnob) && (pnob + phos) > 0, pnob / (pnob + phos), NA_real_)
  
  list(
    no_spawners_above = no_spawners_above,
    ho_spawners_above_uncapt = ho_spawners_above_uncapt,
    no_spawners_below = no_spawners_below,
    ho_spawners_below = ho_spawners_below,
    current_spawners_without_HO_captured_above = current_spawners_without_HO_captured_above,
    spawner_deficit = spawner_deficit,
    ho_needed_pre_survival = ho_needed_pre_survival,
    ho_captured_above = ho_captured_above,
    ho_captured_removed = ho_captured_removed,
    ho_spawners_above = ho_spawners_above,
    no_spawners_total = no_spawners_total,
    ho_spawners_total = ho_spawners_total,
    system_spawners_total = system_spawners_total,
    phos = phos,
    pnob = pnob,
    pni = pni
  )
}
