fmt <- function(x) round(x, disp_digits)

round_fish <- function(x) {
  floor(x + 0.5)
}

y2 <- function(y) sprintf("%02d", as.integer(y) %% 100)

first5_non_na <- function(values, years) {
  idx <- which(!is.na(values))
  idx <- idx[seq_len(min(5, length(idx)))]
  tibble(Year = years[idx], Value = values[idx])
}

check_median <- function(df, name) {
  if (nrow(df) < 3) {
    stop(
      sprintf(
        "Median '%s' has fewer than 3 valid years. Please update input_chs.csv.",
        name
      ),
      call. = FALSE
    )
  }
}

years_used_str <- function(df) {
  paste(y2(df$Year), collapse = ",")
}


# Prefer Post -> In -> Pre when optionally defaulting abundance values from CSV.
pick_est <- function(post, inseason, pre) {
  if (!is.na(post)) {
    post
  } else if (!is.na(inseason)) {
    inseason
  } else {
    pre
  }
}


# Helper for year-specific actual values.
# If the selected year has an actual measured value, use it.
# If the selected year is missing that value, fall back to the median-derived value.
pick_actual_or_median <- function(
    actual_value,
    median_value,
    actual_label,
    median_label
) {
  
  if (!is.na(actual_value)) {
    
    list(
      value = actual_value,
      source = actual_label
    )
    
  } else {
    
    list(
      value = median_value,
      source = median_label
    )
  }
}


#------------------------------------------------------------
# Evaluate Anchor sensitivity to spawner goal
#------------------------------------------------------------

evaluate_anchor_spawner_goal <- function(
    spawner_goals,
    anchor_inputs,
    accounting_params,
    start_est = 0
) {
  
  tibble(
    spawner_goal = spawner_goals,
    
    anchor = purrr::map_dbl(
      spawner_goal,
      function(goal) {
        
        tmp_inputs <- anchor_inputs
        tmp_inputs$spawner_goal <- goal
        
        find_smallest_feasible_anchor(
          anchor_inputs = tmp_inputs,
          accounting_params = accounting_params,
          start_est = start_est
        )
      }
    )
  )
}


#------------------------------------------------------------
# Create common scenario inputs
#------------------------------------------------------------

create_scenario_inputs <- function(
    no_manarea_est,
    ho_manarea_est,
    brood_need = 160,
    spawner_goal = 800,
    wl_scaling = 1.4,
    utilization_no = 1,
    utilization_ho = 1
) {
  
  list(
    no_manarea_est = as.numeric(no_manarea_est),
    ho_manarea_est = as.numeric(ho_manarea_est),
    brood_need = as.numeric(brood_need),
    spawner_goal = as.numeric(spawner_goal),
    wl_scaling = as.numeric(wl_scaling),
    utilization_no = as.numeric(utilization_no),
    utilization_ho = as.numeric(utilization_ho)
  )
}


#------------------------------------------------------------
# Calculate Anchor abundance
#------------------------------------------------------------

calculate_anchor_result <- function(
    scenario_inputs,
    anchor_inputs,
    accounting_params,
    start_est = 0
) {
  
  ai <- anchor_inputs
  
  ai$brood_need <-
    scenario_inputs$brood_need
  
  ai$spawner_goal <-
    scenario_inputs$spawner_goal
  
  ai$wl_scaling <-
    scenario_inputs$wl_scaling
  
  anchor_required_no <- find_smallest_feasible_anchor(
    anchor_inputs = ai,
    accounting_params = accounting_params,
    start_est = start_est
  )
  
  anchor_evaluation <- anchor_eval(
    est = anchor_required_no,
    anchor_inputs = ai,
    accounting_params = accounting_params
  )
  
  list(
    anchor_required_no = anchor_required_no,
    brood_need = ai$brood_need,
    spawner_goal = ai$spawner_goal,
    wl_scaling = ai$wl_scaling,
    full_anchor_solution = anchor_evaluation
  )
}


#------------------------------------------------------------
# Run complete Anchor accounting scenario
#------------------------------------------------------------

run_accounting_scenario <- function(
    scenario_inputs,
    accounting_params,
    anchor_inputs,
    start_est = 0
) {
  
  # ------------------------------------------------------------
  # 1. Calculate Anchor
  # ------------------------------------------------------------
  
  anchor_result <- calculate_anchor_result(
    scenario_inputs = scenario_inputs,
    anchor_inputs = anchor_inputs,
    accounting_params = accounting_params,
    start_est = start_est
  )
  
  
  # ------------------------------------------------------------
  # 2. Calculate common fisheries
  # ------------------------------------------------------------
  
  # Fisheries are independent of the weir management framework.
  # Anchor, Recovery Phase, and Sliding Scale scenarios all use
  # the same fishery calculation.
  
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
  # 4. Calculate Anchor pNOB relationship
  # ------------------------------------------------------------
  
  pnob_slope_result <- calculate_pnob_slope(
    scenario_inputs = scenario_inputs,
    anchor_result = anchor_result
  )
  
  
  # ------------------------------------------------------------
  # 5. Allocate broodstock
  # ------------------------------------------------------------
  
  brood_result <- allocate_final_brood(
    scenario_inputs,
    weir_result,
    pnob_slope_result
  )
  
  
  # ------------------------------------------------------------
  # 6. Allocate natural spawners
  # ------------------------------------------------------------
  
  final_metrics <- calculate_final_spawners(
    scenario_inputs = scenario_inputs,
    accounting_params = accounting_params,
    weir_result = weir_result,
    brood_result = brood_result
  )
  
  
  # ------------------------------------------------------------
  # Return complete accounting
  # ------------------------------------------------------------
  
  list(
    scenario_inputs = scenario_inputs,
    anchor_result = anchor_result,
    fishery_result = fishery_result,
    weir_result = weir_result,
    pnob_slope_result = pnob_slope_result,
    brood_result = brood_result,
    final_metrics = final_metrics
  )
}


#------------------------------------------------------------
# Tidy Anchor output
#------------------------------------------------------------

tidy_accounting_result <- function(result) {
  
  ar <- result$anchor_result
  fr <- result$fishery_result
  br <- result$brood_result
  fm <- result$final_metrics
  
  tibble(
    
    # Anchor framework
    anchor_required_no =
      ar$anchor_required_no,
    
    
    # Fisheries
    sport_NO_impact =
      fr$sport_no_impacts,
    
    treaty_NO_impact =
      fr$treaty_no_impacts,
    
    total_NO_impact =
      fr$total_no_impacts,
    
    HO_sport_harvest =
      fr$ho_sport_harvest,
    
    HO_treaty_harvest =
      fr$ho_treaty_harvest,
    
    total_HO_harvest =
      fr$total_ho_harvest,
    
    
    # Broodstock
    NO_brood =
      br$no_brood_actual,
    
    HO_brood =
      br$ho_brood_actual,
    
    total_brood =
      br$no_brood_actual +
      br$ho_brood_actual,
    
    
    # Natural spawning
    HO_removed =
      fm$ho_captured_removed,
    
    NO_spawners =
      fm$no_spawners_total,
    
    HO_spawners =
      fm$ho_spawners_total,
    
    system_spawners =
      fm$system_spawners_total,
    
    
    # Realized hatchery influence
    pHOS =
      fm$phos,
    
    pNOB =
      fm$pnob,
    
    PNI =
      fm$pni
  )
}


#------------------------------------------------------------
# Evaluate Anchor scenario from abundance inputs
#------------------------------------------------------------

evaluate_accounting_scenario <- function(
    no_manarea_est,
    ho_manarea_est,
    accounting_params,
    anchor_inputs,
    brood_need = 160,
    spawner_goal = 800,
    wl_scaling = 1.4,
    utilization_no = 1,
    utilization_ho = 1,
    start_est = 0
) {
  
  scenario_inputs <- create_scenario_inputs(
    no_manarea_est = no_manarea_est,
    ho_manarea_est = ho_manarea_est,
    brood_need = brood_need,
    spawner_goal = spawner_goal,
    wl_scaling = wl_scaling,
    utilization_no = utilization_no,
    utilization_ho = utilization_ho
  )
  
  result <- run_accounting_scenario(
    scenario_inputs = scenario_inputs,
    accounting_params = accounting_params,
    anchor_inputs = anchor_inputs,
    start_est = start_est
  )
  
  tidy_accounting_result(result)
}