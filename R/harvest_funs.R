# Purpose: Functions for calculating Lostine River sport and treaty
# fishery impacts and harvest.
#
# Author: Ryan N. Kinzer
#
# Fishery calculations are common to all weir management frameworks.
# Anchor, Recovery Phase, and Sliding Scale management rules are applied
# only after fisheries and weir accounting.


#------------------------------------------------------------
# Sport fishery NO impacts
#------------------------------------------------------------

impact_sport_wl <- function(NO_WL) {
  
  if (is.na(NO_WL)) {
    return(NA_real_)
  }
  
  # Historical Wallowa/Lostine sport fishery impact schedule.
  if (NO_WL < 300) {
    
    impact <- 0
    
  } else if (NO_WL < 1000) {
    
    impact <- 0.03 * (NO_WL - 300)
    
  } else if (NO_WL < 1500) {
    
    impact <-
      0.03 * (1000 - 300) +
      0.06 * (NO_WL - 1000)
    
  } else if (NO_WL < 2000) {
    
    impact <-
      0.03 * (1000 - 300) +
      0.06 * (1500 - 1000) +
      0.06 * (NO_WL - 1500)
    
  } else {
    
    impact <-
      0.03 * (1000 - 300) +
      0.06 * (1500 - 1000) +
      0.06 * (2000 - 1500) +
      0.12 * (NO_WL - 2000)
  }
  
  impact
}


#------------------------------------------------------------
# Treaty fishery NO impacts
#------------------------------------------------------------

impact_treaty_wl_direct <- function(NO_WL) {
  
  if (is.na(NO_WL)) {
    return(NA_real_)
  }
  
  # Historical Wallowa/Lostine treaty fishery impact schedule.
  if (NO_WL < 300) {
    
    impact <- 0.01 * NO_WL
    
  } else if (NO_WL < 1000) {
    
    impact <-
      0.01 * 300 +
      0.08 * (NO_WL - 300)
    
  } else if (NO_WL < 1500) {
    
    impact <-
      0.01 * 300 +
      0.08 * (1000 - 300) +
      0.16 * (NO_WL - 1000)
    
  } else if (NO_WL < 2000) {
    
    impact <-
      0.01 * 300 +
      0.08 * (1000 - 300) +
      0.16 * (1500 - 1000) +
      0.19 * (NO_WL - 1500)
    
  } else {
    
    impact <-
      0.01 * 300 +
      0.08 * (1000 - 300) +
      0.16 * (1500 - 1000) +
      0.19 * (2000 - 1500) +
      0.28 * (NO_WL - 2000)
  }
  
  impact
}


#------------------------------------------------------------
# Treaty HO harvest balance
#------------------------------------------------------------

calculate_ho_treaty_harvest_balance <- function(
    sport_impact_wl,
    ho_sport_harvest,
    treaty_no_impacts
) {
  
  # Total sport harvest burden represented by the sport fishery.
  #
  # sport_impact_wl represents natural-origin mortality from the
  # sport fishery. At 10% hooking mortality, the corresponding
  # natural-origin handle is 10 times the impact.
  sport_no_handle <-
    10 * sport_impact_wl
  
  sport_total_burden <-
    sport_no_handle +
    ho_sport_harvest
  
  # Treaty harvest is balanced against the sport fishery burden,
  # after accounting for treaty impacts to natural-origin fish.
  ho_treaty_harvest <- round_fish(
    max(
      0,
      sport_total_burden -
        treaty_no_impacts
    )
  )
  
  ho_treaty_harvest
}


#------------------------------------------------------------
# Scale values proportionally to a cap
#------------------------------------------------------------

scale_to_cap <- function(
    values,
    cap
) {
  
  values[is.na(values)] <- 0
  
  values <-
    pmax(
      0,
      values
    )
  
  total <-
    sum(values)
  
  if (
    is.na(cap) ||
    cap <= 0
  ) {
    return(
      rep(
        0,
        length(values)
      )
    )
  }
  
  if (
    total <= cap
  ) {
    return(values)
  }
  
  values *
    cap /
    total
}


#------------------------------------------------------------
# Common fishery calculation
#------------------------------------------------------------

#------------------------------------------------------------
# Common fishery calculation
#------------------------------------------------------------

calculate_fisheries <- function(
    scenario_inputs
) {
  
  # ------------------------------------------------------------
  # Starting abundance
  # ------------------------------------------------------------
  
  no_available <- round_fish(
    max(
      0,
      scenario_inputs$no_manarea_est
    )
  )
  
  ho_available <- round_fish(
    max(
      0,
      scenario_inputs$ho_manarea_est
    )
  )
  
  
  # ------------------------------------------------------------
  # Wallowa/Lostine NO abundance
  # ------------------------------------------------------------
  
  no_wl <-
    scenario_inputs$wl_scaling *
    scenario_inputs$no_manarea_est
  
  
  # ------------------------------------------------------------
  # Sport NO impacts
  # ------------------------------------------------------------
  
  sport_impact_wl <-
    max(
      0,
      impact_sport_wl(no_wl)
    )
  
  # Convert Wallowa/Lostine sport NO impacts back to the
  # Lostine management-area scale.
  sport_no_impacts_calculated <- round_fish(
    max(
      0,
      sport_impact_wl /
        scenario_inputs$wl_scaling
    )
  )
  
  # Sport impacts cannot exceed available NO adults.
  sport_no_impacts <- min(
    sport_no_impacts_calculated,
    no_available
  )
  
  
  # ------------------------------------------------------------
  # Treaty NO impacts
  # ------------------------------------------------------------
  
  treaty_no_impacts_calculated <- round_fish(
    max(
      0,
      impact_treaty_wl_direct(no_wl)
    )
  )
  
  # Treaty impacts occur after sport impacts and therefore cannot
  # exceed the NO adults remaining after the sport fishery.
  no_remaining_after_sport <- max(
    0,
    no_available -
      sport_no_impacts
  )
  
  treaty_no_impacts <- min(
    treaty_no_impacts_calculated,
    no_remaining_after_sport
  )
  
  total_no_impacts <-
    sport_no_impacts +
    treaty_no_impacts
  
  
  # ------------------------------------------------------------
  # HO sport harvest
  # ------------------------------------------------------------
  
  # Sport NO impacts represent 10% hooking mortality.
  #
  # Therefore:
  #
  # NO handle = NO impact / 0.10
  #
  # Expected HO harvest is then calculated from the HO:NO
  # abundance ratio.
  
  if (
    is.na(scenario_inputs$ho_manarea_est) ||
    scenario_inputs$ho_manarea_est < 20 ||
    is.na(sport_no_impacts) ||
    sport_no_impacts <= 0 ||
    scenario_inputs$no_manarea_est <= 0
  ) {
    
    ho_sport_harvest_calculated <- 0
    ho_sport_harvest <- 0
    
    sport_closed <- TRUE
    
  } else {
    
    allowed_no_handle <-
      10 * sport_no_impacts
    
    ho_sport_harvest_calculated <- round_fish(
      max(
        0,
        allowed_no_handle *
          (
            scenario_inputs$ho_manarea_est /
              scenario_inputs$no_manarea_est
          )
      )
    )
    
    # Sport harvest cannot exceed available HO adults.
    ho_sport_harvest <- min(
      ho_sport_harvest_calculated,
      ho_available
    )
    
    sport_closed <- FALSE
  }
  
  
  # ------------------------------------------------------------
  # HO treaty harvest
  # ------------------------------------------------------------
  
  ho_treaty_harvest_calculated <-
    calculate_ho_treaty_harvest_balance(
      sport_impact_wl = sport_impact_wl,
      ho_sport_harvest = ho_sport_harvest,
      treaty_no_impacts = treaty_no_impacts
    )
  
  # Treaty harvest occurs after sport harvest and therefore cannot
  # exceed the HO adults remaining after the sport fishery.
  ho_remaining_after_sport <- max(
    0,
    ho_available -
      ho_sport_harvest
  )
  
  ho_treaty_harvest <- min(
    ho_treaty_harvest_calculated,
    ho_remaining_after_sport
  )
  
  
  # ------------------------------------------------------------
  # Total HO harvest
  # ------------------------------------------------------------
  
  total_ho_harvest <-
    ho_sport_harvest +
    ho_treaty_harvest
  
  
  # ------------------------------------------------------------
  # Return common fishery accounting
  # ------------------------------------------------------------
  
  list(
    
    no_wl =
      no_wl,
    
    sport_impact_wl =
      sport_impact_wl,
    
    sport_no_impacts =
      sport_no_impacts,
    
    treaty_no_impacts =
      treaty_no_impacts,
    
    total_no_impacts =
      total_no_impacts,
    
    ho_sport_harvest =
      ho_sport_harvest,
    
    ho_treaty_harvest =
      ho_treaty_harvest,
    
    total_ho_harvest =
      total_ho_harvest,
    
    sport_closed =
      sport_closed
  )
}