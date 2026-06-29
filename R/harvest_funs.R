# Purpose: The functions within this script calculate sport and tribal harvest
# impacts.
# Original Author: Kyle Bratcher
# Modified: Ryan Kinzer updated typos within the natural origin impact
# calculations.


# Impact to state (sport) fishery on natural-origin Wallowa/Lostine fish
# Based on Table 7, FMEP (ODFW/WDFW 2012)

# Critical threshold = 300, MAT = 1,000

# Rate structure (state column):
#   Scenario A: < 300        -> 0%
#   Scenario B: 300-1000     -> 3% of margin above critical (300)
#   Scenario C: 1000-1500    -> B + 6% of margin above MAT (1000)
#   Scenario D: 1500-2000    -> C + 6% of margin above MAT (1000) [same rate as C; explicit bracket retained for clarity]
#   Scenario E: > 2000       -> D + 12% of margin above 2x MAT (2000)


impact_sport_wl <- function(NO_WL) {
  if (is.na(NO_WL)) return(0)
  if (NO_WL < 300)  return(0)
  if (NO_WL < 1000) return(round((NO_WL - 300) * 0.03))
  if (NO_WL < 1500) return(round(21 + (NO_WL - 1000) * 0.06))
  if (NO_WL < 2000) return(round(51 + (NO_WL - 1500) * 0.06))
  round(81 + (NO_WL - 2000) * 0.12)
}

# Rate structure (tribal column):
#   Scenario A: < 300        -> 1% of all fish
#   Scenario B: 300-1000     -> A + 8% of margin above critical (300)
#   Scenario C: 1000-1500    -> B + 16% of margin above MAT (1000)
#   Scenario D: 1500-2000    -> C + 19% of margin above 1.5x MAT (1500)
#   Scenario E: > 2000       -> D + 28% of margin above 2x MAT (2000)

# Treaty NO impacts are calculated from WL-scaled abundance.
impact_treaty_wl_direct <- function(NO_WL) {
  if (is.na(NO_WL)) return(0)
  if (NO_WL < 300)  return(round(NO_WL * 0.01))
  if (NO_WL < 1000) return(round(3 + (NO_WL - 300) * 0.08))
  if (NO_WL < 1500) return(round(59 + (NO_WL - 1000) * 0.16))
  if (NO_WL < 2000) return(round(139 + (NO_WL - 1500) * 0.19))
  round(234 + (NO_WL - 2000) * 0.28)
}


# HO treaty harvest is a catch-balancing residual, not a separate HO
# treaty harvest curve.
#
# Background:
# The co-manager fishery framework operates on a 50/50 total harvest share
# between tribal and non-tribal fisheries (FMEP section 2.1.3, Table 7).
# The sport fishery (non-tribal) generates two types of fish encounters:
#   - NO fish caught and released incidentally (sport_impact_wl, on WL scale)
#   - HO fish harvested by sport anglers (ho_sport_harvest)
# Together these define the total non-tribal fishery burden on WL-scale fish.
#
# The tribal fishery has already been accounted for in NO fish through
# treaty_no_impacts. To balance the allocation, the tribal side must also
# receive an equivalent HO harvest. That equivalent is the residual:
#   HO treaty harvest = (sport_impact_wl + ho_sport_harvest) - treaty_no_impacts
#
# In plain terms: if the non-tribal sport fishery totals X fish (NO impacts +
# HO harvest combined), and the tribal fishery has already taken Y NO fish,
# then the tribal side is owed (X - Y) additional HO fish to maintain balance.
#
# Capped at zero so harvest cannot become negative in years when treaty NO
# impacts already exceed the total sport fishery value.

calculate_ho_treaty_harvest_balance <- function(sport_impact_wl, ho_sport_harvest, treaty_no_impacts) {
  max(0, (sport_impact_wl + ho_sport_harvest) - treaty_no_impacts)
}

# Scale a set of projected impacts/harvests down to a cap while preserving
# the projected split. Values and cap are constrained to be non-negative.
scale_to_cap <- function(values, cap) {
  values[is.na(values)] <- 0
  values <- pmax(0, values)
  
  if (is.na(cap)) cap <- 0
  cap <- max(0, cap)
  
  total <- sum(values)
  
  if (total <= 0 || cap <= 0) {
    values[] <- 0
    return(values)
  }
  
  if (total <= cap) {
    return(values)
  }
  
  values * (cap / total)
}


# uses the above functions to calculate fisheries outcomes
calculate_scenario_fisheries <- function(scenario_inputs, accounting_params, anchor_result) {
  no_wl <- scenario_inputs$wl_scaling * scenario_inputs$no_manarea_est
  
  # Projected NO fishery values from the existing curves.
  sport_impact_wl_projected <- max(0, impact_sport_wl(no_wl))
  sport_no_impacts_projected <- max(
    0,
    sport_impact_wl_projected / scenario_inputs$wl_scaling
  )
  treaty_no_impacts_projected <- max(0, impact_treaty_wl_direct(no_wl))
  total_no_impacts_projected <-
    sport_no_impacts_projected + treaty_no_impacts_projected
  
  # Broodstock is the only priority protected ahead of fisheries in Part B.
  # The Anchor is used here only to calculate the pNOB line, not as a fishery cap.
  pnob_slope_for_cap <- min(
    1,
    scenario_inputs$no_manarea_est / anchor_result$anchor_required_no
  )
  
  no_brood_target <- pnob_slope_for_cap * scenario_inputs$brood_need
  
  required_no_after_fisheries <- if (
    !is.na(accounting_params$trap_prop_no) &&
    accounting_params$trap_prop_no > 0
  ) {
    no_brood_target / accounting_params$trap_prop_no
  } else {
    Inf
  }
  
  max_no_impacts_allowed <- max(
    0,
    scenario_inputs$no_manarea_est - required_no_after_fisheries
  )
  
  no_impacts_allowed <- scale_to_cap(
    values = c(sport_no_impacts_projected, treaty_no_impacts_projected),
    cap = max_no_impacts_allowed
  )
  
  sport_no_impacts <- no_impacts_allowed[1]
  treaty_no_impacts <- no_impacts_allowed[2]
  total_no_impacts_allowed <- sport_no_impacts + treaty_no_impacts
  
  # Convert the allowed Lostine sport NO impact back to WL scale for the
  # HO treaty harvest balance calculation.
  sport_impact_wl <- max(0, sport_no_impacts * scenario_inputs$wl_scaling)
  
  # HO sport harvest is derived from the allowed NO sport impact using two steps.
  #
  # Step 1 — Infer total allowed angler encounters with NO fish.
  #   The sport fishery targets adipose-clipped hatchery fish; wild (NO) fish
  #   are encountered incidentally and released. The FMEP (section 1.4.3)
  #   establishes 10% as the assumed catch-and-release hooking mortality rate
  #   for this fishery. sport_no_impacts represents the NO fish expected to die
  #   from incidental encounters. Working backwards:
  #     allowed_no_handle = sport_no_impacts / 0.10  =  sport_no_impacts * 10
  #   This gives the total number of NO fish that can be handled (encountered
  #   and released) while staying within the allowed NO mortality.
  #
  # Step 2 — Scale to HO fish using the HO:NO abundance ratio.
  #   Fish encounter rates in the fishery are assumed proportional to abundance.
  #   If the HO:NO ratio is R, then for every NO fish encountered, R HO fish
  #   are expected to be encountered and harvested (HO fish are kept, not
  #   released). Therefore:
  #     ho_sport_harvest_projected = allowed_no_handle * (ho_manarea_est / no_manarea_est)
  #
  # This approach anchors HO sport harvest to the allowed NO impact rather than
  # running a separate HO harvest curve, which keeps the method internally
  # consistent with the existing FMEP fishery framework.
  
  if (
    is.na(scenario_inputs$ho_manarea_est) ||
    scenario_inputs$ho_manarea_est < 20 ||
    is.na(sport_no_impacts) ||
    sport_no_impacts <= 0
  ) {
    ho_sport_harvest_projected <- 0
    sport_closed <- TRUE
  } else {
    allowed_no_handle <- 10 * sport_no_impacts
    ho_sport_harvest_projected <- max(
      0,
      round(
        allowed_no_handle *
          (scenario_inputs$ho_manarea_est / scenario_inputs$no_manarea_est)
      )
    )
    sport_closed <- FALSE
  }
  
  # HO treaty harvest is a catch-balancing residual from the total sport
  # fishery on the WL-scale framework. This is projected before the HO brood cap.
  ho_treaty_harvest_projected <- calculate_ho_treaty_harvest_balance(
    sport_impact_wl = sport_impact_wl,
    ho_sport_harvest = ho_sport_harvest_projected,
    treaty_no_impacts = treaty_no_impacts
  )
  total_ho_harvest_projected <-
    ho_sport_harvest_projected + ho_treaty_harvest_projected
  
  # Estimate remaining brood need after allowed NO fisheries.
  no_after_fisheries_for_brood_check <-
    scenario_inputs$no_manarea_est - sport_no_impacts - treaty_no_impacts
  
  no_captured_for_brood_check <-
    accounting_params$trap_prop_no * no_after_fisheries_for_brood_check
  
  no_brood_actual_for_brood_check <- min(
    no_brood_target,
    no_captured_for_brood_check
  )
  
  ho_brood_needed <- max(
    0,
    scenario_inputs$brood_need - no_brood_actual_for_brood_check
  )
  
  required_ho_after_fisheries <- if (
    !is.na(accounting_params$trap_prop_ho) &&
    accounting_params$trap_prop_ho > 0
  ) {
    ho_brood_needed / accounting_params$trap_prop_ho
  } else {
    Inf
  }
  
  max_ho_harvest_allowed <- max(
    0,
    scenario_inputs$ho_manarea_est - required_ho_after_fisheries
  )
  
  ho_harvest_allowed <- scale_to_cap(
    values = c(ho_sport_harvest_projected, ho_treaty_harvest_projected),
    cap = max_ho_harvest_allowed
  )
  
  ho_sport_harvest <- ho_harvest_allowed[1]
  ho_treaty_harvest <- ho_harvest_allowed[2]
  total_ho_harvest_allowed <- ho_sport_harvest + ho_treaty_harvest
  
  list(
    no_wl = no_wl,
    sport_closed = sport_closed,
    
    sport_impact_wl_projected = sport_impact_wl_projected,
    sport_no_impacts_projected = sport_no_impacts_projected,
    treaty_no_impacts_projected = treaty_no_impacts_projected,
    total_no_impacts_projected = total_no_impacts_projected,
    
    pnob_slope_for_cap = pnob_slope_for_cap,
    no_brood_target = no_brood_target,
    required_no_after_fisheries = required_no_after_fisheries,
    max_no_impacts_allowed = max_no_impacts_allowed,
    
    sport_impact_wl = sport_impact_wl,
    sport_no_impacts = sport_no_impacts,
    treaty_no_impacts = treaty_no_impacts,
    total_no_impacts_allowed = total_no_impacts_allowed,
    
    ho_sport_harvest_projected = ho_sport_harvest_projected,
    ho_treaty_harvest_projected = ho_treaty_harvest_projected,
    total_ho_harvest_projected = total_ho_harvest_projected,
    
    ho_brood_needed = ho_brood_needed,
    required_ho_after_fisheries = required_ho_after_fisheries,
    max_ho_harvest_allowed = max_ho_harvest_allowed,
    
    ho_sport_harvest = ho_sport_harvest,
    sport_total_for_balance = sport_impact_wl + ho_sport_harvest,
    ho_treaty_harvest = ho_treaty_harvest,
    total_ho_harvest_allowed = total_ho_harvest_allowed
  )
}