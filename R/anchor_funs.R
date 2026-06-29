# Purpose: The functions within this script calculate the "anchor".
# Original Author: Kyle Bratcher
# Modified: Ryan Kinzer


anchor_eval <- function(est, anchor_inputs, accounting_params) {
  no_wl <- anchor_inputs$wl_scaling * est
  sport_wl <- impact_sport_wl(no_wl)
  sport_no <- sport_wl / anchor_inputs$wl_scaling # is the scaling factor the same for Wallowa expansion and harvest impacts
  treaty_no <- impact_treaty_wl_direct(no_wl)
  
  no_after_fisheries <- est - sport_no - treaty_no
  
  no_captured_weir <- accounting_params$trap_prop_no * no_after_fisheries # why is this not the weir efficiency
  
  no_above_weir_uncaptured <- calculate_uncaptured_above(
    captured = no_captured_weir,
    weir_efficiency = accounting_params$weir_efficiency
  )
  
  no_below_weir <- no_after_fisheries - (no_captured_weir + no_above_weir_uncaptured)
  
  # Anchor assumes broodstock is supplied entirely by NO fish.
  no_brood_at_anchor <- min(anchor_inputs$brood_need, no_captured_weir)
  no_captured_available_above_anchor <- no_captured_weir - no_brood_at_anchor
  
  no_spawners_above_anchor <-   # NO fish released above and surviving to spawning
    (no_above_weir_uncaptured + no_captured_available_above_anchor) *
    accounting_params$survival_above_weir
  
  no_spawners_below_anchor <- no_below_weir * accounting_params$survival_below_weir  # NO fish below the weir surviving to spawn
  no_spawners_total_anchor <- no_spawners_above_anchor + no_spawners_below_anchor # total NO successful spawners
  
  brood_deficit <- anchor_inputs$brood_need - no_captured_weir  #this does or does not include the weir efficiency 
  spawner_deficit <- anchor_inputs$spawner_goal - no_spawners_total_anchor
  
  list(
    candidate_no_est = est,
    sport_wl = sport_wl,
    sport_no = sport_no,
    treaty_no = treaty_no,
    no_after_fisheries = no_after_fisheries,
    no_captured_weir = no_captured_weir,
    no_above_weir_uncaptured = no_above_weir_uncaptured,
    no_below_weir = no_below_weir,
    no_brood_at_anchor = no_brood_at_anchor,
    no_captured_available_above_anchor = no_captured_available_above_anchor,
    no_spawners_above_anchor = no_spawners_above_anchor,
    no_spawners_below_anchor = no_spawners_below_anchor,
    no_spawners_total_anchor = no_spawners_total_anchor,
    brood_deficit = brood_deficit,
    spawner_deficit = spawner_deficit
  )
}

anchor_feasibility <- function(est, anchor_inputs, accounting_params) {
  ev <- anchor_eval(est, anchor_inputs, accounting_params)
  # Returns a value whose SIGN indicates feasibility:
  #   > 0  both objectives met (brood need and spawner goal satisfied)
  #   = 0  exactly at the feasibility boundary
  #   < 0  infeasible (at least one objective not yet met)
  #
  # Constructed as min(-brood_deficit, -spawner_deficit):
  #   A deficit > 0 means a shortfall exists, so negating it gives < 0.
  #   Both deficits must be <= 0 (objectives met) for this to return >= 0.
  #
  # This sign convention is required by uniroot() in find_smallest_feasible_anchor(),
  # which searches for the crossing point where feasibility transitions from
  # negative to non-negative — i.e., the minimum NO abundance that satisfies
  # both objectives.
  min(-ev$brood_deficit, -ev$spawner_deficit)
}

find_smallest_feasible_anchor <- function(
    anchor_inputs,
    accounting_params,
    start_est = 0,
    step = 5,
    max_est = 5000
) {
  est <- max(0, start_est)
  
  while (anchor_feasibility(est, anchor_inputs, accounting_params) < 0 && est < max_est) {
    est <- est + step
  }
  
  if (anchor_feasibility(est, anchor_inputs, accounting_params) < 0) {
    stop("Anchor solver: no feasible solution found up to max_est.", call. = FALSE)
  }
  
  lower <- max(est - 100, 0)
  upper <- est
  
  while (anchor_feasibility(lower, anchor_inputs, accounting_params) >= 0 && lower > 0) {
    lower <- max(0, lower - 10)
  }
  
  if (anchor_feasibility(lower, anchor_inputs, accounting_params) < 0 &&
      anchor_feasibility(upper, anchor_inputs, accounting_params) >= 0) {
    rt <- uniroot(
      function(x) anchor_feasibility(x, anchor_inputs, accounting_params),
      interval = c(lower, upper)
    )
    return(rt$root)
  }
  
  est
}