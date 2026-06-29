##############################################
# Lostine Management Area Accounting Toolkit #
# Reviewer Refactor: Anchor + Accounting      #
##############################################

# ============================================================
# 0) CONFIG, LIBRARIES, AND USER OPTIONS
# ============================================================
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(knitr)   # for kable prints
})

# Display rounding only. Calculations retain full precision.
disp_digits <- 2
fmt <- function(x) round(x, disp_digits)

# Historical input file used to calculate median-derived parameters.
input_file <- "input_chs.csv"
# Optional: default scenario values from a CSV year row.
# For reviewer workflow, leave FALSE and edit the manual values below.
use_csv_scenario_values <- FALSE
csv_scenario_year <- 2023

# Print reviewer tables in R.
print_tables <- TRUE

# Utility: 2-digit year string for table display.
y2 <- function(y) sprintf("%02d", as.integer(y) %% 100)

message(sprintf("Lostine Toolkit — Historical file: %s", input_file))


# ============================================================
# 1) HISTORICAL DATA AND MEDIAN-DERIVED PARAMETERS
# ============================================================
# The historical CSV is used here only to calculate median-derived
# accounting parameters. Scenario/in-year values are entered later.

dat <- read_csv(input_file, show_col_types = FALSE) %>%
  filter(.data$Stock == "LST") %>%
  arrange(desc(.data$Year))

first5_non_na <- function(values, years) {
  idx <- which(!is.na(values))
  idx <- idx[seq_len(min(5, length(idx)))]
  tibble(Year = years[idx], Value = values[idx])
}

check_median <- function(df, name) {
  if (nrow(df) < 3) {
    stop(sprintf("Median '%s' has fewer than 3 valid years. Please update input_chs.csv.", name), call. = FALSE)
  }
}

years_used_str <- function(df) paste(y2(df$Year), collapse = ",")

# Build trap ratio columns safely. Division by zero becomes NA.
dat <- dat %>%
  mutate(
    trap_ratio_no = if_else(Post.NO == 0, NA_real_, Trap.NO / Post.NO),
    trap_ratio_ho = if_else(Post.HO == 0, NA_real_, Trap.HO / Post.HO)
  )

m_trap_no <- first5_non_na(dat$trap_ratio_no, dat$Year)
m_trap_ho <- first5_non_na(dat$trap_ratio_ho, dat$Year)
m_weir    <- first5_non_na(dat$Weir.Eff,     dat$Year)
m_abv     <- first5_non_na(dat$PSS.ABV,      dat$Year)
m_blw     <- first5_non_na(dat$PSS.BLW,      dat$Year)

purrr::walk2(
  list(m_trap_no, m_trap_ho, m_weir, m_abv, m_blw),
  c("trap_prop_no", "trap_prop_ho", "weir_efficiency", "survival_above_weir", "survival_below_weir"),
  check_median
)

median_params <- list(
  trap_prop_no        = median(m_trap_no$Value),
  trap_prop_ho        = median(m_trap_ho$Value),
  weir_efficiency     = median(m_weir$Value),
  survival_above_weir = median(m_abv$Value),
  survival_below_weir = median(m_blw$Value),
  years_trap_no       = years_used_str(m_trap_no),
  years_trap_ho       = years_used_str(m_trap_ho),
  years_weir          = years_used_str(m_weir),
  years_abv           = years_used_str(m_abv),
  years_blw           = years_used_str(m_blw)
)


# ============================================================
# 2) MANUAL SCENARIO / ANCHOR INPUTS
# ============================================================
# For review, edit these values directly. The scenario accounting is not
# tied to a specific management year unless use_csv_scenario_values is TRUE.

scenario_inputs <- list(
  no_manarea_est = 382,  # current/scenario natural-origin abundance estimate
  ho_manarea_est = 536,   # current/scenario hatchery-origin abundance estimate
  brood_need     = 160,   # broodstock need
  spawner_goal   = 800,   # adult spawner goal
  wl_scaling     = 1.4,    # Wild-Lostine scaling factor
  utilization_no  = 1.0,   # proportion of allowed NO impacts actually taken (1.0 = full utilization)
  utilization_ho  = 1.0    # proportion of allowed HO harvest actually taken (1.0 = full utilization)
)

# Prefer Post -> In -> Pre when optionally defaulting abundance values from CSV.
pick_est <- function(post, inseason, pre) {
  if (!is.na(post)) post else if (!is.na(inseason)) inseason else pre
}

# Helper for year-specific actual values.
# If the selected year has an actual measured value, use it.
# If the selected year is missing that value, fall back to the median-derived value.
pick_actual_or_median <- function(actual_value, median_value, actual_label, median_label) {
  if (!is.na(actual_value)) {
    list(value = actual_value, source = actual_label)
  } else {
    list(value = median_value, source = median_label)
  }
}

# Default accounting parameters are the median-derived values.
# These are used for normal planning/review scenarios when current-year actuals
# are unavailable or when use_csv_scenario_values is FALSE.
accounting_params <- list(
  trap_prop_no        = median_params$trap_prop_no,
  trap_prop_ho        = median_params$trap_prop_ho,
  weir_efficiency     = median_params$weir_efficiency,
  survival_above_weir = median_params$survival_above_weir,
  survival_below_weir = median_params$survival_below_weir,
  source_trap_no      = paste0("Median years: ", median_params$years_trap_no),
  source_trap_ho      = paste0("Median years: ", median_params$years_trap_ho),
  source_weir         = paste0("Median years: ", median_params$years_weir),
  source_abv          = paste0("Median years: ", median_params$years_abv),
  source_blw          = paste0("Median years: ", median_params$years_blw)
)

if (use_csv_scenario_values) {
  row_scenario <- dat %>% filter(.data$Year == csv_scenario_year) %>% slice(1)

  if (nrow(row_scenario) == 0) {
    stop("No LST row found for csv_scenario_year.", call. = FALSE)
  }

  # Scenario values from the selected year where available.
  scenario_inputs$no_manarea_est <- pick_est(row_scenario$Post.NO, row_scenario$In.NO, row_scenario$Pre.NO)
  scenario_inputs$ho_manarea_est <- pick_est(row_scenario$Post.HO, row_scenario$In.HO, row_scenario$Pre.HO)
  scenario_inputs$brood_need     <- row_scenario$BS.Need

  # spawner_goal and wl_scaling remain manual unless those values are added to the CSV.

  actual_year_label <- paste0("Actual year: ", csv_scenario_year)

  trap_no_choice <- pick_actual_or_median(
    actual_value = row_scenario$trap_ratio_no,
    median_value = median_params$trap_prop_no,
    actual_label = actual_year_label,
    median_label = paste0("Median fallback years: ", median_params$years_trap_no)
  )

  trap_ho_choice <- pick_actual_or_median(
    actual_value = row_scenario$trap_ratio_ho,
    median_value = median_params$trap_prop_ho,
    actual_label = actual_year_label,
    median_label = paste0("Median fallback years: ", median_params$years_trap_ho)
  )

  weir_choice <- pick_actual_or_median(
    actual_value = row_scenario$Weir.Eff,
    median_value = median_params$weir_efficiency,
    actual_label = actual_year_label,
    median_label = paste0("Median fallback years: ", median_params$years_weir)
  )

  abv_choice <- pick_actual_or_median(
    actual_value = row_scenario$PSS.ABV,
    median_value = median_params$survival_above_weir,
    actual_label = actual_year_label,
    median_label = paste0("Median fallback years: ", median_params$years_abv)
  )

  blw_choice <- pick_actual_or_median(
    actual_value = row_scenario$PSS.BLW,
    median_value = median_params$survival_below_weir,
    actual_label = actual_year_label,
    median_label = paste0("Median fallback years: ", median_params$years_blw)
  )

  accounting_params <- list(
    trap_prop_no        = trap_no_choice$value,
    trap_prop_ho        = trap_ho_choice$value,
    weir_efficiency     = weir_choice$value,
    survival_above_weir = abv_choice$value,
    survival_below_weir = blw_choice$value,
    source_trap_no      = trap_no_choice$source,
    source_trap_ho      = trap_ho_choice$source,
    source_weir         = weir_choice$source,
    source_abv          = abv_choice$source,
    source_blw          = blw_choice$source
  )
}

# The Anchor method uses these same policy/scenario assumptions, but only NO fish.
anchor_inputs <- list(
  brood_need   = scenario_inputs$brood_need,
  spawner_goal = scenario_inputs$spawner_goal,
  wl_scaling   = scenario_inputs$wl_scaling
)

message(sprintf(
  "Scenario inputs — NO: %.3f  HO: %.3f  BroodNeed: %.1f  WL scaling: %.3f  SpawnerGoal: %.1f",
  scenario_inputs$no_manarea_est,
  scenario_inputs$ho_manarea_est,
  scenario_inputs$brood_need,
  scenario_inputs$wl_scaling,
  scenario_inputs$spawner_goal
))

message(sprintf(
  "Accounting parameters — trap NO: %.3f [%s]; trap HO: %.3f [%s]; weir: %.3f [%s]; above survival: %.3f [%s]; below survival: %.3f [%s]",
  accounting_params$trap_prop_no, accounting_params$source_trap_no,
  accounting_params$trap_prop_ho, accounting_params$source_trap_ho,
  accounting_params$weir_efficiency, accounting_params$source_weir,
  accounting_params$survival_above_weir, accounting_params$source_abv,
  accounting_params$survival_below_weir, accounting_params$source_blw
))


# ============================================================
# 3) FISHERY CURVES AND COMMON ACCOUNTING HELPERS
# ============================================================
# These functions are shared by Part A and Part B.

# Sport impacts are calculated on Wild-Lostine scale and then converted
# back to Lostine-only subtraction.
impact_sport_wl <- function(NO_WL) {
  if (is.na(NO_WL)) return(0)
  if (NO_WL < 300)  return(0)
  if (NO_WL < 1000) return(round((NO_WL - 300) * 0.03))
  if (NO_WL < 1500) return(round(21 + (NO_WL - 1000) * 0.06))
  if (NO_WL < 2000) return(round(51 + (NO_WL - 1500) * 0.06))
  round(81 + (NO_WL - 1500) * 0.12)
}

# Treaty NO impacts are calculated from WL-scaled abundance.
# Unlike sport NO impacts, treaty impacts are not scaled back to Lostine because
# treaty impacts occur within the Lostine River itself. Sport impacts occur in the
# Wallowa River, where fish from outside the Lostine are also encountered.
# Therefore, treaty NO impacts are subtracted directly from Lostine NO abundance
# by management convention.
impact_treaty_wl_direct <- function(NO_WL) {
  if (is.na(NO_WL)) return(0)
  if (NO_WL < 300)  return(round(NO_WL * 0.01))
  if (NO_WL < 1000) return(round(3 + (NO_WL - 300) * 0.08))
  if (NO_WL < 1500) return(round(59 + (NO_WL - 1000) * 0.16))
  if (NO_WL < 2000) return(round(139 + (NO_WL - 1500) * 0.19))
  round(234 + (NO_WL - 1999) * 0.28)
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

# Weir efficiency identity:
#   weir_efficiency = captured / (captured + uncaptured_above)
# Rearranged:
#   uncaptured_above = (captured / weir_efficiency) - captured
calculate_uncaptured_above <- function(captured, weir_efficiency) {
  if (is.na(captured) || is.na(weir_efficiency) || weir_efficiency <= 0) return(NA_real_)
  (captured / weir_efficiency) - captured
}


# ============================================================
# PART A — SET THE ANCHOR
# ============================================================
# Purpose:
# Determine the NO-only Anchor abundance: the minimum natural-origin
# abundance required to satisfy both broodstock need and spawner goal.
#
# This is a standalone methodological section. It does not calculate HO
# values, pHOS, pNOB, PNI, or final mixed-origin accounting.

anchor_eval <- function(est, anchor_inputs, accounting_params) {
  no_wl <- anchor_inputs$wl_scaling * est
  sport_wl <- impact_sport_wl(no_wl)
  sport_no <- sport_wl / anchor_inputs$wl_scaling
  treaty_no <- impact_treaty_wl_direct(no_wl)

  no_after_fisheries <- est - sport_no - treaty_no
  no_captured_weir <- accounting_params$trap_prop_no * no_after_fisheries
  no_above_weir_uncaptured <- calculate_uncaptured_above(
    captured = no_captured_weir,
    weir_efficiency = accounting_params$weir_efficiency
  )
  no_below_weir <- no_after_fisheries - (no_captured_weir + no_above_weir_uncaptured)

  # Anchor assumes broodstock is supplied entirely by NO fish.
  no_brood_at_anchor <- min(anchor_inputs$brood_need, no_captured_weir)
  no_captured_available_above_anchor <- no_captured_weir - no_brood_at_anchor

  no_spawners_above_anchor <-
    (no_above_weir_uncaptured + no_captured_available_above_anchor) *
    accounting_params$survival_above_weir

  no_spawners_below_anchor <- no_below_weir * accounting_params$survival_below_weir
  no_spawners_total_anchor <- no_spawners_above_anchor + no_spawners_below_anchor

  brood_deficit <- anchor_inputs$brood_need - no_captured_weir
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

anchor_required_no <- find_smallest_feasible_anchor(
  anchor_inputs = anchor_inputs,
  accounting_params = accounting_params,
  start_est = 0
)

anchor_at_solution <- anchor_eval(
  est = anchor_required_no,
  anchor_inputs = anchor_inputs,
  accounting_params = accounting_params
)

anchor_result <- list(
  anchor_required_no = anchor_required_no,
  brood_need = anchor_inputs$brood_need,
  spawner_goal = anchor_inputs$spawner_goal,
  wl_scaling = anchor_inputs$wl_scaling,
  trap_prop_no = accounting_params$trap_prop_no,
  weir_efficiency = accounting_params$weir_efficiency,
  survival_above_weir = accounting_params$survival_above_weir,
  survival_below_weir = accounting_params$survival_below_weir,
  no_captured_at_anchor = anchor_at_solution$no_captured_weir,
  no_spawners_at_anchor = anchor_at_solution$no_spawners_total_anchor,
  sport_impact_wl_at_anchor = anchor_at_solution$sport_wl,
  sport_no_impacts_at_anchor = anchor_at_solution$sport_no,
  treaty_no_impacts_at_anchor = anchor_at_solution$treaty_no,
  full_anchor_solution = anchor_at_solution
)


# ============================================================
# PART B — APPLY THE ANCHOR TO SCENARIO ACCOUNTING
# ============================================================
# Purpose:
# Use the Anchor result to perform actual in-year or scenario accounting.
# This is where HO fish, pHOS, pNOB, PNI, final brood allocation, and
# final mixed-origin spawner accounting are calculated.
# Part B priority order:
#   1. Protect broodstock
#   2. Allow fisheries from remaining available fish
#   3. Allow spawning with fish remaining after broodstock and fisheries
#   4. Remove excess captured HO
#
# The Anchor is used to calculate the pNOB line.
# It is not used as a fishery cap in Part B.

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

fishery_result <- calculate_scenario_fisheries(scenario_inputs, accounting_params, anchor_result)
weir_result <- calculate_scenario_weir_accounting(scenario_inputs, accounting_params, fishery_result)
pnob_slope_result <- calculate_pnob_slope(scenario_inputs, anchor_result)
brood_result <- allocate_final_brood(scenario_inputs, weir_result, pnob_slope_result)
final_metrics <- calculate_final_spawners(scenario_inputs, accounting_params, weir_result, brood_result)

accounting_result <- list(
  scenario_inputs = scenario_inputs,
  fishery_result = fishery_result,
  weir_result = weir_result,
  brood_result = brood_result,
  final_metrics = final_metrics
)


# ============================================================
# REVIEW TABLES — PRINT ONLY, NO CSV EXPORTS
# ============================================================

anchor_table <- tibble(
  Section = c(
    rep("Median parameter", 4),
    rep("Anchor input", 3),
    rep("Anchor result", 7)
  ),
  Name = c(
    "trap_prop_no",
    "weir_efficiency",
    "survival_above_weir",
    "survival_below_weir",
    "brood_need",
    "spawner_goal",
    "wl_scaling",
    "anchor_required_no",
    "sport_impact_wl_at_anchor",
    "no_captured_at_anchor",
    "no_spawners_at_anchor",
    "sport_no_impacts_at_anchor",
    "treaty_no_impacts_at_anchor",
    "anchor_feasibility"
  ),
  Value = c(
    accounting_params$trap_prop_no,
    accounting_params$weir_efficiency,
    accounting_params$survival_above_weir,
    accounting_params$survival_below_weir,
    anchor_result$brood_need,
    anchor_result$spawner_goal,
    anchor_result$wl_scaling,
    anchor_result$anchor_required_no,
    anchor_result$sport_impact_wl_at_anchor,
    anchor_result$no_captured_at_anchor,
    anchor_result$no_spawners_at_anchor,
    anchor_result$sport_no_impacts_at_anchor,
    anchor_result$treaty_no_impacts_at_anchor,
    anchor_feasibility(anchor_result$anchor_required_no, anchor_inputs, accounting_params)
  ),
  YearsUsed = c(
    accounting_params$source_trap_no,
    accounting_params$source_weir,
    accounting_params$source_abv,
    accounting_params$source_blw,
    "", "", "", "", "", "", "", "", "", ""
  )
) %>%
  mutate(Value = fmt(Value))

partB_rows <- list(
  tibble(
    Section = "Scenario input",
    Name = c("no_manarea_est", "ho_manarea_est", "brood_need", "spawner_goal",
             "wl_scaling", "utilization_no", "utilization_ho"),
    Value = c(
      scenario_inputs$no_manarea_est,
      scenario_inputs$ho_manarea_est,
      scenario_inputs$brood_need,
      scenario_inputs$spawner_goal,
      scenario_inputs$wl_scaling,
      scenario_inputs$utilization_no,
      scenario_inputs$utilization_ho
    )
  ),
  tibble(
    Section = "Accounting parameter",
    Name = c("trap_prop_no", "trap_prop_ho", "weir_efficiency", "survival_above_weir", "survival_below_weir"),
    Value = c(
      accounting_params$trap_prop_no,
      accounting_params$trap_prop_ho,
      accounting_params$weir_efficiency,
      accounting_params$survival_above_weir,
      accounting_params$survival_below_weir
    )
  ),
  tibble(
    Section = "Projected fishery",
    Name = c(
      "sport_closed",
      "no_wl",
      "sport_impact_wl_projected",
      "sport_no_impacts_projected",
      "treaty_no_impacts_projected",
      "total_no_impacts_projected",
      "ho_sport_harvest_projected",
      "ho_treaty_harvest_projected",
      "total_ho_harvest_projected"
    ),
    Value = c(
      as.numeric(fishery_result$sport_closed),
      fishery_result$no_wl,
      fishery_result$sport_impact_wl_projected,
      fishery_result$sport_no_impacts_projected,
      fishery_result$treaty_no_impacts_projected,
      fishery_result$total_no_impacts_projected,
      fishery_result$ho_sport_harvest_projected,
      fishery_result$ho_treaty_harvest_projected,
      fishery_result$total_ho_harvest_projected
    )
  ),
  tibble(
    Section = "Brood-protection caps",
    Name = c(
      "no_brood_target",
      "required_no_after_fisheries",
      "max_no_impacts_allowed",
      "ho_brood_needed",
      "required_ho_after_fisheries",
      "max_ho_harvest_allowed"
    ),
    Value = c(
      fishery_result$no_brood_target,
      fishery_result$required_no_after_fisheries,
      fishery_result$max_no_impacts_allowed,
      fishery_result$ho_brood_needed,
      fishery_result$required_ho_after_fisheries,
      fishery_result$max_ho_harvest_allowed
    )
  ),
  tibble(
    Section = "Allowed fishery after brood caps",
    Name = c(
      "sport_no_impacts",
      "treaty_no_impacts",
      "total_no_impacts_allowed",
      "ho_sport_harvest",
      "ho_treaty_harvest",
      "total_ho_harvest_allowed"
    ),
    Value = c(
      fishery_result$sport_no_impacts,
      fishery_result$treaty_no_impacts,
      fishery_result$total_no_impacts_allowed,
      fishery_result$ho_sport_harvest,
      fishery_result$ho_treaty_harvest,
      fishery_result$total_ho_harvest_allowed
    )
  ),
  tibble(
    Section = "Post-fishery abundance",
    Name = c("no_after_fisheries", "ho_after_fisheries"),
    Value = c(weir_result$no_after_fisheries, weir_result$ho_after_fisheries)
  ),
  tibble(
    Section = "Weir partition",
    Name = c(
      "no_captured_weir", "ho_captured_weir",
      "no_above_weir_uncaptured", "ho_above_weir_uncaptured",
      "no_below_weir", "ho_below_weir"
    ),
    Value = c(
      weir_result$no_captured_weir,
      weir_result$ho_captured_weir,
      weir_result$no_above_weir_uncaptured,
      weir_result$ho_above_weir_uncaptured,
      weir_result$no_below_weir,
      weir_result$ho_below_weir
    )
  ),
  tibble(
    Section = "pNOB and brood",
    Name = c(
      "pnob_slope", "no_brood_desired", "no_brood_actual", "ho_brood_actual",
      "no_captured_available_above", "ho_captured_available_above"
    ),
    Value = c(
      brood_result$pnob_slope,
      brood_result$no_brood_desired,
      brood_result$no_brood_actual,
      brood_result$ho_brood_actual,
      brood_result$no_captured_available_above,
      brood_result$ho_captured_available_above
    )
  ),
  tibble(
    Section = "Spawner accounting",
    Name = c(
      "no_spawners_above",
      "ho_spawners_above_uncapt",
      "no_spawners_below",
      "ho_spawners_below",
      "current_spawners_without_HO_captured_above"
    ),
    Value = c(
      final_metrics$no_spawners_above,
      final_metrics$ho_spawners_above_uncapt,
      final_metrics$no_spawners_below,
      final_metrics$ho_spawners_below,
      final_metrics$current_spawners_without_HO_captured_above
    )
  ),
  tibble(
    Section = "HO placement",
    Name = c(
      "spawner_deficit",
      "ho_needed_pre_survival",
      "ho_captured_above",
      "ho_captured_removed"
    ),
    Value = c(
      final_metrics$spawner_deficit,
      final_metrics$ho_needed_pre_survival,
      final_metrics$ho_captured_above,
      final_metrics$ho_captured_removed
    )
  ),
  tibble(
    Section = "Final totals and indices",
    Name = c("no_spawners_total", "ho_spawners_total", "system_spawners_total", "phos", "pnob", "pni"),
    Value = c(
      final_metrics$no_spawners_total,
      final_metrics$ho_spawners_total,
      final_metrics$system_spawners_total,
      final_metrics$phos,
      final_metrics$pnob,
      final_metrics$pni
    )
  )
)

partB_table <- bind_rows(partB_rows) %>%
  mutate(
    YearsUsed = case_when(
      Name == "trap_prop_no"        ~ accounting_params$source_trap_no,
      Name == "trap_prop_ho"        ~ accounting_params$source_trap_ho,
      Name == "weir_efficiency"     ~ accounting_params$source_weir,
      Name == "survival_above_weir" ~ accounting_params$source_abv,
      Name == "survival_below_weir" ~ accounting_params$source_blw,
      TRUE ~ ""
    ),
    Value = fmt(Value)
  )

final_dispositions_table <- tibble(
  Section = c(
    rep("Brood disposition", 4),
    rep("HO captured disposition", 2),
    rep("Spawning disposition", 3),
    rep("Final indices", 3)
  ),
  Name = c(
    "no_brood_actual",
    "ho_brood_actual",
    "no_captured_available_above",
    "ho_captured_available_above",
    "ho_captured_above",
    "ho_captured_removed",
    "no_spawners_total",
    "ho_spawners_total",
    "system_spawners_total",
    "phos",
    "pnob",
    "pni"
  ),
  Value = c(
    brood_result$no_brood_actual,
    brood_result$ho_brood_actual,
    brood_result$no_captured_available_above,
    brood_result$ho_captured_available_above,
    final_metrics$ho_captured_above,
    final_metrics$ho_captured_removed,
    final_metrics$no_spawners_total,
    final_metrics$ho_spawners_total,
    final_metrics$system_spawners_total,
    final_metrics$phos,
    final_metrics$pnob,
    final_metrics$pni
  ),
  YearsUsed = ""
) %>%
  mutate(Value = fmt(Value))

if (print_tables) {
  cat("\n========================================\n")
  cat("Lostine Toolkit — Reviewer Refactor\n")
  cat("========================================\n\n")

  cat("# Table A — Anchor Method\n")
  print(kable(anchor_table, align = "llcr"))

  cat("\n# Table B — Scenario Accounting\n")
  print(kable(partB_table, align = "llcr"))

  cat("\n# Table C — Final Dispositions and Metrics\n")
  print(kable(final_dispositions_table, align = "llcr"))
}


# ============================================================
# TABLE D — SENSITIVITY SWEEP ACROSS NO ABUNDANCE
# ============================================================
# Purpose:
# Show how key outputs change across a range of natural-origin
# abundance values, holding all other inputs at scenario defaults.
# This helps co-managers understand method behavior without
# having to run the tool repeatedly.
#
# HO abundance is scaled proportionally with NO across the sweep
# using the HO:NO ratio from the current scenario, so the sweep
# reflects a plausible range of conditions rather than holding
# HO fixed at an arbitrary value.

ho_no_ratio <- scenario_inputs$ho_manarea_est / scenario_inputs$no_manarea_est

no_sweep_values <- c(100, 200, 300, 400, 500, 600, 700, 800,
                     anchor_result$anchor_required_no,
                     1000, 1200, 1500, 2000)

sweep_rows <- purrr::map_dfr(sort(no_sweep_values), function(no_val) {
  si <- scenario_inputs
  si$no_manarea_est <- no_val
  si$ho_manarea_est <- round(no_val * ho_no_ratio)
  
  fr <- calculate_scenario_fisheries(si, accounting_params, anchor_result)
  wr <- calculate_scenario_weir_accounting(si, accounting_params, fr)
  ps <- calculate_pnob_slope(si, anchor_result)
  br <- allocate_final_brood(si, wr, ps)
  fm <- calculate_final_spawners(si, accounting_params, wr, br)
  
  tibble(
    NO_est                  = round(no_val),
    HO_est                  = si$ho_manarea_est,
    pnob_slope              = fmt(ps),
    sport_NO_impact         = fmt(fr$sport_no_impacts),
    treaty_NO_impact        = fmt(fr$treaty_no_impacts),
    HO_sport_harvest        = fmt(fr$ho_sport_harvest),
    HO_treaty_harvest       = fmt(fr$ho_treaty_harvest),
    NO_brood                = fmt(br$no_brood_actual),
    HO_brood                = fmt(br$ho_brood_actual),
    HO_removed              = fmt(fm$ho_captured_removed),
    NO_spawners             = fmt(fm$no_spawners_total),
    HO_spawners             = fmt(fm$ho_spawners_total),
    system_spawners         = fmt(fm$system_spawners_total),
    pHOS                    = fmt(fm$phos),
    pNOB                    = fmt(fm$pnob),
    PNI                     = fmt(fm$pni),
    anchor_flag             = ifelse(
      abs(no_val - anchor_result$anchor_required_no) < 1,
      "<-- Anchor", ""
    )
  )
})

if (print_tables) {
  cat("\n# Table D — Sensitivity Sweep (NO abundance range, HO scaled proportionally)\n")
  cat(sprintf(
    "# HO:NO ratio held at %.2f  |  All other inputs at scenario defaults\n",
    ho_no_ratio
  ))
  print(kable(sweep_rows, align = "rrrrrrrrrrrrrrrrr"))
}