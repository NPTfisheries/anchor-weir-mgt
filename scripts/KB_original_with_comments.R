##############################################
# Lostine Management Area Accounting Toolkit #
# Reviewer Refactor: Anchor + Accounting      #
##############################################

# ============================================================
# 0) CONFIG, LIBRARIES, AND USER OPTIONS
# ============================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(knitr)
})

source('./R/helpers.R')
source('./R/anchor_funs.R')
source('./R/harvest_funs.R')
source('./R/weir_funs.R')


# Display rounding only. Calculations retain full precision.
disp_digits <- 2

# Historical input file used to calculate median-derived parameters.
input_file <- "input_chs.csv"
file_path <- file.path('./data',input_file)
# Optional: default scenario values from a CSV year row.
# For reviewer workflow, leave FALSE and edit the manual values below.
use_csv_scenario_values <- FALSE
csv_scenario_year <- 2023

# Print reviewer tables in R.
print_tables <- TRUE

# Utility: 2-digit year string for table display.

message(sprintf("Lostine Toolkit — Historical file: %s", input_file))

# ============================================================
# 1) HISTORICAL DATA AND MEDIAN-DERIVED PARAMETERS
# ============================================================
# The historical CSV is used here only to calculate median-derived
# accounting parameters. Scenario/in-year values are entered later.

dat <- read_csv(file_path, show_col_types = FALSE) %>%
  filter(.data$Stock == "LST") %>%
  arrange(desc(.data$Year))


# Build trap ratio columns safely. Division by zero becomes NA.
dat <- dat %>%
  mutate(
    trap_ratio_no = if_else(Post.NO == 0, NA_real_, Trap.NO / Post.NO), # trap_ratio is the number trapped out of total return?
    trap_ratio_ho = if_else(Post.HO == 0, NA_real_, Trap.HO / Post.HO)
  )

m_trap_no <- first5_non_na(dat$trap_ratio_no, dat$Year)
m_trap_ho <- first5_non_na(dat$trap_ratio_ho, dat$Year)
m_weir    <- first5_non_na(dat$Weir.Eff,     dat$Year)
m_abv     <- first5_non_na(dat$PSS.ABV,      dat$Year) #what is PSS?
m_blw     <- first5_non_na(dat$PSS.BLW,      dat$Year)

# runs a check to ensure at least three years of data
purrr::walk2(
  list(m_trap_no, m_trap_ho, m_weir, m_abv, m_blw),
  c("trap_prop_no", "trap_prop_ho", "weir_efficiency", "survival_above_weir", "survival_below_weir"),
  check_median
)

median_params <- list(
  trap_prop_no        = median(m_trap_no$Value), # median proportion of NO captured at the weir; is this equal to weir efficiency?
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
# PART A — SET THE ANCHOR
# ============================================================
# Purpose:
# Determine the NO-only Anchor abundance: the minimum natural-origin
# abundance required to satisfy both broodstock need and spawner goal.
#
# This is a standalone methodological section. It does not calculate HO
# values, pHOS, pNOB, PNI, or final mixed-origin accounting.


#debugonce(find_smallest_feasible_anchor)

anchor_required_no <- find_smallest_feasible_anchor(
  anchor_inputs = anchor_inputs,
  accounting_params = accounting_params,
  start_est = 0)

# ---------
# RK added this function to better understand the sensitivity of the spawner goal
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

test_spawner_goal <- seq(100, 2000, by = 50)

test_anchor_results <- evaluate_anchor_spawner_goal(
  spawner_goals = test_spawner_goal,
  anchor_inputs = anchor_inputs,
  accounting_params = accounting_params
)

# Fit linear model
fit <- lm(anchor ~ spawner_goal, data = test_anchor_results)

coefs <- coef(fit)
r2 <- summary(fit)$r.squared

eqn <- sprintf(
  "Anchor = %.2f + %.3f × Spawner Goal\nR² = %.3f",
  coefs[1], coefs[2], r2
)

ggplot(test_anchor_results,
       aes(x = spawner_goal, y = anchor)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  geom_smooth(method = "lm",
              se = FALSE,
              linetype = 2,
              colour = "red") +
  annotate(
    "label",
    x = min(test_anchor_results$spawner_goal),
    y = max(test_anchor_results$anchor),
    hjust = 0,
    vjust = 1,
    label = eqn,
    size = 5
  ) +
  labs(
    x = "Spawner Goal",
    y = "Minimum Feasible NO Anchor",
    title = "Sensitivity of Anchor to Spawner Goal"
  ) +
  theme_bw()

#----------- back to Kyle's
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
