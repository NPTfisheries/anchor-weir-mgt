
fmt <- function(x) round(x, disp_digits)

y2 <- function(y) sprintf("%02d", as.integer(y) %% 100)

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