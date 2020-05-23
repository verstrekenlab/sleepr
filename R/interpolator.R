interpolate <- function(d_small, time_window_length) {
  # Generate a table with all the time windows
  # our data would contain
  # if there were no missing values
  time_map <- data.table::data.table(
  t = seq(from = d_small[1,t], to = d_small[.N,t], by = time_window_length),
  key = "t"
  )
  # merge time_map so only windows that are missing in d_small stay
  # i.e. subset time_map with key values NOT present in d_small
  missing_val <- time_map[!d_small]

  d_small <- d_small[time_map, roll = T]
  d_small[, is_interpolated := FALSE]
  d_small[missing_val, is_interpolated := TRUE]
  d_small[is_interpolated == T, moving := FALSE]
}
