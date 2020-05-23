#' @export
#' @import data.table
distance_sum <- function(d, time_window_length) {
  d <- prepare_data_for_motion_detector(d,
                                        c("t", "xy_dist_log10x1000"),
                                        time_window_length)
  d[, t := NULL]
  d <- d[, .(dist_sum = sum(10**(xy_dist_log10x1000/1000))),  by = 't_round']
  d
}

attr(distance_sum, "needed_columns") <- function(...) {
  c("t", "dist_sum")
}
#' @export
#' @import data.table
velocity_avg <- function(d, time_window_length) {

  d2 <- copy(d)
  d2 <- prepare_data_for_motion_detector(d2,
                                        c("t", "xy_dist_log10x1000"),
                                        time_window_length)

  d2[, dt := c(NA, diff(t))]
  d2[, dist := 10**((xy_dist_log10x1000)/1e3)]
  d2[, velocity :=  dist / dt]
  # t_round is included in the columns because it is in the by arg
  d2 <- d2[, .(vel_avg = mean(velocity)),  by = 't_round']
  d2
}
attr(velocity_avg, "needed_columns") <- function(...) {
  c("t", "vel_avg")
}
#' Custom annotation from the dt_raw file
#'
#' This function gives aggregates a variable of interest in a custom way
#' All datapoints in every time_window_length seconds is aggregated into a single datapoint
#'
#' @param data  [data.table] containing behavioural variable from or one multiple animals.
#' When it has a key, unique values, are assumed to represent unique individuals (e.g. in a [behavr] table).
#' Otherwise, it analysis the data as coming from a single animal. `data` must have a column `t` representing time.
#' @param time_window_length number of seconds to be used by the motion classifier.
#' This corresponds to the sampling period of the output data.
#' @param custom_function function used to produce the custom annotation
#' @param ... extra arguments to be passed to `custom_function`.
#' @return a [behavr] table similar to `data` with additional variables/annotations.
#' The resulting data will only have one data point every `time_window_length` seconds.
#' @details
#' The default `time_window_length` is 300 seconds -- it is also known as the "5-minute rule".
#' @seealso
#' @import logging
#' @export
custom_annotation_closure <- function(custom_function) {


  custom_annotation <- function(data,
                                time_window_length = 10, #s
                                ...
  ){
    moving = .N = is_interpolated  = .SD = asleep = NULL
    # all columns likely to be needed.
    columns_to_keep <- c("t", attr(custom_function, 'needed_columns')())


    wrapped <- function(d){
      if(nrow(d) < 100)
        return(NULL)
      # todo if t not unique, stop

      d_small <- custom_function(d, time_window_length, ...)
      data.table::setnames(d_small, "t_round", "t")

      if(key(d_small) != "t")
        stop("Key in output of motion_classifier_FUN MUST be `t'")

      if(nrow(d_small) < 1)
        return(NULL)
      # the times to  be queried
      time_map <- data.table::data.table(t = seq(from=d_small[1,t], to=d_small[.N,t], by=time_window_length),
                                         key = "t")
      missing_val <- time_map[!d_small]

      d_small <- d_small[time_map,roll=T]
      d_small[,is_interpolated := FALSE]
      d_small[missing_val,is_interpolated:=TRUE]
      d_small <- stats::na.omit(d[d_small,
                                  on=c("t"),
                                  roll=T])
      d_small[, intersect(columns_to_keep, colnames(d_small)), with=FALSE]
    }

    if(is.null(key(data)))
      return(wrapped(data))
    data[,
         wrapped(.SD),
         by=key(data)]
  }

  attr(custom_annotation, "needed_columns") <- attr(custom_function, 'needed_columns')()

  return(custom_annotation)
}


# distance_sum_annotation <- custom_annotation_closure(custom_function = distance_sum)
# velocity_average_annotation <- custom_annotation_closure(custom_function = velocity_avg())

#
# attr(velocity_average_annotation, "needed_columns") <- c("t", "vel_avg")
