log10x1000_inv <- function(x) { return(10 ^ (x / 1000))}

#' A function to compute the distance traversed by an animal
#' Preprocess a raw ethoscope dataset by computing the sum of the number of pixels
#' traversed by an animal on each time bin
#' @return data.table of columns t and dist_sum
#' @export
#' @import data.table
distance_sum_enclosed <- function(d, time_window_length) {

  . <- xy_dist_log10x1000 <- NULL
  d <- prepare_data_for_motion_detector(d,
                                        c("t", "xy_dist_log10x1000"),
                                        time_window_length)
  d[, t := NULL]
  d <- d[, .(dist_sum = sum(10**(xy_dist_log10x1000/1000))),  by = 't_round']
  return(d)
}

attr(distance_sum_enclosed, "needed_columns") <- function(...) {
  c("t", "dist_sum")
}


#' Compute velocity aggregates using xy_dist_log10x1000
velocity_avg_enclosed <- function(data, time_window_length=10) {

  dt <- dist <- velocity <- vel_avg <- t <- dt <- xy_dist_log10x1000 <- . <- NULL

  data <- prepare_data_for_motion_detector(
    data,
    c("t", "xy_dist_log10x1000"),
    time_window_length
  )

  data[, dt := c(NA, diff(t))]
  data[, dist := 10**((xy_dist_log10x1000)/1e3)]
  data[, velocity :=  dist / dt]
  #t_round is included in the columns because it is in the by arg
  data <- data[, .(vel_avg = mean(velocity)), by = 't_round']
  return(data)
}
attr(velocity_avg_enclosed, "needed_columns") <- function(...) {
  c("t", "vel_avg")
}

#' Generic function to aggregate movement with some statistic
#' @param func Aggregating function (max, min, median, mean, etc)
#' @param feature Name of a column in the sqlite3 file e.g. body_movement
#' @param statistic Name of the column resulting from aggregation e.g. max_movement
#' @param score Name of the column providing a score i.e. category to the statistic e.g. micromovement
#' score is usually a binary variable i.e. TRUE/FALSE
#' @param preproc_FUN Optional, function to preprocess the input before computing the feature
#' (if the data needs some transformation like reverting xy_dist_log10x1000 back to a distance)
#' @param time_window_length Size of non overlapping time bins, in seconds
#' @param threshold If the statistic is greater than this value, the score is TRUE, and 0 otherwise
movement_detector_enclosed <- function(func, feature, statistic, score, preproc_FUN=NULL) {

  dt <- . <- NULL

  closure <- function(data, time_window_length=10, threshold=1) {

    # data$body_movement <- data$xy_dist_log10x1000
    d <- prepare_data_for_motion_detector(data,
                                          c("t", feature, "x"),
                                          time_window_length,
                                          "has_interacted")

    d[, dt := c(NA, diff(t))]
    #d[,surface_change := xor_dist * 1e-3]

    setnames(d, feature, "feature")
    # restore the distance from the log-transformed variable
    if (! is.null(preproc_FUN)) d[, feature := preproc_FUN(feature)]

    # Get a central summary value for variables of interest
    # for each window given by t_round
    # See prepare_data_for_motion_detector to learn
    # how is t_round computed
    # velocity_corrected -> max
    # has_interacted -> sum
    # beam_cross -> sum
    d_small <- d[, .(
      statistic = func(feature[1:.N])
    ), by = "t_round"]

    # Gist of the program!!
    # Score movement as TRUE/FALSE value for every window
    # Score is TRUE if max_velocity of the window is > 1
    # Score FALSE otherwise
    d_small[, score :=  ifelse(statistic > threshold, TRUE,FALSE)]
    setnames(d_small, "score", score)
    setnames(d_small, "statistic", statistic)
    # setnames(d_small, "feature", feature)

    # Set t_round as the representative time of the window
    # i.e. t becomes the begining of the window and not the t
    # of the first frame in the window
    return(d_small)
  }

  attr(closure, "needed_columns") <- function(...) {
    c("t", statistic, score)
  }

  return(closure)
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
#' @export
custom_annotation_wrapper <- function(custom_function) {

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
      d_small <- d_small[, intersect(columns_to_keep, colnames(d_small)), with=FALSE]
      return(d_small)
    }

    if(is.null(key(data)))
      return(wrapped(data))
    data[,
         wrapped(.SD),
         by=key(data)]
  }

  attr(custom_annotation, "needed_columns") <- function() {attr(custom_function, 'needed_columns')()}

  return(custom_annotation)
}

#' @export
#' @rdname velocity_avg_enclosed
velocity_avg <- function()  {}
velocity_avg <- custom_annotation_wrapper(velocity_avg_enclosed)

#' @export
#' @rdname movement_detector_enclosed
max_movement_detector <- function()  {}
max_movement_detector <- custom_annotation_wrapper(movement_detector_enclosed(function(x) x, max, "body_movement", "max_movement", "micromovement"))

#' @export
#' @rdname movement_detector_enclosed
median_movement_detector <- function() {}
median_movement_detector <- custom_annotation_wrapper(movement_detector_enclosed(function(x) x, median, "body_movement", "median_movement", "micromovement"))

#' @export
#' @rdname movement_detector_enclosed
sum_movement_detector <- function() {}
sum_movement_detector <- custom_annotation_wrapper(movement_detector_enclosed(function(x) x, sum, "body_movement", "sum_movement", "micromovement"))



