log10x1000_inv <- function(x) { return(10 ^ (x / 1000))}

#' A function to compute the distance traversed by an animal
#' Preprocess a raw ethoscope dataset by computing the sum of the number of pixels
#' traversed by an animal on each time bin
#' @return data.table of columns t and dist_sum
#' @inheritParams sleep_annotation
#' @export
#' @import data.table
distance_sum_enclosed <- function(data, time_window_length=10) {

  . <- xy_dist_log10x1000 <- NULL
  d <- prepare_data_for_motion_detector(data,
                                        c("t", "xy_dist_log10x1000"),
                                        time_window_length)
  d[, t := NULL]
  d <- d[, .(dist_sum = sum(log10x1000_inv(xy_dist_log10x1000))),  by = 't_round']
  return(d)
}

attr(distance_sum_enclosed, "needed_columns") <- function() {
  c("t", "dist_sum")
}


#' Compute velocity aggregates using xy_dist_log10x1000
#' @inheritParams movement_detector_enclosed
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
attr(velocity_avg_enclosed, "needed_columns") <- function() {
  c("t", "vel_avg")
}

#' Generic function to aggregate movement with some statistic
#' @param data  [data.table] containing behavioural variable from or one multiple animals.
#' When it has a key, unique values, are assumed to represent unique individuals (e.g. in a [behavr] table).
#' Otherwise, it analysis the data as coming from a single animal. `data` must have a column `t` representing time.
#' @param time_window_length number of seconds to be used by the motion classifier.
#' This corresponds to the sampling period of the output data.
#' @param func Aggregating function (max, min, median, mean, etc)
#' @param feature Name of a column in the sqlite3 file e.g. xy_dist_log10x1000
#' @param statistic Name of the column resulting from aggregation e.g. max_movement
#' @param score Name of the column providing a score i.e. category to the statistic e.g. micromovement
#' score is usually a binary variable i.e. TRUE/FALSE
#' @param preproc_FUN Optional, function to preprocess the input before computing the feature
#' (if the data needs some transformation like reverting xy_dist_log10x1000 back to a distance)
#' @param time_window_length Size of non overlapping time bins, in seconds
# @param threshold If the statistic is greater than this value, the score is TRUE, and 0 otherwise
#' @rdname custom_annotation_wrapper
#' @details movement_detector_enclosed takes:
#' \itemize{
#' \item{the name of an R summary function (mean, max, etc)}
#' \item{the name of a column in the future datasets to apply the function to}
#' \item{the name of the resulting summary column}
#' \item{the name of an alternative boolean column, which is set to TRUE if the summary column has a value greater than a threshold (default 1)}
#' \item{a preprocessing function to be applied to the column before the summary function is applied to it}
#' }
#' @example
#' max_movement_detector <- custom_annotation_wrapper(movement_detector_enclosed("max", "xy_dist_log10x1000", "max_movement", "micromovement", log10x1000_inv))
movement_detector_enclosed <- function(func, feature, statistic, score, preproc_FUN=NULL) {

  dt <- . <- NULL

  closure <- function(data, time_window_length=10, threshold=1) {

    message(paste0("Movement detector - ", func, " running.\ntime_window_length = ", time_window_length))
    func <- match.fun(func)
    # data$body_movement <- data$xy_dist_log10x1000
    d <- prepare_data_for_motion_detector(data,
                                          c("t", feature, "x"),
                                          time_window_length,
                                          "has_interacted")

    d[, dt := c(NA, diff(t))]

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

    # Set t_round as the representative time of the window
    # i.e. t becomes the begining of the window and not the t
    # of the first frame in the window
    return(d_small)
  }

  attr(closure, "needed_columns") <- function() {
    c("t", statistic, score)
  }
  attr(closure, "parameters") <- function() {
    return(names(formals(func)))
  }

  attr(closure, "variables") <- function() {
    statistic
  }

  return(closure)
}


#' Custom annotation from the dt_raw file
#'
#' This function gives aggregates a variable of interest in a custom way
#' All datapoints in every time_window_length seconds is aggregated into a single datapoint
#'
#' @param custom_function function used to produce the custom annotation
#' @param ... Extra arguments to be passed to `custom_function`.
#' @return a [behavr] table similar to `data` with additional variables/annotations.
#' The resulting data will only have one data point every `time_window_length` seconds.
#' @details
#' The default `time_window_length` is 300 seconds -- it is also known as the "5-minute rule".
#' custom_annotation_wrapper simplifies writing new annotation functions by leaving the shared functionality here
#' and the dedicated functionality to the new function.
#' This function adds to the functionality in the annotation function:
#' \itemize{
#' \item{Check a minimal amount of data is available and quit otherwise}
#' \item{Restore the name of the time column to remove the effects of binning}
#' \item{Check the amount of data after annotation is also enough (at least 1)}
#' \item{Apply a rolling interpolation of the labels to the data (assume the last available data point)}
#' }
#' It implements 3 attributes:
#' \itemize{
#' \item{needed_columns: A function that returns the columns needed by the function in its passed data}
#' \item{parameters: A function that returns the name of the parameters used by the function (including other functions' called by it)}
#' \item{variables: A function that returns the name of the newly produced columns by the function}
#' }
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
  attr(custom_annotation, "parameters") <- function() {
    args <- names(formals(custom_annotation))
    args <- c(args, attr(custom_function, "parameters")())
    args <- unique(args)
    args <- args[args != "..."]
    args <- args[args != "data"]
    return(args)
  }

  attr(custom_annotation, "variables") <- function() {
    attr(custom_function, "variables")()
  }

  return(custom_annotation)
}

#' @export
#' @rdname velocity_avg_enclosed
velocity_avg <- function(data, time_window_length)  {}
velocity_avg <- custom_annotation_wrapper(velocity_avg_enclosed)

#' Find the maximum distance traversed by the animal
#' @rdname max_movement_detector
#' @inheritParams sleep_annotation
#' @param threshold numeric, a value that splits a continuous variable into two states
#' @export
max_movement_detector <- function(data, time_window_length=10, threshold=1)  {}
max_movement_detector <- custom_annotation_wrapper(movement_detector_enclosed("max", "xy_dist_log10x1000", "max_movement", "micromovement", log10x1000_inv))

#' Find the median distance traversed by the animal
#' @rdname median_movement_detector
#' @inheritParams max_movement_detector
#' @export
median_movement_detector <- function(data, time_window_length=10, threshold=1) {}
median_movement_detector <- custom_annotation_wrapper(movement_detector_enclosed("median", "xy_dist_log10x1000", "median_movement", "micromovement", log10x1000_inv))

#' Find the total distance traversed by the animal
#' @rdname sum_movement_detector
#' @export
#' @inheritParams max_movement_detector
sum_movement_detector <- function(data, time_window_length=10, threshold=1) {}
sum_movement_detector <- custom_annotation_wrapper(movement_detector_enclosed("sum", "xy_dist_log10x1000", "sum_movement", "micromovement", log10x1000_inv))

#' @export
#' @inheritParams sleep_annotation
distance_annotation <- function(data, time_window_length=10) {}
distance_annotation <- custom_annotation_wrapper(sleepr::distance_sum_enclosed)

