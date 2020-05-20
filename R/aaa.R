#' @importFrom data.table ":="
#' @importFrom data.table "key"
#' @importFrom data.table "%between%"
#' @import fslbehavr
NULL


#' Score sleep behaviour from immobility
#'
#' This function first uses a motion classifier to decide whether an animal is moving during a given time window.
#' Then, it defines sleep as contiguous immobility for a minimum duration.
#'
#' @param data  [data.table] containing behavioural variable from or one multiple animals.
#' When it has a key, unique values, are assumed to represent unique individuals (e.g. in a [behavr] table).
#' Otherwise, it analysis the data as coming from a single animal. `data` must have a column `t` representing time.
#' @param time_window_length number of seconds to be used by the motion classifier.
#' This corresponds to the sampling period of the output data.
#' @param min_time_immobile Minimal duration (in s) of a sleep bout.
#' Immobility bouts longer or equal to this value are considered as sleep.
#' @param motion_detector_FUN function used to classify movement
#' @param ... extra arguments to be passed to `motion_classifier_FUN`.
#' @return a [behavr] table similar to `data` with additional variables/annotations (i.e. `moving` and `asleep`).
#' The resulting data will only have one data point every `time_window_length` seconds.
#' @details
#' The default `time_window_length` is 300 seconds -- it is also known as the "5-minute rule".
#' `sleep_annotation` is typically used for ethoscope data, whilst `sleep_dam_annotation` only works on DAM2 data.
#' These functions are *rarely used directly*, but rather passed as an argument to a data loading function,
#' so that analysis can be performed on the go.
#' @examples
# # We start by making toy data for one animal:
#' dt_one_animal <- toy_ethoscope_data(seed=2)
#' ####### Ethoscope, corrected velocity classification #########
#' sleep_dt <-  sleep_annotation(dt_one_animal, masking_duration=0)
#' print(sleep_dt)
#' # We could make a sleep `barecode'
#' \dontrun{
#' library(ggplot2)
#' ggplot(sleep_dt, aes(t,y="Animal 1",fill=asleep)) +
#'                                    geom_tile() + scale_x_time()
#' }
#' ####### Ethoscope, virutal beam cross classification #########
#' sleep_dt2 <-  sleep_annotation(dt_one_animal,
#'                              motion_detector_FUN=virtual_beam_cross_detector)
#' \dontrun{
#' library(ggplot2)
#' ggplot(sleep_dt2, aes(t,y="Animal 1",fill=asleep)) +
#'                                    geom_tile() + scale_x_time()
#' }
#' ####### DAM data, de facto beam cross classification ######
#' dt_one_animal <- toy_dam_data(seed=7)
#' sleep_dt <- sleep_dam_annotation(dt_one_animal)
#' \dontrun{
#' library(ggplot2)
#' ggplot(sleep_dt, aes(t,y="Animal 1",fill=asleep)) +
#'                                    geom_tile() + scale_x_time()
#' }
#' @seealso
#' * [motion_detectors] -- options for the `motion_detector_FUN` argument
#' * [bout_analysis] -- to further analyse sleep bouts in terms of onset and length
#' @references
#' * The relevant [rethomic tutorial section](https://rethomics.github.io/sleepr) -- on sleep analysis
#' @import logging
#' @export
sleep_annotation_wrapper <- function(time_window_length.=10, min_time_immobile.=300, velocity_correction_coef.=0.003, motion_detector_FUN. = max_velocity_detector, verbose. = TRUE, ...) {

  sleep_annotation <- function(data,
                               time_window_length = time_window_length.,
                               min_time_immobile = min_time_immobile.,
                               velocity_correction_coef = velocity_correction_coef.,
                               motion_detector_FUN = motion_detector_FUN.,
                               verbose = verbose,
                               ...
  ){

    moving = .N = is_interpolated  = .SD = asleep = NULL
    # all columns likely to be needed.

    if (verbose) {
      logging::loginfo('Enclosed annotation evironment')
      logging::loginfo(glue::glue('velocity_correction_coef: {velocity_correction_coef}'))
      logging::loginfo(glue::glue('min_time_immobile: {min_time_immobile}'))
      logging::loginfo(glue::glue('time_window_length: {time_window_length}'))
    }

    columns_to_keep <- c("t", "x", "y", "max_velocity", "interactions",
                         "beam_crosses", "moving","asleep", "is_interpolated")

    data_copy <- copy(data)

    data_copy <- data_copy[, phase := ifelse(t %% hours(24) > hours(12), 'D', 'L')]

    wrapped <- function(d, ...) {

      if(nrow(d) < 100)
        return(NULL)
      # todo if t not unique, stop

      logging::loginfo("Running motion_detector_FUN")
      d_small <- motion_detector_FUN(d, time_window_length., ...)

      if(key(d_small) != "t")
        stop("Key in output of motion_classifier_FUN MUST be `t'")

      if(nrow(d_small) < 1)
        return(NULL)

      logging::loginfo("Computing time_map")
      # the times to  be queried
      time_map <- data.table::data.table(t = seq(from=d_small[1,t], to=d_small[.N,t], by=time_window_length),
                                         key = "t")
      missing_val <- time_map[!d_small]

      d_small <- d_small[time_map,roll=T]
      d_small[,is_interpolated := FALSE]
      d_small[missing_val,is_interpolated:=TRUE]
      d_small[is_interpolated == T, moving := FALSE]

      logging::loginfo("Computing sleep")
      d_small[,asleep := sleep_contiguous(moving,
                                          1/time_window_length,
                                          min_valid_time = min_time_immobile)]

      logging::loginfo("Removing missing data")
      d_small <- stats::na.omit(d[d_small,
                                  on=c("t"),
                                  roll=T])

      logging::loginfo("Subsetting result")
      d_small[, intersect(columns_to_keep, colnames(d_small)), with=FALSE]

      return(d_small)
    }

    if(is.null(key(data_copy))) {
      return(wrapped(data_copy))
    }

    data_copy <- data_copy[,
         wrapped(.SD, ...),
         by=key(data_copy)]

    return(data_copy)
  }

  attr(sleep_annotation, "needed_columns") <- function(motion_detector_FUN = max_velocity_detector,
                                                       ...){
    needed_columns <- attr(motion_detector_FUN, "needed_columns")
    if(!is.null(needed_columns))
      needed_columns(...)
  }

  return(sleep_annotation)

}


#' @export
#' @import logging
# @rdname
# TODO Make sure this function is ok with receiving sleep_annotation (ethoscope)
# specific arguments even though they are not used
sleep_dam_annotation_wrapper <- function(min_time_immobile.=300, verbose.=TRUE) {

  sleep_dam_annotation <- function(data, min_time_immobile = min_time_immobile., verbse = verbose.) {

    asleep = moving = activity = duration = .SD = . = NULL

    if (verbose) {
      logging::loginfo('Enclosed annotation evironment')
      logging::loginfo(glue::glue('velocity_correction_coef: {velocity_correction_coef}'))
      logging::loginfo(glue::glue('min_time_immobile: {min_time_immobile}'))
      logging::loginfo(glue::glue('time_window_length: {time_window_length}'))
    }


    data_copy <- copy(data)
    data_copy <- data_copy[, phase := ifelse(t %% hours(24) > hours(12), 'D', 'L')]

    wrapped <- function(d){
      if(! all(c("activity", "t") %in% names(d)))
        stop("data from DAM should have a column named `activity` and one named `t`")

      out <- data.table::copy(d)
      col_order <- c(colnames(d),"moving", "asleep")
      out[, moving := activity > 0]
      bdt <- bout_analysis(moving, out)
      bdt[, asleep := duration >= min_time_immobile & !moving]
      out <- bdt[,.(t, asleep)][out, on = "t", roll=TRUE]
      data.table::setcolorder(out, col_order)
      return(out)
    }

    if(is.null(key(data_copy)))
      return(wrapped(data_copy))

    data_copy <- data_copy[,
         wrapped(.SD),
         by=key(data_copy)]

    metadata <- data_copy[,meta=T]
    metadata$machine_name <- metadata$file_info %>% lapply(function(x) stringr::str_split(string = x$file, pattern = "\\.")[[1]][1]) %>% unlist
    setmeta(data_copy, metadata)
    return(data_copy)

  }

  return(sleep_dam_annotation)
}
