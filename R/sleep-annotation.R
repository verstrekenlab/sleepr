#' Compute euclidean distance between two points
#'
#' @param x numeric vector
#' @param y numeric vector
#' @export
euclidean_distance <- function(x, y) {
  square_diffs_x <- (x[-1]-x[1:(length(x)-1)])**2 # horizontal side of each triangle squared
  square_diffs_y <- (y[-1]-y[1:(length(y)-1)])**2 # vertical side of each triangle squared
  result <- c(NA, sqrt(square_diffs_x +  square_diffs_y)) # sqrt of the sum of squares
  return(result)
}



#' Score sleep behaviour from immobility
#'
#' This function first uses a motion classifier to decide whether an animal is moving during a given time window.
#' Then, it defines sleep as contiguous immobility for a minimum duration.
#'
#' @param data [data.table] containing behavioural variable from or one multiple animals.
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
sleep_annotation <- function() {
}

#' @export
sleep_dam_annotation <- function() {
  
}

sleep_dam_annotation <- sleep_dam_annotation_closure()

sleep_annotation <- sleep_annotation_closure()

