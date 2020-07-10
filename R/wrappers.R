#' #' Score sleep behaviour from immobility
#' #'
#' #' This function first uses a motion classifier to decide whether an animal is moving during a given time window.
#' #' Then, it defines sleep as contiguous immobility for a minimum duration.
#' #'
#' #' @param data  [data.table] containing behavioural variable from or one multiple animals.
#' #' When it has a key, unique values, are assumed to represent unique individuals (e.g. in a [behavr] table).
#' #' Otherwise, it analysis the data as coming from a single animal. `data` must have a column `t` representing time.
#' #' @param time_window_length number of seconds to be used by the motion classifier.
#' #' This corresponds to the sampling period of the output data.
#' #' @param min_time_immobile Minimal duration (in s) of a sleep bout.
#' #' Immobility bouts longer or equal to this value are considered as sleep.
#' #' @param motion_detector_FUN function used to classify movement
#' #' @param ... extra arguments to be passed to `motion_classifier_FUN`.
#' #' @return a [behavr] table similar to `data` with additional variables/annotations (i.e. `moving` and `asleep`).
#' #' The resulting data will only have one data point every `time_window_length` seconds.
#' #' @details
#' #' The default `time_window_length` is 300 seconds -- it is also known as the "5-minute rule".
#' #' `sleep_annotation` is typically used for ethoscope data, whilst `sleep_dam_annotation` only works on DAM2 data.
#' #' These functions are *rarely used directly*, but rather passed as an argument to a data loading function,
#' #' so that analysis can be performed on the go.
#' #' @examples
#' # # We start by making toy data for one animal:
#' #' dt_one_animal <- toy_ethoscope_data(seed=2)
#' #' ####### Ethoscope, corrected velocity classification #########
#' #' sleep_dt <-  sleep_annotation(dt_one_animal, masking_duration=0)
#' #' print(sleep_dt)
#' #' # We could make a sleep `barecode'
#' #' \dontrun{
#' #' library(ggplot2)
#' #' ggplot(sleep_dt, aes(t,y="Animal 1",fill=asleep)) +
#' #'                                    geom_tile() + scale_x_time()
#' #' }
#' #' ####### Ethoscope, virutal beam cross classification #########
#' #' sleep_dt2 <-  sleep_annotation(dt_one_animal,
#' #'                              motion_detector_FUN=virtual_beam_cross_detector)
#' #' \dontrun{
#' #' library(ggplot2)
#' #' ggplot(sleep_dt2, aes(t,y="Animal 1",fill=asleep)) +
#' #'                                    geom_tile() + scale_x_time()
#' #' }
#' #' ####### DAM data, de facto beam cross classification ######
#' #' dt_one_animal <- toy_dam_data(seed=7)
#' #' sleep_dt <- sleep_dam_annotation(dt_one_animal)
#' #' \dontrun{
#' #' library(ggplot2)
#' #' ggplot(sleep_dt, aes(t,y="Animal 1",fill=asleep)) +
#' #'                                    geom_tile() + scale_x_time()
#' #' }
#' #' @seealso
#' #' * [motion_detectors] -- options for the `motion_detector_FUN` argument
#' #' * [bout_analysis] -- to further analyse sleep bouts in terms of onset and length
#' #' @references
#' #' @export
#'
#' sleep_annotation_wrapper <- function(time_window_length = 10,
#'                                      min_time_immobile = 300,
#'                                      motion_detector_FUN = max_velocity_detector
#'                                      extra_columns = NULL,
#'                                      ...
#' ) {
#'
#'   sleep_annotation <- function(data,
#'                                time_window_length = time_window_length,
#'                                min_time_immobile = min_time_immobile,
#'                                motion_detector_FUN = motion_detector_FUN,
#'                                extra_columns = extra_columns,
#'                                ...
#'
#'   ){
#'
#'     moving = .N = is_interpolated  = .SD = asleep = NULL
#'     message(sprintf("Time window length: %f", time_window_length))
#'     message(sprintf("Minimum time immobile: %f", min_time_immobile))
#'
#'     # all columns likely to be needed.
#'     columns_to_keep <- c(
#'       extra_columns,
#'       c(
#'         "t", "x", "y", "max_velocity", "interactions",
#'         "beam_crosses", "moving","asleep", "is_interpolated"
#'       )
#'     )
#'
#'     wrapped <- function(d) {
#'       if(nrow(d) < 100)
#'         return(NULL)
#'       # todo if t not unique, stop
#'
#'       message("Running motion_detector_FUN")
#'       d_small <- motion_detector_FUN(d, time_window_length, ...)
#'
#'       if(key(d_small) != "t")
#'         stop("Key in output of motion_classifier_FUN MUST be `t'")
#'
#'       if(nrow(d_small) < 1)
#'         return(NULL)
#'
#'       # Set moving to FALSE for windows where no data is available
#'       d_small <- interpolate(d_small, time_window_length)
#'
#'       message("Computing sleep")
#'       # Actually annotate sleep
#'       # * Use moving as input
#'       # * Sampling frequency tells the algorithm
#'       # how frequently is moving sampled i.e.
#'       # how many seconds pass between each moving score
#'       # i.e. this frequency MUST be 1 / time_window_length
#'       # since time_window_length is the size of the windows
#'       # used when binning the time series in the moving annotation
#'       # * min_time_immobile tells the algorithm how many seconds must pass
#'       # for a run of non moving states to be considered as sleep
#'       # default of 300 s -> 5 minutes
#'       d_small[,asleep := sleep_contiguous(
#'         moving,
#'         1/time_window_length,
#'         min_valid_time = min_time_immobile
#'       )]
#'
#'       # TODO What is this doing exactly?
#'       message("Removing missing data")
#'       d_small <- stats::na.omit(
#'         d[d_small, on = c("t"), roll = T]
#'       )
#'
#'       message("Subsetting result")
#'       # keep the columns_to_keep only
#'       d_small[, intersect(columns_to_keep, colnames(d_small)), with=FALSE]
#'
#'       return(d_small)
#'     } # end of wrapped
#'
#'     # call wrapped on the whole dataframe at once
#'     # if there is no key
#'     if(is.null(key(data))) {
#'       return(wrapped(data))
#'     }
#'
#'     # however, if there is a key, call wrapped separately
#'     # for each block with the same key
#'     # (split-apply-combine)
#'     # this way we annotate for each fly separately
#'     data <- data[, wrapped(.SD), by = key(data)]
#'
#'     return(data)
#'   }
#'
#'   attr(sleep_annotation, "needed_columns") <- function(
#'     motion_detector_FUN = max_velocity_detector,
#'     ...
#'   ) {
#'     needed_columns <- attr(motion_detector_FUN, "needed_columns")
#'     if(!is.null(needed_columns)) needed_columns(...)
#'   }
#'
#'   return(sleep_annotation)
#' }
#'
#' #' @export
#' # @rdname
#'
#' sleep_dam_annotation <- function(data, min_time_immobile=300) {
#'
#'   asleep = moving = activity = duration = .SD = . = NULL
#'
#'   loginfo(sprintf('Minimum time immobile set to %s', min_time_immobile))
#'
#'   data_copy <- copy(data)
#'   data_copy <- data_copy[, phase := ifelse(t %% hours(24) > hours(12), 'D', 'L')]
#'
#'   wrapped <- function(d){
#'     if(! all(c("activity", "t") %in% names(d)))
#'       stop("data from DAM should have a column named `activity` and one named `t`")
#'
#'     out <- data.table::copy(d)
#'     col_order <- c(colnames(d),"moving", "asleep")
#'     out[, moving := activity > 0]
#'     bdt <- bout_analysis(moving, out)
#'     bdt[, asleep := duration >= min_time_immobile & !moving]
#'     out <- bdt[,.(t, asleep)][out, on = "t", roll=TRUE]
#'     data.table::setcolorder(out, col_order)
#'     return(out)
#'   }
#'
#'   if(is.null(key(data_copy)))
#'     return(wrapped(data_copy))
#'
#'   data_copy <- data_copy[,
#'                          wrapped(.SD),
#'                          by=key(data_copy)]
#'
#'   metadata <- data_copy[,meta=T]
#'   metadata$machine_name <- metadata$file_info %>% lapply(function(x) stringr::str_split(string = x$file, pattern = "\\.")[[1]][1]) %>% unlist
#'   setmeta(data_copy, metadata)
#'   return(data_copy)
#'
#' }
#'
#'
