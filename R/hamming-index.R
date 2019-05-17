#' Compute hamming index with a sliding window
#'
#' Takes a sequence of asleep/awake states in logical format (TRUE-> asleep)
#' and computes the distance to the closes fully consolidated sequence of equal
#' length and proportion of asleep/awake states
#' @param asleep_sequence A logical vector stating the sequence of asleep/awake states
#' @details The hamming index of a sleep run i the minimal number pairs of positions
#' that need to be swapped to get to a continuous run of asleep states i.e.
#' a run with only one sleeping bout
#' @return A double (numeric) which results from taking the hamming index
#' and dividing by the length of the sequence
#' @export
hamming_index <- function(asleep_sequence) {

  # count the number of bins where the fly is asleep
  x <- sum(asleep_sequence)
  if(x == 0) {
    return(as.double(NA))
  }

  # count how many bins we have
  n <- length(asleep_sequence)

  # initialize a counter storing the max number of asleep
  # states found within all the windows of length x analyzed
  # so far
  max_ones <- 0
  # initialize a counter storing
  # the number of ones of the current window
  nones <- 0

  # number of windows to be analyzed
  niters <- (n - x + 1)

  # for every window
  for (i in 1:niters) {
    # get the number of bins where the fly was asleep
    nones <- sum(asleep_sequence[i:(i+x-1)])
    # check if this number is > than the recorded max so far
    found <- nones > max_ones

    if(is.na(found)) {
      print(i)
      print(i+x-1)
      print(asleep_sequence)
    }
    # if it is >, set the max to this number
    if (found) {
      max_ones <- nones
    }
  }
  # the hamming index will be the difference between
  # the total number of asleep bins (x)
  # and the maximum number of asleep bins found within the windows
  # (max ones)
  # normalized by the number of bins
  # i.e. the numerator is the number of awake states in the window
  # with the greatest number of asleep states
  hi <- (x - max_ones) / n
  return(hi)
}

