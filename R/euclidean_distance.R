#' Compute euclidean distance between two points
#'
#' @param x numeric of length n giving X coordinates of a list of points
#' @param y numeric of same length as x giving Y coordinates of a list of points
#' @return A numeric of length n giving the distance between point i and i-1.
#' Its initial value is NA
#' @details Distance between 1,1 and 2,2 (should be square root of 2)
#' euclidean_distance(c(1,2), c(1,2))
#' NA 1.414214
#' @export
euclidean_distance <- function(x, y) {
  square_diffs_x <- (x[-1]-x[1:(length(x)-1)])**2 # horizontal side of each triangle squared
  square_diffs_y <- (y[-1]-y[1:(length(y)-1)])**2 # vertical side of each triangle squared
  result <- c(NA, sqrt(square_diffs_x +  square_diffs_y)) # sqrt of the sum of squares
  return(result)
}
