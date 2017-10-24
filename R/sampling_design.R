#' Sample locations in a study area of rectangular shape
#'
#' @param x_range A numeric vector holding minimum and maximum coordinate in
#'   x-direction.
#' @param y_range A numeric vector holding minimum and maximum coordinate in
#'   y-direction.
#' @param n Sample size, i.e. number of sample locations.
#'
#' @return A \code{\link[data.table]{data.table}} object with \code{n} rows
#'   holding an identifier and xy-coordinates.
#' @export
#'
xy_sample <- function(x_range, y_range, n) {
  u_x <- runif(n);
  u_y <- runif(n);
  s <- data.table(id = 1:n,
                  x_sl = min(x_range) + u_x*max(x_range),
                  y_sl = min(y_range) + u_y*max(y_range));
  return(s);
}