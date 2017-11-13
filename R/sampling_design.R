#' Sample locations in a study area of rectangular shape
#'
#' @param x_range A numeric vector holding minimum and maximum coordinate in
#'   x-direction.
#' @param y_range A numeric vector holding minimum and maximum coordinate in
#'   y-direction.
#' @param n Sample size, i.e. number of sample locations.
#' @param M Number of independent samples of size \code{n}.
#'
#' @return A \code{\link[data.table]{data.table}} object with \code{n} rows
#'   holding an identifier and xy-coordinates.
#' @export
#'
xy_sample <- function(x_range, y_range, n, M = 1) {
  x_u <- runif(n*M);
  y_u <- runif(n*M);
  dt_s <- data.table(id_set = rep(1:M, each = n),
                     id_point = rep(1:n, M),
                     x_pt = min(x_range) + x_u*max(x_range),
                     y_pt = min(y_range) + y_u*max(y_range));
  return(dt_s);
}
