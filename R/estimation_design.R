#' Simple random sampling estimator of population means
#'
#' @param data A data.table object with plot-level information of stand
#'   characteristics as comin from the \code{plot_data} function.
#'
#' @return A data.table object with mean and variance estimates of variables of
#'   interest.
#' @import data.table
#' @export
est_srs <- function(data) {
  dt_est <- data[,
                 list(variable = names(.SD),
                      y_hat = sapply(.SD, mean),
                      v_hat = sapply(.SD, function(x) var(x)/length(x))),
                 by = id_set,
                 .SDcols = 3:ncol(data)]

  return(dt_est);
}
