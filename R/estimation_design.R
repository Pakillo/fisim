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
  y_hat <- as.matrix(data[, lapply(.SD, mean), .SDcols = 2:ncol(data)]);
  v_hat <- as.matrix(data[, lapply(.SD, var), .SDcols = 2:ncol(data)]/nrow(data));
  return(data.table(var = colnames(y_hat),
                    y_hat = y_hat[1, ],
                    v_hat = v_hat[1, ]));
}