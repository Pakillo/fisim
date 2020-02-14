#' Simple random sampling estimator of population means
#'
#' @param data A data.table object with plot-level information of stand
#'   characteristics as coming from the \code{sum_data} function.
#' @param est_col A character vector of column names for which estimates are
#'   produced. Used when there is more than one way to expand observations from
#'   individual sample locations.
#' @details Variance estimation is done following the simple random sampling
#'   estimator for replicated sampling.
#' @return A data.table object with mean and variance estimates for the target
#'  variables.
#' @import data.table
#' @export
est_srs <- function(point_data, est_col = NULL) {
  if (!is.point_data(point_data)) {
    stop("'point_data' is not an object of class point_data!")
  }
  if (is.null(est_col)) {
    if (point_data$r_design == "k_tree") {
      est_col <- c("Prodan", "Eberhardt", "d_mean", "a_mean");
    } else if (point_data$r_design != "k_tree") {
      est_col <- "ef_ana";
    }
  }

  dt_est <- point_data$data[,
                            list(est_appr = names(.SD),
                                 y_hat = sapply(.SD, mean),
                                 v_hat = sapply(.SD, stats::var)/.N),
                            by = list(id_sample, variable),
                            .SDcols = est_col];
  return(dt_est);
}
