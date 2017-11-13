#' Helper function to extract population elements into a sample using a list of
#' indices
#'
#' @param dt_pop The population from which to extract sampled elements.
#' @param s An integer matrix where rows represent sample locations and columns
#'   represent indices of the trees in \code{dt_pop} that are selected at the
#'   individual sample locations. Zeroes are used to indicate no neighbours and
#'   to ensure a rectangular data format.
#'
#' @return A \code{\link[data.table]{data.table}} object containing the population
#'   elements from data that were selected into the sample. The data.table
#'   object has a new variable \code{id_plot} for grouping single sample
#'   locations.
#' @export
#' @import data.table
extract_data <- function(dt_pop, s) {
  if (!is.data.table(dt_pop)) stop("dt_pop is not a data.table object!");

  dt_s_tree <- dt_pop[s[, "s"]];
  dt_s_tree[, ':='(id_set = s[, "id_set"],
                   id_point = s[, "id_point"],
                   ef = attributes(s)$ef)];

  # Make sure that all plot IDs are present in the output
  dt_plot_id <- data.table(id_point = 1:attributes(s)$sample_size);
  dt_s_tree <- dt_s_tree[dt_plot_id, on = "id_point"];

  setattr(dt_s_tree, "response_design", attributes(s)$response_design);
  return(dt_s_tree);
}


#' Helper function to summarize tree-level sample data to the plot-level for
#' estimation of poplation parameters
#'
#' @param dt_s_tree A data.table object as coming from \code{\link{sample_data}}
#'   function
#' @param target_var A character vector of variable names which should be
#'   summarized at the plot level
#'
#' @return A data.table object with stand characteristics at the plot level.
#'   Currently, basal area per hectare and the number of trees per hectare are
#'   included.
#' @export
#' @import data.table
sum_data <- function(dt_s_tree, target_vars) {

  if (!is.data.table(dt_s_tree)) stop("data is not a data.table object!");

  # Make sure that plots with no trees receive entries of zero
  for (c in c(target_vars, "f_edge", "ef")) {
    dt_s_tree[is.na(get(c)), (c) := 0];
  }
  dt_s_plot <- dt_s_tree[,
                         lapply(.SD, function(x) sum(ef*f_edge*x)),
                         keyby = c("id_set", "id_point"),
                         .SDcols = target_vars];
  return(dt_s_plot);
}


#' Helper function to estimate density of tree locations per unit area using a
#' kernel smoothed intensity function
#'
#' @param tree_loc A matrix where each row is an individual observation and the
#'   first column is the x-coordinate and the second column is the y-coordinate.
#'
#' @details An object of class \code{\link[spatstat]{ppp}} is created from the
#'   coordinates in \code{tree_loc} and an pixel image of intensity values is
#'   created using \code{\link[spatstat]{density.ppp}}.
#'
#' @return A summary of the values from pixel image of intensity values.
#'
#' @seealso \link[spatstat]{ppp}, \link[spatstat]{density.ppp}
#' @export
est_density <- function(tree_loc) {
  ppp_tree <- spatstat::ppp(x = unique(tree_loc[, 1]), y = unique(tree_loc[, 2]),
                            window = spatstat::owin(xrange = range(tree_loc[, 1]),
                                                    yrange = range(tree_loc[, 2])));
  dens <- spatstat::density.ppp(ppp_tree);
  return(summary(dens));
}

#' Performs the \emph{walkthrough} edge-correction method for a set of sample
#' trees and a given boundary
#'
#' @param dt_s_tree A data.table object as coming from
#'   \code{\link{extract_data}} function
#' @param dt_s_loc A data.table object as coming from \code{\link{xy_sample}}
#'   function
#' @param boundary A \code{\link[sp]{SpatialPolygons}} object containing one
#'   polygon, describing the boundary of the population.
#' @param on Indicate on which columns \code{dt_s_tree} and \code{dt_s_loc}
#'   should be joined, see \code{\link[data.table]{data.table}}.
#'
#' @return A vector with row indices identifying trees in dt_s_tree that should
#'   be counted twice.
#' @references Ducey, M.J., Gove, J.H., Valentine, H.T., 2004. A walkthrough
#'   solution to the boundary overlap problem. Forest Science 50, 427–435.
#' @import data.table
#' @export
edge_corr_wt <- function(dt_s_tree,
                         dt_s_loc,
                         boundary,
                         on = c("id_set", "id_point")) {
  dt_s_tree <- dt_s_tree[dt_s_loc, on = on];

  # Distance to sample location
  dt_s_tree[,
            ':='(d_x = x_rel - x_pt,
                   d_y = y_rel - y_pt),
            by = on];
  # Walkthrough points
  dt_s_tree[,
            ':='(x_wt = x_rel + d_x,
                   y_wt = y_rel + d_y),
            by = on];

  wt_points <- sp::SpatialPoints(as.matrix(dt_s_tree[!is.na(x_rel), list(x_wt, y_wt)]));
  return(which(rgeos::gWithin(wt_points, boundary, byid = TRUE) == FALSE));
}
