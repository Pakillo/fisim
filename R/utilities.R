#' Helper function to extract population elements into a sample using a list of
#' indices
#'
#' @param tree_pop The population from which to extract sampled elements.
#' @param response An integer matrix where rows represent sample locations and columns
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
extract_data <- function(tree_pop, response) {
  if (!is.tree_pop(tree_pop)) {
    stop("'tree_pop' is not an object of class tree_pop");
  }
  if (!is.response(response)) {
    stop("'response' is not an object of class response");
  }
  dt_s_tree <- tree_pop$data[response$data[, s]];
  dt_s_tree[, ':='(id_sample = response$data[, id_sample],
                   id_point = response$data[, id_point],
                   ef = response$data[, ef])];
  if (response$r_design == "k_tree") {
    dt_s_tree[, ':='(ef_alt1 = response$data[, ef_alt1],
                     ef_alt2 = response$data[, ef_alt2])];
  }

  return(tree_sample(data = dt_s_tree,
                     r_design = response$r_design,
                     r_design_parm = response$r_design_parm));
}


#' Summarizing of tree-level
#'
#' Helper function to summarize tree-level sample data to the plot-level for
#' estimation of population parameters
#'
#' @param tree_sample A \code{\link{tree_sample}} object as coming from
#'   the \code{\link{extract_data}} function.
#' @param target_vars A character vector of variable names which should be
#'   summarized at the plot level
#'
#' @return A data.table object with stand characteristics at the plot level.
#'   Currently, basal area per hectare and the number of trees per hectare are
#'   included.
#' @export
#' @import data.table
sum_data <- function(tree_sample, target_vars) {
  if (!is.tree_sample(tree_sample)) {
    stop("'tree_sample' is not an object of class tree_sample!");
    }

  # Make sure that plots with no trees receive entries of zero
  for (c in c(target_vars, "f_edge", "ef")) {
    tree_sample$data[is.na(get(c)), (c) := 0];
  }

  if (tree_sample$r_design == "k_tree") {
    dt_s_plot <- sum_k_tree(tree_sample = tree_sample, target_vars = target_vars);
  } else if (tree_sample$r_design != "k_tree") {
    dt_s_plot <- sum_ef(tree_sample = tree_sample, target_vars = target_vars);
  }

  dt_edge <- tree_sample$data[,
                              list(n = sum(n), n_corr = sum(f_edge)),
                              by = list(id_sample, id_point)];
  dt_edge[, is_edge := 0];
  dt_edge[n < n_corr, is_edge := 1];
  return(point_data(data = dt_s_plot[dt_edge[,
                                             list(id_sample, id_point, is_edge)],
                                     on = c("id_sample", "id_point")],
         r_design = tree_sample$r_design,
         r_design_parm = tree_sample$r_design_parm));
}


#' Summarizing of tree-level data for k-tree sampling
#'
#' Helper function to summarize tree-level sample data to the plot-level, where
#' k-tree sampling (\code{\link{k_tree}}) was used for selecting trees at sample
#' locations.
#'
#'
#' @param tree_sample A \code{\link{tree_sample}} object as coming from the
#'   \code{\link{extract_data}} function.
#' @param target_vars A character vector of variable names which should be
#'   summarized at the plot level
#'
#' @details For k-tree sampling, individual tree inlcusion probabilities can
#'   only be derived if the tree stem positions of a substantially larger number
#'   of trees are available than are actually selected at the single sample
#'   locations (see Kleinn & Vilcko 2006a). In practical timber surveys or
#'   forest inventories, such extensive numbers of tree positions around the
#'   actual sample of trees are typically not available. Instead approximations
#'   for expanding the locally observed forest stand characteristics are
#'   applied. This functions provides a set of four approximations that are
#'   presented in Kleinn & Vilcko (2006b).
#'
#' @return A data.table object with stand characteristics at the plot level.
#'   Currently, basal area per hectare and the number of trees per hectare are
#'   included. The four different approximations are provided under the
#'   following columns:
#'   \itemize{
#'   \item \code{Prodan} - The approach of Prodan, where only half of the k-th
#'   tree is considered (Prodan, 1968).
#'   \item \code{Eberhardt} - Expansion is based on the distance to the k-th
#'   tree and corrected with factor (k-1)/k (Eberhardt 1967)
#'   \item \code{d_mean} - A circle, where the radius is the average between
#'   the distances to tree k and tree k+1
#'   \item \code{a_mean} - The radius of the average area of the circles that
#'   are built from the distances to tree k and tree k+1
#'   }
#' @references
#' \itemize{
#' \item Kleinn, C., Vilčko, F., 2006a. Design-unbiased estimation for
#' point-to-tree distance sampling. Canadian Journal of Forest Research 36,
#' 1407–1414. https://doi.org/10.1139/x06-038
#' \item Kleinn, C., Vilčko, F., 2006b. A
#' new empirical approach for estimation in k-tree sampling. Forest Ecology and
#' Management 237, 522–533. https://doi.org/10.1016/j.foreco.2006.09.072
#' \item Prodan, M., 1968. Punktstichprobe für die Forsteinrichtung (A point sample for
#' forest management planning). Forst Holzwirt 23 (11), 225–226.
#' \item Eberhardt, L.L., 1967. Some developments in distance sampling. Biometrics 23,
#' 207–216.
#' }
#' @import data.table
sum_k_tree <- function(tree_sample, target_vars) {
  if (!is.tree_sample(tree_sample)) {
    stop("'tree_sample' is not an object of class tree_sample!");
  }

  n_prodan <- rep(1, tree_sample$r_design_parm);
  n_prodan[tree_sample$r_design_parm] <- 0.5;
  cf_eberhardt <- (tree_sample$r_design_parm - 1)/tree_sample$r_design_parm;
  dt_s_plot <- tree_sample$data[,
                                list(variable = target_vars,
                                     Prodan = sapply(.SD, function(x) sum(ef*n_prodan*x)),
                                     Eberhardt = sapply(.SD, function(x) sum(ef*x)*cf_eberhardt),
                                     d_mean = sapply(.SD, function(x) sum(ef_alt1*x)),
                                     a_mean = sapply(.SD, function(x) sum(ef_alt2*x))),
                                keyby = c("id_sample", "id_point"),
                                .SDcols = target_vars];
  return(dt_s_plot);
}


#' Expansion-factor-based summarizing of tree-level data
#'
#' Helper function to summarize tree-level sample data to the plot-level for
#' estimation of population parameters, where the exact expansion factors are
#' assumed to be known.
#'
#' @param tree_sample A \code{\link{tree_sample}} object as coming from
#'   the \code{\link{extract_data}} function.
#' @param target_vars A character vector of variable names which should be
#'   summarized at the plot level
#'
#' @return A data.table object with stand characteristics at the plot level.
#'   Currently, basal area per hectare and the number of trees per hectare are
#'   included.
#' @import data.table
sum_ef <- function(tree_sample, target_vars) {
  if (!is.tree_sample(tree_sample)) {
    stop("'tree_sample' is not an object of class tree_sample!");
  }

  dt_s_plot <- tree_sample$data[,
                                list(variable = target_vars,
                                     ef_ana = sapply(.SD, function(x) sum(ef*f_edge*x))),
                                keyby = c("id_sample", "id_point"),
                                .SDcols = target_vars];
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
  ppp_tree <- spatstat::ppp(x = tree_loc[, 1], y = tree_loc[, 2],
                            window = spatstat::owin(xrange = range(tree_loc[, 1]),
                                                    yrange = range(tree_loc[, 2])));
  dens <- spatstat::density.ppp(ppp_tree);
  return(summary(dens));
}

#' Performs the \emph{walkthrough} edge-correction method for a set of sample
#' trees and a given boundary
#'
#' @param tree_pop A \code{\link{tree_pop}} object.
#' @param tree_sample A \code{\link{tree_sample}} object as coming from the
#'   \code{\link{extract_data}} function.
#' @param sample_loc A \code{\link{sample_loc}} object as created by the
#'   \code{\link{xy_sample}} function or provided by the user.
#'
#' @return The \code{tree_sample} object with modified \code{f_edge} column
#'   indicating which trees are counted twice.
#' @references Ducey, M.J., Gove, J.H., Valentine, H.T., 2004. A walkthrough
#'   solution to the boundary overlap problem. Forest Science 50, 427–435.
#' @import data.table
#' @export
edge_corr_wt <- function(tree_pop, tree_sample, sample_loc) {
  if (!is.tree_pop(tree_pop)) {
    stop("'tree_pop' is not an object of class tree_pop");
  }
  if (!is.sample_loc(sample_loc)) {
    stop("'sample_loc' is not an object of class sample_loc");
  }
  if (!is.tree_sample(tree_sample)) {
    stop("'tree_sample' is not an object of class tree_sample");
  }
  if (tree_sample$r_design == "k_tree") {
    warning("Edge correction for k-tree sampling not meaningful! Consider to
            run simulations for k-tree sampling without edge correction.");
  }

  dt_s_tree <- tree_sample$data[sample_loc$data,
                                list(id_sample, id_point, id_tree, id_stem, x_s, y_s, x_tree, y_tree),
                                on = c("id_sample", "id_point")];

  # Walkthrough points
  dt_s_tree[,
            ':='(x_wt = 2*x_tree - x_s,
                   y_wt = 2*y_tree - y_s),
            by = list(id_sample, id_point)];

  # Identify points outside border
  wt_points <- sp::SpatialPoints(as.matrix(dt_s_tree[!is.na(x_tree), list(x_wt, y_wt)]));
  idx <- which(rgeos::gWithin(wt_points, tree_pop$boundary, byid = TRUE) == FALSE);

  return(tree_sample(tree_sample$data[dt_s_tree[!is.na(x_wt)][idx, list(id_sample, id_point, id_tree, id_stem)],
                                      f_edge := 2,
                                      on = c("id_sample", "id_point", "id_tree", "id_stem")],
         r_design = tree_sample$r_design,
         r_design_parm = tree_sample$r_design_parm));
}


#' Helper function to extract the area slots from
#' \code{\link[sp]{SpatialPolygons}} objects
#'
#' @param sp_poly An object of class \code{\link[sp]{SpatialPolygons}}
#'
#' @return A numeric vector with the areas of all polygons in the object. The
#'   areas for polygons that represent holes are stored as negative values.
#' @examples
#' library(sp);
#' library(maptools);
#' data("state.vbm");
#' a <- extract_area(state.vbm);
#' @export
extract_area <- function(sp_poly) {
  return(unlist(lapply(sp_poly@polygons,
                       function(x) lapply(x@Polygons,
                                          function(x) x@area*x@ringDir))));
}
