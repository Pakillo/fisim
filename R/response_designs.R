#' Select trees at sample locations using fixed-area plots
#'
#' @param tree_pop A tree population object (\code{\link{tree_pop}})
#' @param sample_loc A sample location object (\code{\link{sample_loc}}) as
#'   provided by the \code{\link{xy_sample}} function.
#' @param r The radius of the fixed-area sample plot in the same units as the
#'   coorinates in \code{tree_loc} and \code{sample_loc}.
#' @param k Number of neighbors to search for within \code{r}. See details.
#'
#' @details The RANN package is used to find all trees within the specified
#'   radius. The RANN package is a wrapper for the Approximate Near Neighbor
#'   (ANN) C++ library allowing for fast nearest-neighbor searches using
#'   kd-trees. Here, radius search is used, meaning that the \code{k}-nearest
#'   neighbors within \code{r} are searched for. To find all trees within
#'   \code{r}, \code{k} should be large enough (max the number of rows in
#'   \code{tree_loc}). If not specified by the user, \code{k} is estimated by
#'   multiplying the average tree density with plot area. Tree density is
#'   inferred using the \code{density.ppp} function from the spatstat package.
#'   Note that the estimation costs some extra computation time and it is thus
#'   recommended to provide a large enough \code{k} when the function is used in
#'   simulations.
#' @return An object of class \code{\link{response}} with indices of trees in
#'   \code{tree_pop} that are selected at the individual sample locations.
#' @export
fixed_area <- function(tree_pop, sample_loc, r, k = NULL) {
  if (!is.tree_pop(tree_pop)) {
    stop("'tree_pop' is not an object of class tree_pop");
  }
  if (!is.sample_loc(sample_loc)) {
    stop("'sample_loc' is not an object of class sample_loc");
  }
  if (is.null(k)) {
    dens <- est_density(as.matrix(tree_pop$data[, list(x_tree, y_tree)]));
    k <- ceiling(dens$max*r^2*pi*2); # Increase by 100% to make sure all neighbors within r are included
    k <- min(k, nrow(tree_pop$data));
  }

  nn <- RANN::nn2(data = tree_pop$data[, list(x_tree, y_tree)],
                  query = sample_loc$data[, list(x_s, y_s)],
                  k = k,
                  searchtype = "radius",
                  radius = r);

  dt_s <- data.table(id_sample = rep(sample_loc$data[, id_sample], each = k),
                     id_point = rep(sample_loc$data[, id_point], each = k),
                     s = as.vector(t(nn$nn.idx)));
  dt_s <- dt_s[s != 0];

  ha <- 10000; # One hectare -> 10000 square meter
  dt_s[, ef := ha/(pi*r^2)];

  # Add no response
  dt_s <- dt_s[sample_loc$data[, list(id_sample, id_point)],
               list(id_sample, id_point, s, ef),
               on = c("id_sample", "id_point")];

  return(response(data = dt_s,
                  r_design = "fixed_area",
                  r_design_parm = r));
}


#' The \code{k} trees closest to a sample location are selected into the sample
#'
#' @param tree_pop A tree population object (\code{\link{tree_pop}})
#' @param sample_loc A sample location object (\code{\link{sample_loc}}) as
#' provided by the \code{\link{xy_sample}} function.
#' @param k Number of nearets neighbors.
#' @details The RANN package is used to find the \code{k} trees clostest to the
#'   respective sample locations in \code{sample_loc}. The RANN package is a
#'   wrapper for the Approximate Near Neighbor (ANN) C++ library allowing for
#'   fast nearest-neighbor searches using kd-trees.
#' @return An object of class \code{\link{response}} with indices of trees in
#'   \code{tree_pop} that are selected at the individual sample locations.
#' @export
k_tree <- function(tree_pop, sample_loc, k = 6) {
  if (!is.tree_pop(tree_pop)) {
    stop("'tree_pop' is not an object of class tree_pop");
  }
  if (!is.sample_loc(sample_loc)) {
    stop("'sample_loc' is not an object of class sample_loc");
  }

  nn <- RANN::nn2(data = tree_pop$data[, list(x_tree, y_tree)],
                  query = sample_loc$data[, list(x_s, y_s)],
                  k = k + 1,
                  searchtype = "standard");

  dt_s <- data.table(id_sample = rep(sample_loc$data[, id_sample], each = k),
                     id_point = rep(sample_loc$data[, id_point], each = k),
                     s = as.vector(t(nn$nn.idx[, 1:k])));

  d_k <- nn$nn.dists[, k];
  d_k1 <- nn$nn.dists[, k + 1];
  ha <- 10000; # One hectare -> 10000 square meter
  dt_ef <- data.table(id_sample = sample_loc$data[, id_sample],
                      id_point = sample_loc$data[, id_point],
                      ef = ha/(pi*d_k^2),
                      ef_alt1 = ha/(pi*(0.5*(d_k + d_k1))^2),
                      ef_alt2 = ha/(pi*(sqrt(0.5*(d_k^2 + d_k1^2)))^2));

  return(response(data = dt_s[dt_ef, on = c("id_sample", "id_point")],
                  r_design = "k_tree",
                  r_design_parm = k));
}


#' Select trees proportional to size using Bitterlich/angle-count method
#'
#' @param tree_pop A tree population object (\code{\link{tree_pop}})
#' @param sample_loc A sample location object (\code{\link{sample_loc}}) as
#'   provided by the \code{\link{xy_sample}} function.
#' @param baf Basal area factor. Determines together with the trunk diameter at
#'   1.3 m height up to which distance trees are included in the sample. Typical
#'   values for temperate forest conditions range between 1 and 4. The number of
#'   trees included at a sample location multiplied with the \code{baf} gives a
#'   direct estimate of basal area in square meter per hectare.
#' @param k Helper variable needed in the \code{\link{nn2}} function for the
#'   kd-tree search of candidate trees that are possibly included into the
#'   sample. Maximum the numer of rows in \code{tree_loc}. See details.
#' @details The selection of trees according to the Bitterlich/ange-count method
#'   was implemented in several steps. First the critical distance up to which a
#'   tree is selected from a specific sample location is calculated using
#'   \code{tree_dbh} and \code{baf}. The larger \code{tree_dbh} and \code{baf},
#'   the larger the critical distance of inclusion. In the second step the RANN
#'   package is used to identify candidate trees for inlcusion into the sample
#'   and their distances to the sample locations. The RANN package is a wrapper
#'   for the Approximate Near Neighbor (ANN) C++ library allowing for fast
#'   nearest-neighbor searches using kd-trees. Here, radius search is used,
#'   meaning that the \code{k}-nearest neighbors within a specific radius are
#'   searched for. As a radius the maximum critial distance is used. In
#'   combination with \code{k}, a larger than necessary sample of trees is
#'   selected at each sample location. If not specified by the user, \code{k} is
#'   estimated by multiplying the average tree density with the area of a circle
#'   with the maximum critical distance as radius. Tree density is inferred
#'   using the \code{density.ppp} function from the spatstat package. Note that
#'   the estimation costs some extra computation time and it is thus recommended
#'   to provide a large enough \code{k} when the function is used in
#'   simulations. In the third and last step, the set of candidate trees is
#'   querried for trees, where the actual distance to the sample location is
#'   smaller than the individual critical distance.
#' @return An integer matrix where rows represent sample locations and columns
#'   represent indices of the trees in \code{tree_loc} that are selected at the
#'   individual sample locations. Zeroes are used to indicate no neighbours and
#'   to ensure a rectangular data format. See the output of the
#'   \code{\link[RANN]{nn2}} function.
#' @export
angle_count <- function(tree_pop, sample_loc, baf = 1, k = NULL) {
  c <- sqrt(2500/baf);
  d_crit <- c*tree_pop$data[, dbh]/100;

  if (is.null(k)) {
    dens <- est_density(as.matrix(tree_pop$data[, list(x_tree, y_tree)]));
    k <- ceiling(dens$max*max(d_crit)^2*pi);
    k <- min(k, nrow(tree_pop$data));
  }

  nn <- RANN::nn2(data = tree_pop$data[, list(x_tree, y_tree)],
                  query = sample_loc$data[, list(x_s, y_s)],
                  k = k,
                  searchtype = "radius",
                  radius = max(d_crit));

  dt_s <- data.table(id_sample = rep(sample_loc$data[, id_sample], each = k),
                     id_point = rep(sample_loc$data[, id_point], each = k),
                     s = as.vector(t(nn$nn.idx)),
                     d = as.vector(t(nn$nn.dists)));
  dt_s <- dt_s[s != 0];
  dt_s[, d_crit := d_crit[s]];
  dt_s <- dt_s[d <= d_crit];

  ha <- 10000; # One hectare -> 10000 square meter
  dt_s[, ef := ha/(dt_s[, d_crit^2*pi])]

  # Add no response
  dt_s <- dt_s[sample_loc$data[, list(id_sample, id_point)],
               list(id_sample, id_point, s, ef),
               on = c("id_sample", "id_point")];

  return(response(data = dt_s,
                  r_design = "angle_count",
                  r_design_parm = baf));
}
