#' Select trees at sample locations using fixed-area plots
#'
#' @param tree_loc A matrix where each row is an individual observation and the
#'   first column is the x-coordinate and the second column is the y-coordinate.
#' @param sample_loc A matrix with x-coordinates in the first column and
#'   y-coordinates in the second. Each row represents one sample location.
#' @param r The radius of the fixed-area sample plot in the same units as the
#'   coorinates in \code{tree_loc} and \code{sample_loc}.
#' @param k Number of neighbors to search for within \code{r}. See details.
#' @param n Sample size. The number of point locations selected at which trees
#'   should be sampled.
#' @param M Number of independent samples of size \code{n}.
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
#' @return An integer matrix where rows represent sample locations and columns
#'   representing indices of trees in \code{tree_loc} that are selected at the
#'   individual sample locations. Zeroes are used to indicate no neighbours and
#'   to ensure a rectangular data format. See the output of the
#'   \code{\link[RANN]{nn2}} function.
#' @export
fixed_area <- function(tree_loc, sample_loc, r, k = NULL, n, M = 1) {
  if (is.null(k)) {
    dens <- est_density(tree_loc);
    k <- ceiling(dens$max*r^2*pi*2); # Multiply by 2 to make sure all neighbors within r are included
    k <- min(k, nrow(tree_loc));
  }

  nn <- RANN::nn2(data = tree_loc,
                  query = sample_loc,
                  k = k,
                  searchtype = "radius",
                  radius = r);

  s_init <- as.vector(t(nn$nn.idx));
  id_set <- rep(1:M, each = n*ncol(nn$nn.idx)); # sample location id
  id_point <- rep(rep(1:n, each = ncol(nn$nn.idx)), M);
  i_sel <- s_init != 0; # index of selected elements
  s <- s_init[i_sel];
  id_set <- id_set[i_sel];
  id_point <- id_point[i_sel];
  s <- cbind(id_set, id_point, s);
  attr(s, "sample_size") <- n;
  attr(s, "response_design") <- "fixed_area";
  attr(s, "plot_radius") <- r;
  attr(s, "ef") <- 10000/(pi*r^2);
  return(s);
}

#' The \code{k} trees closest to a sample location are selected into the sample
#'
#' @param tree_loc A matrix where each is an individual observation and the
#'   first column is the x-coordinate and the second column is the y-coordinate.
#' @param sample_loc A matrix with x-coordinates in the first column and
#'   y-coordinates in the second. Each row represents one sample location.
#' @param k Number of nearets neighbors.
#' @details The RANN package is used to find the \code{k} trees clostest to the
#'   respective sample locations in \code{sample_loc}. The RANN package is a
#'   wrapper for the Approximate Near Neighbor (ANN) C++ library allowing for
#'   fast nearest-neighbor searches using kd-trees.
#' @return An integer matrix where rows represent sample locations and columns
#'   represent indices of the trees in \code{tree_loc} that are selected at the
#'   individual sample locations. See the output of the \code{\link[RANN]{nn2}}
#'   function.
#' @export
k_tree <- function(tree_loc, sample_loc, k = 6) {

  nn <- RANN::nn2(data = tree_loc,
                  query = sample_loc,
                  k = k,
                  searchtype = "standard");

  s <- as.vector(t(nn$nn.idx));
  names(s) <- rep(1:nrow(nn$nn.idx), each = k);
  attr(s, "sample_size") <- nrow(sample_loc);
  attr(s, "response_design") <- "k_tree";
  attr(s, "k") <- k;
  return(s);
}

#' Select trees proportional to size using Bitterlich/angle-count method
#'
#' @param tree_loc A matrix where each row is an individual observation and the
#'   first column is the x-coordinate and the second column is the y-coordinate.
#' @param tree_dbh A vector with trunk diameters in cm measured at a height of
#'   1.3 m for each observation in \code{tree_loc} (dbh - diameter at breast
#'   height).
#' @param sample_loc A data.table object as returned from
#'   \code{\link{xy_sample}}. The first column is an identifier, the second and
#'   third store sample location coordinates, respectively. The number of rows
#'   corresponds to the size of the sample.
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
angle_count <- function(tree_loc, tree_dbh, sample_loc, baf = 1, k = NULL) {
  c <- sqrt(2500/baf);
  d_crit <- c*tree_dbh/100;

  if (is.null(k)) {
    dens <- est_density(tree_loc);
    k <- ceiling(dens$max*max(d_crit)^2*pi);
    k <- min(k, nrow(tree_loc));
  }

  nn <- RANN::nn2(data = tree_loc,
                  query = sample_loc,
                  k = k,
                  searchtype = "radius",
                  radius = max(d_crit));

  s_init <- as.vector(t(nn$nn.idx));
  sl_id <- rep(1:nrow(nn$nn.idx), each = ncol(nn$nn.idx)); # sample location id
  i_init <- s_init != 0; # index of non-zero elements
  s_init <- s_init[i_init];
  sl_id <- sl_id[i_init];

  d_crit <- d_crit[s_init];
  d_init <- as.vector(t(nn$nn.dists));
  d_init <- d_init[i_init];
  i <- d_init <= d_crit;
  s <- s_init[i];
  sl_id <- sl_id[i];
  names(s) <- sl_id;

  attr(s, "sample_size") <- nrow(sample_loc);
  attr(s, "response_design") <- "angle_count";
  attr(s, "baf") <- baf;
  attr(s, "ef") <- 10000/(d_crit[i]^2*pi);
  return(s);
}
