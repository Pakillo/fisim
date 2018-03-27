#' Objects of class sample_loc
#'
#' The function \code{sample_loc} creates sample location objects, named lists
#' holding coordinates of sample locations spread over some area and
#' accompanying information about how the samples were generated.
#'
#' @param data A \code{\link[data.table]{data.table}} object holding the
#'   cooridnates of sample locations and two identifiers. The first identifier
#'   'id_sample' is used to separate repeated samples of size \code{n} from each
#'   other. There should be \code{M} distinct integers. The second identifier
#'   'id_point' is used to identify individual sample locations in a specific
#'   sample. Both identifiers should be provided as integers. The next two
#'   columns - 'x_s' and 'y_s' - are used for the coordinates and should be
#'   provided as doubles.
#' @param n Sample size
#' @param M Number of repeated samples
#' @param method Were the points spread out randomly or regular
#'
#' @details Before object initialization, a couple of checks are performed,
#'   e.g., column names, data types, data consitency, etc.
#'
#'   In case random sampling was performed, each sample in \code{data} should
#'   have exactly \code{n} sample locations. If systematic sampling was applied,
#'   sample size varies from sample to sample depending on the shape of the
#'   study region and sample size. If the expected sample size (\code{n})
#'   provided by the user, deviates by more than 1 from the average sample size
#'   inferred from \code{data}, a warning is thrown rather than an error.
#'
#'   The prefered way for creating objects of type sample_loc is via the
#'   \code{\link{xy_sample}} function.
#'
#' @return An object of class \code{sample_loc}. A named list with the following
#'   four entries:
#'
#'   \enumerate{
#'   \item \code{$data} A \code{\link[data.table]{data.table}}
#'   object containing two identifier columns and two columns for storing
#'   coordinates
#'   \item \code{$n} Sample size
#'   \item \code{$M} Number of repeated
#'   samples
#'   \item \code{$method} Random or regular placement of sample
#'   locations
#'   }
#' @export
#'
#' @examples
#' library(spatstat);
#' xy_sample(letterR, n = 10, M = 10, method = "random");
sample_loc <- function(data, n, M, method) {
  if (!is.data.table(data)) {
    stop("'data' is not a data.table object!");
  }
  if (!identical(names(data), c("id_sample", "id_point", "x_s", "y_s"))) {
    stop("Column names of 'data' do not meet class requirements!");
  }
  if (M != data[, length(unique(id_sample))]) {
    stop("Number of repated samples in 'data' <length(unique(id_sample))> and 'M' differ!");
  }
  if (is.na(match(method, c("random", "regular")))) {
    stop("'method' is not one of 'random' or 'regular'");
  }
  if (method == "random") {
    if (n != data[, .N, by = id_sample][, mean(N)]) {
      stop("Sample size 'n' differs from 'data'!");
    }
  } else if (method == "regular") {
    n_data <- data[, .N, by = id_sample][, mean(N)];
    if (abs(n - n_data > 1)) {
      warning("Expected sample size 'n' differs from data average by more than 1!");
    }
  }

  sample_loc <- list(data = data,
                     n = n,
                     M = M,
                     method = method);
  class(sample_loc) <- "sample_loc";
  return(sample_loc);
}

#' Check if an object is a sample locations object.
#'
#' @param x An \code{R} object
#'
#' @return TRUE if \code{x} is of class \code{\link{sample_loc}};
#' @export
is.sample_loc <- function(x) {
  inherits(x, "sample_loc");
}

#' Objects of class tree_pop
#'
#' The function \code{tree_pop} creates tree population objects holding data on
#' individual trees forming a forest and a polygon that defines the border of
#' the forest. The objects are used as an input for simulating forest
#' inventories.
#'
#' @param data A \code{\link[data.table]{data.table}} object that describes the
#'   individual trees in a forest. Each row represents a tree, while tree
#'   attributes are stored in columns. A minimum set of attributes with specific
#'   columns names is required; see Details section.
#' @param boundary A \code{\link[sp]{SpatialPolygons}} object that defines the
#'   border of the forest.
#'
#' @details The following columns are needed as a minimum requirement in
#'   \code{data}:
#'   \itemize{
#'   \item \code{id_tree} - integer values used to
#'   identify individual trees
#'   \item \code{id_stem} - integer values to identify
#'   individual stems. A tree can fork into multiple stems and stems are treated
#'   as individual record under certain circumstances (exact definitions depend
#'   on the particular surveys at which trees were measured). If a tree has two
#'   stems, for example, \code{id_tree} occurs twice for that tree, while
#'   \code{id_stem} is unique for each stem.
#'   \item \code{x_tree} The x-coordinates of the trees, defined at the stem
#'   center at a height of 1.3 meter above ground level. The data type is
#'   numeric (double) and the units are usually in meter.
#'   \item \code{y_tree} The y-coordinates of the trees, defined at the stem
#'   center at a height of 1.3 meter above ground level. The data type is
#'   numeric (double) and the units are usually in meter.
#'   \item \code{is_dead} A logical indicating whether the stem is dead
#'   (\code{TRUE}) or alive (\code{FALSE}).
#'   \item \code{dbh} The diameter of the stem, measured at a height of
#'   1.3 meter above ground in cm (diameter at breast height, \code{numeric}).
#'   }
#'
#'   During object initialization, three variables required for simulations will
#'   be added:
#'   \itemize{
#'   \item \code{ba} The cross-section of the stems at a height of 1.3 meter
#'   above ground level calculated from the provided \code{dbh} values assuming
#'   a circular stem.
#'   \item \code{n} The number of stems; initialized with one for each entry.
#'   This is simply a helper variable for simualtions and estimation.
#'   \item \code{f_edge} Edge factor, used for edge correction of trees close
#'   to the forest border. Since such trees have a lower probability of being
#'   included into the sample due to the restriction of the sampling to the
#'   study area defined by \code{boundary}, edge correction must be performed
#'   to achieve unbiase estimation. The variable is initialized with one and
#'   will be modified during the simulation according to appropriate edge
#'   correction methods chosen by the user. After edge correction, some trees
#'   are counted twice to compensate for the lower inlcusion probabilities of
#'   border trees.
#'   }
#'
#'   This minimum set of attributes is required for the estimation of stem
#'   counts, basal area (a measure for site occupation) and diameter
#'   distributions. Users may add additional variables like species,
#'   tree volume, aboveground biomass, etc. It will be possible to include such
#'   addtional variables in the simulation.
#'
#' @return An object of class \code{tree_pop}. A named list with the following
#'   two entries:
#'
#'   \enumerate{
#'   \item \code{$data} A \code{\link[data.table]{data.table}}
#'   object as described above.
#'   \item \code{$boundary} A \code{\link[sp]{SpatialPolygons}} object that
#'   defines the border of the forest.
#'   }
#' @export
tree_pop <- function(data, boundary) {
  if (!is.data.table(data)) {
    stop("'data' is not a data.table object!");
  }
  sp_class <- c("SpatialPolygons", "SpatialPolygonsDataFrame");
  if (is.na(match(class(boundary), sp_class))) {
    stop("'boundary' neither of class 'SpatialPolygons' nor of class 'SpatialPolygonsDataFrame'");
  }
  col_names <- c("id_tree", "id_stem", "x_tree", "y_tree", "is_dead", "dbh");
  if (sum(is.na(match(col_names, names(data)))) > 0) {
    stop("Column names of 'data' do not meet class requirements!");
  }
  x_range <- data[, range(x_tree)];
  y_range <- data[, range(y_tree)];
  ext <- bbox(boundary);
  if (x_range[1] < ext["x", "min"]) {
    stop("x-coordinate in 'data' smaller than extent of 'boundary'");
  }
  if (x_range[2] > ext["x", "max"]) {
    stop("x-coordinate in 'data' larger than extent of 'boundary'");
  }
  if (y_range[1] < ext["x", "min"]) {
    stop("y-coordinate in 'data' smaller than extent of 'boundary'");
  }
  if (y_range[2] > ext["x", "max"]) {
    stop("y-coordinate in 'data' larger than extent of 'boundary'");
  }

  tree_pop <- list(data = data,
                   boundary = boundary);
  tree_pop$data[, ':='(ba = (dbh/100)^2*pi/4,
                       n = 1L,
                       f_edge = 1L)];
  class(tree_pop) <- "tree_pop";
  return(tree_pop);
}


#' Check if an object is a tree population.
#'
#' @param x An \code{R} object
#'
#' @return TRUE if \code{x} is of class \code{\link{tree_pop}};
#' @export
is.tree_pop <- function(x) {
  inherits(x, "tree_pop");
}

#' Objects of class response
#'
#' The function \code{response} creates response design objects that store the
#' results of a particular response design that was applied at given sample
#' locations to some tree population. The function is only used internally by
#' the functions that apply the particular response designs:
#' \code{\link{fixed_area}}, \code{\link{k_tree}}, \code{\link{angle_count}}.
#'
#' @param data A \code{\link[data.table]{data.table}} object holding identifiers
#'   for the sample and respective point locations and row indices that indicate
#'   which population elements are selected into the sample.
#' @param r_design A \code{character} string indicating the type of the response
#'   design. One of the following: "fixed_area", "k_tree", "angle_count".
#' @param r_design_parm A response design specific parameter of type
#'   \code{numeric}. For "fixed_area" the plot radius in meter; for "k_tree" the
#'   number of trees closest to the sample location; for "angle_count" the basal
#'   area factor.
#' @param ef The expansion factor, i.e., the factor that is used to prorate tree
#'   attributes to per hectare values (the inverse of the inclusion
#'   probabilities). Provided by the particular response design functions.
#' @param ef_alt1 Alternative expansion factor that may be applied when using
#'   k-tree sampling.
#' @param ef_alt2 A third alternative expansion factor approximation for k-tree
#'   sampling
#'
#' @details The \code{\link[data.table]{data.table}} object (\code{data}) should
#'   have the following columns:
#'   \itemize{
#'   \item \code{id_sample}
#'   \item \code{id_point}
#'   \item \code{s}
#'   \item \code{ef}
#'   }
#'
#'   Some response designs (currently only \code{\link{k_tree}}), may add
#'   additional alternative expansion factors as optional approximations, where
#'   unbiased estimators are missing or difficult to derive.
#'
#' @return An object of class \code{response}. A named list with the following
#' entries:
#'
#'   \itemize{
#'   \item \code{$data} A \code{\link[data.table]{data.table}}
#'   object as described above.
#'   \item \code{$r_design}
#'   \item \code{$r_design_parm}
#'   }
response <- function(data, r_design, r_design_parm) {
  if (!is.data.table(data)) {
    stop("'data' is not a data.table object!");
  }
  col_names <- c("id_sample", "id_point", "s", "ef");
  if (sum(is.na(match(col_names, names(data)))) > 0) {
    stop("Column names of 'data' do not meet class requirements!");
  }
  response <- list(data = data,
                   r_design = r_design,
                   r_design_parm = r_design_parm);
  class(response) <- "response";
  return(response);
}

#' Check if an object is a response design object.
#'
#' @param x An \code{R} object
#'
#' @return TRUE if \code{x} is of class \code{\link{response}};
#' @export
is.response <- function(x) {
  inherits(x, "response");
}

#' The function \code{tree_sample} creates tree sample objects. Such objects
#' assign individual tree information from the population to the sample
#' locations following a particular response design. The function is only used
#' internally by \code{\link{extract_data}}.
#'
#'
#' @param data A \code{\link[data.table]{data.table}} object holding individual
#'   tree data together with identifiers that indicate the sample and the point
#'   location at which trees were selected.
#' @param r_design A \code{character} string indicating the type of the response
#'   design. One of the following: "fixed_area", "k_tree", "angle_count".
#' @param r_design_parm A response design specific parameter of type
#'   \code{numeric}. For "fixed_area" the plot radius in meter; for "k_tree" the
#'   number of trees closest to the sample location; for "angle_count" the basal
#'   area factor.
#'
#' @details The \code{\link[data.table]{data.table}} object (\code{data}) should
#'   have the following columns:
#'   \itemize{
#'   \item \code{id_sample}
#'   \item \code{id_point}
#'   \item \code{id_tree}
#'   \item \code{id_stem}
#'   \item \code{x_tree}
#'   \item \code{y_tree}
#'   \item \code{is_dead}
#'   \item \code{dbh}
#'   \item \code{ba}
#'   \item \code{n}
#'   \item \code{f_edge}
#'   \item \code{ef}
#'   }
#'
#'   Some response designs (currently only \code{\link{k_tree}}), may have
#'   additional alternative expansion factors as optional approximations, where
#'   unbiased estimators are missing or difficult to derive. As the column names
#'   for these alternatives \code{ef_alt1} and \code{ef_alt2} are used.
#'
#' @return
#'  \enumerate{
#'   \item \code{$data} A \code{\link[data.table]{data.table}}
#'   object as described above.
#'   \item \code{$r_design}
#'   \item \code{$r_design_parm}
#'   }
tree_sample <- function(data, r_design, r_design_parm) {
  if (!is.data.table(data)) {
    stop("'data' is not a data.table object!");
  }
  col_names <- c("id_sample", "id_point", "id_tree", "id_stem", "x_tree",
                 "y_tree", "is_dead", "dbh", "ba", "n", "f_edge", "ef");
  if (sum(is.na(match(col_names, names(data)))) > 0) {
    stop("Column names of 'data' do not meet class requirements!");
  }

  tree_sample <- list(data = data,
                      r_design = r_design,
                      r_design_parm = r_design_parm);
  class(tree_sample) <- "tree_sample";
  return(tree_sample);
}


#' Check if an object is a tree sample object.
#'
#' @param x An \code{R} object
#'
#' @return TRUE if \code{x} is of class \code{\link{tree_sample}};
#' @export
is.tree_sample <- function(x) {
  inherits(x, "tree_sample");
}


#' Create objects of class \code{point_data}
#'
#' The function \code{point_data} creates point data objects. Such objects store
#' aggregated tree-level information by sample and sample point location. The
#' function is only used internally by \code{\link{sum_data}}.
#'
#'
#' @param data A \code{\link[data.table]{data.table}} object holding aggregated
#'   tree data together with identifiers that indicate the sample and the point
#'   location.
#' @param r_design A \code{character} string indicating the type of the response
#'   design. One of the following: "fixed_area", "k_tree", "angle_count".
#' @param r_design_parm A response design specific parameter of type
#'   \code{numeric}. For "fixed_area" the plot radius in meter; for "k_tree" the
#'   number of trees closest to the sample location; for "angle_count" the basal
#'   area factor.
#'
#' @details The \code{\link[data.table]{data.table}} object (\code{data}) should
#'   have the following columns:
#'   \itemize{
#'   \item \code{id_sample}
#'   \item \code{id_point}
#'   \item \code{variable}
#'   }
#'
#'   When \code{\link{k_tree}} sampling was applied at the individual sample
#'   locations, there are four approximations for expanding the target variables
#'   listed under the \code{variable} column (see \code{\link{sum_k_tree}}).
#'   The respective column names are:
#'   \itemize{
#'   \item \code{Prodan}
#'   \item \code{Eberhardt}
#'   \item \code{d_mean}
#'   \item \code{a_mean}
#'   }
#'
#'   For other response designs (e.g., angle count, fixed area), expansion to
#'   per hectare values is based on known expansion factors that can be derived
#'   from the response design parameters and not on approximations. In such
#'   cases, the expanded values are listed under \code{ef_ana}, meaning
#'   analytically derived expansion factors.
#'
#'   The last column \code{is_edge} indicates whether some sort of edge
#'   correction was performed at the respective sample point or not.
#'
#' @return
#' \enumerate{
#' \item \code{$data} A \code{\link[data.table]{data.table}} object as
#' described above.
#' \item \code{$r_design}
#' \item \code{$r_design_parm}
#' }
point_data <- function(data, r_design, r_design_parm) {
  if (!is.data.table(data)) {
    stop("'data' is not a data.table object!");
  }
  col_names <- c("id_sample", "id_point", "variable");
  if (sum(is.na(match(col_names, names(data)))) > 0) {
    stop("Column names of 'data' do not meet class requirements!");
  }

  point_data <- list(data = data,
                     r_design = r_design,
                     r_design_parm = r_design_parm);
  class(point_data) <- "point_data";
  return(point_data);
}


#' Check if an object is a point data object.
#'
#' @param x An \code{R} object
#'
#' @return TRUE if \code{x} is of class \code{\link{point_data}};
#' @export
is.point_data <- function(x) {
  inherits(x, "point_data");
}