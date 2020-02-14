#' fisim
#'
#' Forest Inventory Simulator (fisim) is a package for simulating response and
#' samplign designs common in forest inventories.
#'
#' @name fisim
#' @docType package
#' @import data.table spatstat sp rgeos
NULL

# avoid R CMD check warnings due to global variables
utils::globalVariables(c("N", "d", "dbh", "ef", "ef_alt1", "ef_alt2", "f_edge",
                         "id_point", "id_sample", "id_stem", "id_tree",
                         "is_edge", "n", "n_corr", "s", "variable", "x_s",
                         "x_tree", "x_wt", "y_s", "y_tree", "y_wt"))
