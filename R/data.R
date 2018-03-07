#' Full census of a peat swamp forest in Kalimantan on Borneo island
#'
#' A dataset containing tree stem locations and other attributes of 3523 trees
#' measured on 120m by 120m large full census plot in a peat swamp forest on
#' Borneo island.
#'
#' @format A data.table with 3523 rows and 8 columns: \describe{
#'   \item{id_tree}{Tree number - repetitions possible due to multiple stems}
#'   \item{id_stem}{Stem number - one tree can have mulitple stems}
#'   \item{x_tree}{Relative position of the stem center in x-direction}
#'   \item{y_tree}{Relative position of the stem center in y-direction}
#'   \item{is_dead}{Was the tree alive or dead?}
#'   \item{dbh}{Diameter of the stem measured at breast height (1.3m) in cm}
#'   \item{tree_height}{Height of the trees in m}
#'   \item{crownwidth}{Width of the tree cronwns in m}
#'   \item{ba}{Cross-sectional area of the stem at a height of 1.3m in square
#'     meter}
#'   \item{n}{Stem count - helper variable for stem number estimation, unitless}
#'   \item{f_edge}{Edge factor - helper variable for edge_correction, indicates
#'     how often a tree is counted during estimation correcting for edge
#'     effects, unitless}}
"kalimantan_peat"


#' Full census of a beech forest in Göttingen, Germany
#'
#' A dataset containing tree stem locations and other attributes of 270 trees
#' measured on a 120m by 120m large full census plot in a beech forest near
#' Göttingen, Germany.
#'
#' @format An object of class \code{\link{tree_pop}}. The data contains 270
#'   rows and 9 columns: \describe{
#'   \item{id_tree}{Stem number - one tree can have mulitple stems}
#'   \item{id_stem}{Tree number - repetitions possible due to multiple stems}
#'   \item{x_tree}{Relative position of the stem center in x-direction}
#'   \item{y_tree}{Relative position of the stem center in y-direction}
#'   \item{is_dead}{Was the tree alive or dead?}
#'   \item{dbh}{Diameter of the stem measured at breast height (1.3m) in cm}
#'   \item{ba}{Cross-sectional area of the stem at breast height im square
#'     meters}
#'   \item{n}{Stem count - helper variable for stem number estimation, unitless}
#'   \item{f_edge}{Edge factor - helper variable for edge_correction, indicates
#'     how often a tree is counted during estimation correcting for edge
#'     effects, unitless}}
"hberg_beech"
