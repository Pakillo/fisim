#' Launch Shiny App
#'
#' The Shiny App is graphical user interface (GUI) to the functionality of
#' this package. It plots a map of all trees of the considered stand,
#' the sample plots and the trees being part of the sample. The user can select
#' a sampling design and the number of simulations to run. Resulting estimates
#' are displayed in a table and as a histogram.
#'
#' @return Launches the Shiny App and returns a shiny application object.
#'
#' @export
#'
#' @examples
#' \dontrun{launch_shiny_app()}

launch_shiny_app <- function() {
  # wrapper for shiny::shinyApp()
  shiny::shinyApp(ui=shiny_ui, server=shiny_server)
}
