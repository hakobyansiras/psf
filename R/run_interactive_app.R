#' Run KEGG interactive editor and visualization App
#' 
#' @param port the port number to use (Default: 3838)
#' 
#' @return Nothing is returned
#' 
#' @examples
#' \dontrun{
#'     runShinyApp()
#' }
#' 
#' @concept psf
#' @export
run_shiny_app <- function (port=3838) {
  shiny::runApp(system.file('shinyApp', package='psf'), port=port)	
}
