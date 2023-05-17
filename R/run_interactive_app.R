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
run_shiny_app <- function(port=3838) {
  shiny::runApp(system.file('shinyApp', package='psf'), port=port, launch.browser = TRUE)	
}

### need to work on this
#' Run KEGG interactive visualization App for provided pathway
#' 
#' @param pathway list object
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
pathway_shiny_vis <- function(pathway) {
  shiny::runApp(system.file('vis_app', package='psf'), port=port, launch.browser = TRUE)	
}