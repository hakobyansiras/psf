#' Calculates psf for given kegg pathway based on expression matrix and generates pdf report with colored pathways and plots
#' @param kegg_collection list of kegg pathways
#' @param exp_matrix expression fold change matrix with gene entrez id rownames
#' @param file_name name of the pdf report which will be generated
#' @param use_old_images use olde kegg images(for use with curated pathway collection)
#' @import gplots
#' @import ggplot2
#' @import magick
#' @import grDevices
#' @export
calc_psf_and_generate_report_from_collection <- function(kegg_collection, exp_matrix, file_name, use_old_images = F) {
  
  plot_list <- lapply(names(kegg_collection), function(x) {
    
    psf_output <- psf_signal_calculator_and_coloring(entrez_fc = exp_matrix, 
                                                     pathway = kegg_collection[[x]], 
                                                     pathway_name = x, no_color_mode = F)
    
    graphical_data <- graphical_data_generator(kegg_collection[[x]])
    
    if(use_old_images) {
      img_path <- system.file("extdata", "old_imgs", paste0(gsub("path:", "", kegg_collection[[x]]$attrs$name), ".png"), package="psf")
    } else {
      img_path <- system.file("extdata", "pathway_imgs", paste0(gsub("path:", "", kegg_collection[[x]]$attrs$name), ".png"), package="psf")
    }
    
    pathway_img_new <- kegg_designer(group_graphics = kegg_collection[[x]]$group_nodes, node_graphics = graphical_data$node_coords, 
                                     pathway_image = magick::image_read(img_path),
                                     psf_output = psf_output, color_bar_psf_mode = TRUE, 
                                     col_legend_title = "Log PSF value", plot_type = "boxplot")
    
    heatmap_img <- magick::image_graph(width = 1000, height = 800, res = 96)
    
    if(ncol(exp_matrix) > 1) {
      gplots::heatmap.2(psf_output$sink_values_all, col = c(pal1(10), pal2(10)), trace = "none",
                        Rowv = FALSE, Colv = FALSE, ylab = "Sink values", main = paste0(x, "Sink PSV values"),
                        margin = c(10, 8), keysize = 1, key.title = "PSF log value",
                        dendrogram = 'none')
    }
    
    dev.off()
    
    c(pathway_img_new, heatmap_img)
  })
  
  magick::image_write(Reduce(c, plot_list), format = "pdf", file_name, quality = 300)
  
}
