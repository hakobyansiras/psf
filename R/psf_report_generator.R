#' Calculates psf for given kegg pathway based on expression matrix and generates pdf report with colored pathways and plots
#' @param kegg_collection list of kegg pathways
#' @param exp_matrix expression fold change matrix with gene entrez id rownames
#' @param folder_name name of the folder where pdf report(s) will be generated
#' @param use_old_images use olde kegg images(for use with curated pathway collection)
#' @param calculate_significance logical, if true then psf function will also calculate significance for the PSF values by shuffling all the network nodes and checking if the resulted PSF values were calculated by chance. When set to true volcano plot will be generetaed in pdf report.
#' @param coldata table of sample information where first column of the table corresponds to colnames of exp_matrix and the second column is a group information.
#' @param adj two values in between 0 and 1 which specify the x and y adjustment of the node labels, with 0 for left/bottom, 1 for right/top, and 0.5 for centered. Default value is c(0.48, 1). On most devices values outside 0 and 1 will also work.
#' @import gplots
#' @import ggplot2
#' @import ggrepel
#' @import magick
#' @import grDevices
#' @export
calc_psf_and_generate_report_from_collection <- function(kegg_collection, exp_matrix, folder_name, calculate_significance = FALSE, use_old_images = F, coldata = NULL, adj = c(0.48, 1)) {
  
  dir.create(folder_name)
  
  lapply(names(kegg_collection), function(x) {
    
    graphical_data <- graphical_data_generator(kegg_collection[[x]])
    
    if(use_old_images) {
      img_path <- system.file("extdata", "old_imgs", paste0(gsub("path:", "", kegg_collection[[x]]$attrs$name), ".png"), package="psf")
    } else {
      img_path <- system.file("extdata", "pathway_imgs", paste0(gsub("path:", "", kegg_collection[[x]]$attrs$name), ".png"), package="psf")
    }
    
    if(is.null(coldata)) {
      
      psf_output <- psf_signal_calculator_and_coloring(entrez_fc = exp_matrix, 
                                                       pathway = kegg_collection[[x]], 
                                                       pathway_name = x, 
                                                       calculate_significance = calculate_significance, 
                                                       no_color_mode = F)
      
      pathway_img_new <- kegg_designer(group_graphics = kegg_collection[[x]]$group_nodes, node_graphics = graphical_data$node_coords, 
                                       pathway_image = magick::image_read(img_path),
                                       psf_output = psf_output, color_bar_psf_mode = TRUE, 
                                       col_legend_title = "Log PSF value", plot_type = "boxplot", adj = adj)
      
      plots <- magick::image_graph(width = 1000, height = 800, res = 96)
      
      if(ncol(exp_matrix) > 1) {
        gplots::heatmap.2(psf_output$sink_values_all, col = c(pal1(10), pal2(10)), trace = "none",
                          Rowv = FALSE, Colv = FALSE, ylab = "Sink values", main = paste0(x, "Sink PSF values"),
                          margin = c(10, 8), keysize = 1, key.title = "PSF log value",
                          dendrogram = 'none')
      }
      
      dev.off()
      
      if(calculate_significance) {
        ### generating volcano plot
        sink_names <- unname(unlist(graph::nodeData(psf_output$psf_graph[[1]][[x]]$graph, psf_output$psf_graph[[1]][[x]]$sink.nodes, attr = "label")))
        
        volcano_table <- Reduce(rbind,
                                lapply(names(psf_output$psf_graph), function(y) {
                                  
                                  data.frame(log_psf = log(psf_output$psf_graph[[y]][[x]]$signal.at.sink), 
                                             pvalue = psf_output$psf_graph[[y]][[x]]$p.values + 0.00001, 
                                             sink_names = paste0(y, "_", sink_names),
                                             stringsAsFactors = F
                                  )
                                  
                                  
                                })
        )
        
        ## detecting up and down regulated significan sink psf values
        volcano_table$diffpsf <- "NO"
        volcano_table$diffpsf[volcano_table$log_psf > 0.6 & volcano_table$pvalue < 0.05] <- "UP"
        volcano_table$diffpsf[volcano_table$log_psf < -0.6 & volcano_table$pvalue < 0.05] <- "DOWN"
        
        ## adding labels for signifcant sinks
        volcano_table$psflabel <- ""
        volcano_table$psflabel[volcano_table$diffpsf != "NO"] <- volcano_table$sink_names[volcano_table$diffpsf != "NO"]
        
        
        volcano_plot <- ggplot(data=volcano_table, aes(x=log_psf, y=-log10(pvalue), col=diffpsf, label=psflabel)) +
          geom_point() + 
          theme_minimal() +
          ggtitle("PSF with p values") +
          theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5)) +
          geom_text_repel(max.overlaps = Inf) +
          scale_color_manual(values=c("blue", "black", "red")) +
          geom_vline(xintercept=c(-0.6, 0.6), col="red") +
          geom_hline(yintercept=-log10(0.05), col="red")
        
        ggsave("volcano_plot.png", plot = volcano_plot, device = "png", path = NULL,
               scale = 1, width = 1000, height = 800, units = "px",
               dpi = 96, limitsize = TRUE)
        
        volcano_plot_img <- image_read('volcano_plot.png')
        
        file.remove('volcano_plot.png')
        
        plots <- image_append(c(plots, volcano_plot_img))
      }
      
      magick::image_write(c(pathway_img_new, plots), format = "pdf", path = paste0(folder_name, "/", x, ".pdf"), quality = 300) 
      
      
    } else {
      
      exp_by_groups <- lapply(unique(coldata[, 2]), function(x) {
        samples <- coldata[which(coldata[,2] %in% x),1]
        exp_matrix[,samples]
      })
      names(exp_by_groups) <- unique(coldata[, 2])
      
      multiple_group_plots <- lapply(names(exp_by_groups), function(y) {
        exp_matrix <- exp_by_groups[[y]]
        
        
        psf_output <- psf_signal_calculator_and_coloring(entrez_fc = exp_matrix, 
                                                         pathway = kegg_collection[[x]], 
                                                         pathway_name = x, 
                                                         calculate_significance = calculate_significance, 
                                                         no_color_mode = F)
        
        pathway_img_new <- kegg_designer(group_graphics = kegg_collection[[x]]$group_nodes, node_graphics = graphical_data$node_coords, 
                                         pathway_image = magick::image_read(img_path),
                                         psf_output = psf_output, color_bar_psf_mode = TRUE, 
                                         col_legend_title = "Log PSF value", plot_type = "boxplot", adj = adj)
        
        pathway_img_new <- image_annotate(pathway_img_new, y, size = 70)
        
        plots <- magick::image_graph(width = 1000, height = 800, res = 96)
        
        if(ncol(exp_matrix) > 1) {
          gplots::heatmap.2(psf_output$sink_values_all, col = c(pal1(10), pal2(10)), trace = "none",
                            Rowv = FALSE, Colv = FALSE, ylab = "Sink values", main = paste0(x, "Sink PSV values"),
                            margin = c(10, 8), keysize = 1, key.title = "PSF log value",
                            dendrogram = 'none')
        }
        
        dev.off()
        
        if(calculate_significance) {
          ### generating volcano plot
          sink_names <- unname(unlist(graph::nodeData(psf_output$psf_graph[[1]][[x]]$graph, psf_output$psf_graph[[1]][[x]]$sink.nodes, attr = "label")))
          
          volcano_table <- Reduce(rbind,
                                  lapply(names(psf_output$psf_graph), function(y) {
                                    
                                    data.frame(log_psf = log(psf_output$psf_graph[[y]][[x]]$signal.at.sink), 
                                               pvalue = psf_output$psf_graph[[y]][[x]]$p.values + 0.00001, 
                                               sink_names = paste0(y, "_", sink_names),
                                               stringsAsFactors = F
                                    )
                                    
                                    
                                  })
          )
          
          ## detecting up and down regulated significan sink psf values
          volcano_table$diffpsf <- "NO"
          volcano_table$diffpsf[volcano_table$log_psf > 0.6 & volcano_table$pvalue < 0.05] <- "UP"
          volcano_table$diffpsf[volcano_table$log_psf < -0.6 & volcano_table$pvalue < 0.05] <- "DOWN"
          
          ## adding labels for signifcant sinks
          volcano_table$psflabel <- ""
          volcano_table$psflabel[volcano_table$diffpsf != "NO"] <- volcano_table$sink_names[volcano_table$diffpsf != "NO"]
          
          
          volcano_plot <- ggplot(data=volcano_table, aes(x=log_psf, y=-log10(pvalue), col=diffpsf, label=psflabel)) +
            geom_point() + 
            theme_minimal() +
            ggtitle("PSF with p values") +
            theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5)) +
            geom_text_repel(max.overlaps = Inf) +
            scale_color_manual(values=c("blue", "black", "red")) +
            geom_vline(xintercept=c(-0.6, 0.6), col="red") +
            geom_hline(yintercept=-log10(0.05), col="red")
          
          ggsave("volcano_plot.png", plot = volcano_plot, device = "png", path = NULL,
                 scale = 1, width = 1000, height = 800, units = "px",
                 dpi = 96, limitsize = TRUE)
          
          volcano_plot_img <- image_read('volcano_plot.png')
          
          file.remove('volcano_plot.png')
          
          plots <- image_append(c(plots, volcano_plot_img))
        }
        
        c(pathway_img_new, plots)
        
      })
      
      magick::image_write(Reduce(c, multiple_group_plots), format = "pdf", path = paste0(folder_name, "/", x, ".pdf"), quality = 300) 
      
    }
    
  })
  
}