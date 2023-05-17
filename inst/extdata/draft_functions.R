



#### build pathway node coloring function without graphnel object (get only psf results from graphnel and proceed with data frames)



plot_kegg_pathway <- function(graphnel_df, group_graphics, pathway_image, 
                              node_colors = NULL, custom_edge_mapping = FALSE, edge_mapping = FALSE, edge_in_mode = FALSE, highlight_nodes = NULL, highlight_color = "red",
                              col_legend_title = NULL, color_bar_lims, y_adj_sink = 3, node_fill_opacity = 0, present_node_modifications = FALSE, removed_nodes = NULL
                              ) {
  
  ### checking highlight nodes and their border colors
  if(length(highlight_color) > 1) {
    if(length(highlight_color) == length(highlight_nodes)) {
      highlight_color_vector <- setNames(object = highlight_color, nm = highlight_nodes)
    } else {
      stop("Error: highlighted nodes and their colors must be in the same length")
    }
  } else {
    highlight_color_vector <- setNames(object = rep(highlight_color, length(highlight_nodes)), nm = highlight_nodes)
  }
  
  img <- magick::image_draw(pathway_image)
  
  node_graphics <- graphnel_df$node_table
  
  ### node exp coloring
  if(any(node_graphics$node_id %in% node_colors$node_id)) {
    coloring_set <- node_graphics[node_colors$node_id,]
    graphics::rect(coloring_set$x_start,
                   coloring_set$y_start - 1,
                   coloring_set$x_end,
                   coloring_set$y_end,
                   # border = coloring_set$border_color,
                   border = NA,
                   # lty = coloring_set$lty_type,
                   lwd=2,
                   col = grDevices::adjustcolor(node_colors$col, alpha.f = 1)
    )
    
    graphics::text(x = coloring_set$x,
                   y = coloring_set$y_start + y_adj_text,
                   labels = coloring_set$label,
                   col = node_colors$text_col, adj = c(0,0.2) + c(0.48, 1))
    
  }
  
  ### adding sink node labels
  graphics::text(x = node_graphics[node_graphics$sink,"x_end"] + 10,
                 y = node_graphics[node_graphics$sink,"y"] - 30 + y_adj_sink, cex = 3,
                 labels = rep("*", sum(node_graphics$sink)),
                 col = rep("#9ACD32", sum(node_graphics$sink)), adj = c(0,0.2) + c(0.48, 1))
  
  if(!is.null(highlight_nodes)) {
    highlight_set <- node_graphics[which(node_graphics[,"node_id"] %in% highlight_nodes),]
    
    rect( highlight_set$x_start, 
          highlight_set$y_start, 
          highlight_set$x_end, 
          highlight_set$y_end, 
          border = highlight_color_vector[highlight_set$node_id], lty = "solid", lwd=2, col = adjustcolor( "#a3297a", alpha.f = node_fill_opacity))
  }
  
  
  ### scale color bar
  if(!is.null(node_colors)) {
    color_legend_maker(x = magick::image_info(img)$width - 230, y = 50, leg = 200, cols = c(pal1(10), pal2(10)), title = col_legend_title, lims = color_bar_lims, digits=3, prompt=FALSE,
                       lwd=4, outline=TRUE, subtitle = "", fsize = 1.3)
  }
  
  text(x = c(magick::image_info(img)$width - 88, magick::image_info(img)$width - 30),
       y = c(70, 65 + y_adj_sink - 6), cex = c(1.5, 3),
       labels = c("Sink node", "*"),
       col = c("#000000", "#9ACD32"), adj = c(0,0.2) + c(0.48, 1))
  
  
  ## color grop nodes
  if(length(group_graphics) > 0 ) {
    lapply(group_graphics, function(z) {
      rect( z$kegg.gr.x-z$kegg.gr.width*0.5, 
            z$kegg.gr.y+z$kegg.gr.height*0.5, 
            z$kegg.gr.x+z$kegg.gr.width*0.5, 
            z$kegg.gr.y-z$kegg.gr.height*0.5, 
            border = "yellow", lty = "dashed", lwd=2)
    })
  }
  
  if(present_node_modifications) {
    added_nodes <- node_graphics[which(node_graphics$existence == "added"),]
    if(nrow(added_nodes) > 0) {
      graphics::rect(added_nodes$x_start,
                     added_nodes$y_start - 1,
                     added_nodes$x_end,
                     added_nodes$y_end,
                     # border = coloring_set$border_color,
                     border = rep("black", nrow(added_nodes)),
                     # lty = coloring_set$lty_type,
                     lwd=1,
                     col = rep("#BFFFBF", nrow(added_nodes))
      )
      
      graphics::text(x = added_nodes$x,
                     y = added_nodes$y_start + y_adj_text,
                     labels = added_nodes$label,
                     col = rep("black", nrow(added_nodes)), adj = c(0,0.2) + c(0.48, 1))
    }
    
    if(!is.null(removed_nodes)) {
      graphics::rect(removed_nodes$x_start,
                     removed_nodes$y_start - 1,
                     removed_nodes$x_end,
                     removed_nodes$y_end,
                     # border = coloring_set$border_color,
                     border = rep("white", nrow(removed_nodes)),
                     # lty = coloring_set$lty_type,
                     lwd=2,
                     col = rep("white", nrow(removed_nodes))
      )
    }
    
  }
  
  #### continue from here
  if(edge_mapping) {
    
    if(custom_edge_mapping) {
      
      from <- graphnel_df$node_table[which(graphnel_df$node_table$node_id == highlight_nodes[1]),]
      to <- graphnel_df$node_table[which(graphnel_df$node_table$node_id == highlight_nodes[2]),]
      
      shape::Arrows(
        x0 = (from$x_start + from$x_end)/2,
        x1 = (to$x_start + to$x_end)/2,
        y0 = (from$y_start + from$y_end)/2,
        y1 = (to$y_start + to$y_end)/2,
        col = "red", 
        lwd=2, arr.length = 0.2, arr.type = "simple"
      )
        
    } else {
      if(!is.null(graphnel_df$edge_table)) {
        
        if(length(highlight_nodes) == 0) {
          shape::Arrows(
            x0 = graphnel_df$node_table[graphnel_df$edge_table$from,"x"],
            x1 = graphnel_df$node_table[graphnel_df$edge_table$to,"x"],
            y0 = graphnel_df$node_table[graphnel_df$edge_table$from,"y"],
            y1 = graphnel_df$node_table[graphnel_df$edge_table$to,"y"],
            col = graphnel_df$edge_table$color, 
            lty = ifelse(is.na(graphnel_df$edge_table$dashes), FALSE, graphnel_df$edge_table$dashes) + 1,
            lwd=2, arr.length = 0.2, arr.type = "simple"
          )
        } else {
          if(edge_in_mode) {
            index <- unique(which(graphnel_df$edge_table$from %in% highlight_nodes & graphnel_df$edge_table$to %in% highlight_nodes))
          } else {
            index <- unique(c(which(graphnel_df$edge_table$from %in% highlight_nodes), 
                              which(graphnel_df$edge_table$to %in% highlight_nodes)))
          }
          
          shape::Arrows(
            x0 = graphnel_df$node_table[graphnel_df$edge_table$from[index],"x"],
            x1 = graphnel_df$node_table[graphnel_df$edge_table$to[index],"x"],
            y0 = graphnel_df$node_table[graphnel_df$edge_table$from[index],"y"],
            y1 = graphnel_df$node_table[graphnel_df$edge_table$to[index],"y"],
            col = graphnel_df$edge_table$color[index], 
            lty = ifelse(is.na(graphnel_df$edge_table$dashes[index]), FALSE, graphnel_df$edge_table$dashes[index]) + 1,
            lwd=2, arr.length = 0.2, arr.type = "simple"
          )
        }
      }
    }
    
  }
  
  dev.off()
  return(img)
  
}

pathway_image = magick::image_read(system.file("extdata", "old_imgs", paste0(gsub("path:", "", edited_pathways_new$Chemokine_signaling_pathway$attrs$name), ".png"), package="psf"))

plot_kegg_pathway(graphnel_df = pathway_tables, group_graphics = edited_pathways_new$Chemokine_signaling_pathway$group_nodes, pathway_image = pathway_image, 
                              node_colors = NULL, edge_mapping = TRUE, edge_in_mode = TRUE, highlight_nodes = c("6", "20"), highlight_color = "red",
                              col_legend_title = NULL, color_bar_lims)


#### keep for future 
node_color_generator <- function(pathway, sample_id = "mean", log_norm = TRUE, color_nodes = "psf_activities") {
  if("exp_fc" %in% names(pathway)) {
    mapping_data = TRUE
    if(sample_id == "mean") {
      pathway_exp_values <- rowMeans(pathway$exp_fc)
      pathway_psf_values <- rowMeans(pathway$psf_activities)
    } else {
      pathway_exp_values <- pathway$exp_fc[,sample_id]
      pathway_psf_values <- pathway$psf_activities[,sample_id]
    }
    
    ordered_nodes <- rownames(pathway$exp_fc)[order(drop(pathway_exp_values))]
    
    if(log_norm) {
      pathway_exp_values <- round(log(drop(pathway_exp_values)[order(drop(pathway_exp_values))] + 0.00001), digits = 5)
      pathway_psf_values <- round(log(drop(pathway_psf_values)[order(drop(pathway_psf_values))] + 0.00001), digits = 5)
    } else {
      pathway_exp_values <- round(drop(pathway_exp_values)[order(drop(pathway_exp_values))], digits = 5)
      pathway_psf_values <- round(drop(pathway_psf_values)[order(drop(pathway_psf_values))], digits = 5)
    }
  } else {
    mapping_data = FALSE
  }
    
    if(mapping_data) {
      if(color_nodes == "psf_activities") {
        pathway_node_values <- pathway_psf_values
      } else {
        pathway_node_values <- pathway_exp_values
      }
      
      node_colors <- color_code(values = pathway_node_values, pal1 = pal1, pal2 = pal2, log_scale = log_norm)
      
      node_colors <- data.frame(node_id = names(pathway_node_values)[c(which(pathway_node_values <= 0), which(pathway_node_values > 0))], 
                                col = node_colors,
                                text_col = unname(sapply(node_colors, function(x) {c( "black", "white")[  1+(sum( col2rgb(x) *c(299, 587,114))/1000 < 123) ]})),
                                stringsAsFactors = F
      )
      
      node_colors <- node_colors[which(node_colors$node_id %in% graph::nodes(pathway$graph)[which(unlist(graph::nodeData(pathway$graph, attr = "type")) != "map")]),]
      
    } else {
      stop("Please provide pathway with evalueated activity")
    }  
    
  if(log_norm) {
    hover_text <- paste(paste("Log exp FC", pathway_exp_values),
                        paste("Log PSF", pathway_psf_values),
                        sep = "<br>")
  } else {
    hover_text <- paste(paste("Exp FC", pathway_exp_values),
                        paste("PSF", pathway_psf_values),
                        sep = "<br>")
  }
  
  names(hover_text) <- ordered_nodes
  
    if(color_nodes == "psf_activities") {
      col_legend_title = ifelse(log_norm, "Log PSF value", "PSF value")
    } else {
      col_legend_title = ifelse(log_norm, "Log FC value", "FC value")
    }
  
  return(list(node_colors = node_colors, col_legend_title = col_legend_title, color_bar_lims = range(pathway_node_values), hover_text = hover_text))
}









plot_pathway <- function(pathway, plot_type = "visnet", edge_mapping = FALSE,
                         color_nodes = NULL, sample_id = "mean", log_norm = T, 
                         highlight_nodes = NULL, highlight_color = "red", plot_sink_values = F,
                         y_adj_text = 0, y_adj_sink = 0, use_old_images = F) {
  
  ### checking highlight nodes and their border colors
  if(length(highlight_color) > 1) {
    if(length(highlight_color) == length(highlight_nodes)) {
      highlight_color_vector <- setNames(object = highlight_color, nm = highlight_nodes)
    } else {
      stop("Error: highlighted nodes and their colors must be in the same length")
    }
  } else {
    highlight_color_vector <- setNames(object = rep(highlight_color, length(highlight_nodes)), nm = highlight_nodes)
  }
  
  
  if("exp_fc" %in% names(pathway)) {
    mapping_data = TRUE
    if(sample_id == "mean") {
      pathway_exp_values <- rowMeans(pathway$exp_fc)
      pathway_psf_values <- rowMeans(pathway$psf_activities)
    } else {
      pathway_exp_values <- pathway$exp_fc[,sample_id]
      pathway_psf_values <- pathway$psf_activities[,sample_id]
    }
    
    if(log_norm) {
      pathway_exp_values <- round(log(drop(pathway_exp_values)[order(drop(pathway_exp_values))] + 0.00001), digits = 5)
      pathway_psf_values <- round(log(drop(pathway_psf_values)[order(drop(pathway_psf_values))] + 0.00001), digits = 5)
    } else {
      pathway_exp_values <- round(drop(pathway_exp_values)[order(drop(pathway_exp_values))], digits = 5)
      pathway_psf_values <- round(drop(pathway_psf_values)[order(drop(pathway_psf_values))], digits = 5)
    }
  } else {
    mapping_data = FALSE
  }
  
  if(!is.null(color_nodes)) {
    
    if(mapping_data) {
      if(color_nodes == "psf_activities") {
        pathway_node_values <- pathway_psf_values
      } else {
        pathway_node_values <- pathway_exp_values
      }
      
      node_colors <- color_code(values = pathway_node_values, pal1 = pal1, pal2 = pal2, log_scale = log_norm)
      
      node_colors <- data.frame(node_id = names(pathway_node_values)[c(which(pathway_node_values <= 0), which(pathway_node_values > 0))], 
                                col = node_colors,
                                text_col = unname(sapply(node_colors, function(x) {c( "black", "white")[  1+(sum( col2rgb(x) *c(299, 587,114))/1000 < 123) ]})),
                                stringsAsFactors = F
      )
      
      node_colors <- node_colors[which(node_colors$node_id %in% graph::nodes(pathway$graph)[which(unlist(graph::nodeData(pathway$graph, attr = "type")) != "map")]),]
      
    } else {
      stop("Please provide pathway with evalueated activity")
    }  
    
    if(color_nodes == "psf_activities") {
      col_legend_title = ifelse(log_norm, "Log PSF value", "PSF value")
    } else {
      col_legend_title = ifelse(log_norm, "Log FC value", "FC value")
    }
    
  } else {
    node_colors <- NULL
  }
  
  pathway_tables <- graphnel_to_df(pathway, extended = TRUE)
  
  node_graphics <- pathway_tables$node_table
  
  if(plot_type == "kegg") {
    
    if(use_old_images) {
      img_path <- system.file("extdata", "old_imgs", paste0(gsub("path:", "", pathway$attrs$name), ".png"), package="psf")
    } else {
      if(file.exists(system.file("extdata", "pathway_imgs", paste0(gsub("path:", "", pathway$attrs$name), ".png"), package="psf"))) {
        img_path <- system.file("extdata", "pathway_imgs", paste0(gsub("path:", "", pathway$attrs$name), ".png"), package="psf")
      } else {
        # image <- paste0("http://rest.kegg.jp/get/", pathway_id, "/image")
        
        img_path <- tempfile()
        download.file(pathway$attrs$image, destfile = img_path, method = "curl")
      }
      
    }
    
    pathway_image = magick::image_read(img_path)
    
    img <- magick::image_draw(pathway_image)
    
    ### node exp coloring
    rownames(node_colors) <- node_colors$node_id
    
    rownames(node_graphics) <- node_graphics$node_id
    
    sink_node_graphics <- node_graphics[which(node_graphics$sink),]
    
    
    if(any(node_graphics$node_id %in% node_colors$node_id)) {
      coloring_set <- node_graphics[node_colors$node_id,]
      graphics::rect(coloring_set$x_start,
                     coloring_set$y_start - 1,
                     coloring_set$x_end,
                     coloring_set$y_end,
                     # border = coloring_set$border_color,
                     border = NA,
                     # lty = coloring_set$lty_type,
                     lwd=2,
                     col = grDevices::adjustcolor(node_colors$col, alpha.f = 1)
      )
      
      graphics::text(x = coloring_set$x,
                     y = coloring_set$y_start + y_adj_text,
                     labels = coloring_set$label,
                     col = node_colors$text_col, adj = c(0,0.2) + c(0.48, 1))
      
    }
    
    ### adding sink node labels
    graphics::text(x = node_graphics[node_graphics$sink,"x_end"] + 10,
                   y = node_graphics[node_graphics$sink,"y"] - 30 + y_adj_sink, cex = 3,
                   labels = rep("*", sum(node_graphics$sink)),
                   col = rep("#9ACD32", sum(node_graphics$sink)), adj = c(0,0.2) + c(0.48, 1))
    
    if(!is.null(highlight_nodes)) {
      highlight_set <- node_graphics[which(node_graphics[,"node_id"] %in% highlight_nodes),]
      
      rect( highlight_set$x_start, 
            highlight_set$y_start, 
            highlight_set$x_end, 
            highlight_set$y_end, 
            border = highlight_color_vector[highlight_set$node_id], lty = "solid", lwd=2, col = adjustcolor( "#a3297a", alpha.f = 0))
    }
    
    
    ### scale color bar
    if(!is.null(color_nodes)) {
      color_legend_maker(x = magick::image_info(img)$width - 230, y = 50, leg = 200, cols = c(pal1(10), pal2(10)), title = col_legend_title, lims = range(pathway_node_values), digits=3, prompt=FALSE,
                         lwd=4, outline=TRUE, subtitle = "", fsize = 1.3)
    }
    
    text(x = c(magick::image_info(img)$width - 88, magick::image_info(img)$width - 30),
         y = c(70, 65 + y_adj_sink - 6), cex = c(1.5, 3),
         labels = c("Sink node", "*"),
         col = c("#000000", "#9ACD32"), adj = c(0,0.2) + c(0.48, 1))
    
    
    ## color group nodes
    if(length(pathway$group_nodes) > 0 ) {
      lapply(pathway$group_nodes, function(z) {
        graphics::rect( z$kegg.gr.x-z$kegg.gr.width*0.5, 
                        z$kegg.gr.y+z$kegg.gr.height*0.5, 
                        z$kegg.gr.x+z$kegg.gr.width*0.5, 
                        z$kegg.gr.y-z$kegg.gr.height*0.5, 
                        border = "yellow", lty = "dashed", lwd=2)
      })
    }
    
    
    dev.off()
    
    if(plot_sink_values) {
      if(mapping_data) {
        
        sink_names <- unlist(graph::nodeData(pathway$graph, pathway$sink.nodes, attr = "label"))
        
        sink_order_by_coord <- order(as.integer(unname(unlist(graph::nodeData(pathway$graph, pathway$sink.nodes, attr = "kegg.gr.y")))), decreasing = T)
        
        repeating_sinks <- table(sink_names)[which(table(sink_names) > 1)]
        
        if(length(repeating_sinks) > 0) {
          for(i in 1:length(repeating_sinks)) {
            
            sink_names[which(sink_names == names(repeating_sinks[i]))] <- sapply(1:repeating_sinks[i], function(x) {paste0(paste(rep(" ", x - 1), collapse = ""), names(repeating_sinks[i]))})
            
          }
        }
        
        psf_activities_vec <- as.vector(pathway$psf_activities[pathway$sink.nodes,])
        
        names(psf_activities_vec) <- rep(pathway$sink.nodes, ncol(pathway$psf_activities))
        
        if(log_norm) {
          psf_activities_vec <- log(psf_activities_vec)
        }
        
        psf_activities_vec <- psf_activities_vec[order(psf_activities_vec)]
        
        
        sink_signals <- data.frame(sink_name = sink_names[names(psf_activities_vec)], 
                                   signal = psf_activities_vec, 
                                   dot_color = color_code(values = psf_activities_vec, pal1 = pal1, pal2 = pal2, log_scale = log_norm),
                                   stringsAsFactors = F)
        
        sink_signals$sink_name <- factor(sink_signals$sink_name, levels = sink_names[sink_order_by_coord])
        
        if(ncol(pathway$exp_fc) == 1) {
          
          sink_plot <- ggplot(sink_signals, aes(x=sink_name, y=signal)) +
            geom_bar(color="black", lwd=0.2, stat="identity") +
            coord_flip() +
            theme_bw()
          
        } else {
          sink_plot <- ggplot(sink_signals, aes(x=sink_name, y=signal, fill=signal)) +
            geom_boxplot(color="black", lwd=0.2, outlier.shape=NA) +
            geom_point(color = sink_signals$dot_color) +
            geom_jitter(color = sink_signals$dot_color, width = 0.2) +
            guides(fill=FALSE) +
            coord_flip() +
            theme_bw()
        }
        
        ggsave("sink_plot.png", plot = sink_plot, device = "png", path = NULL,
               scale = 1, width = 400, height = magick::image_info(img)$height, units = "px",
               dpi = 96, limitsize = TRUE)
        
        sink_plot_img <- image_read('sink_plot.png')
        
        file.remove('sink_plot.png')
        
        combined_img <- image_append(c(img, sink_plot_img))
        
        return(combined_img)
        
      } else {
        stop("Please provide pathway with evalueated activity")
      }
      
    } else {
      return(img)
    }
    
  } else {
    
    node_table <- pathway_tables$node_table[,c("node_id", "label", "vis_width", "font.color", "size", "font.size", "borderWidth", "color.border", "color.background", "shape", "title", "x", "y")]
    
    colnames(node_table) <- c("id", "label", "widthConstraint.maximum", "font.color", "size", "font.size", "borderWidth", "color.border", "color.background", "shape", "title", "x", "y")
    
    edge_table <- pathway_tables$edge_table[,c("from", "to", "color", "arrows.to.enabled", "arrows.to.type", "label", "dashes")]
    
    if(!is.null(node_colors)) {
      node_table$color.background <- unname(sapply(node_table$id, function(x) {
        if(x %in% node_colors$node_id) {
          node_colors[which(node_colors$node_id == x),"col"]
        } else {
          "#BFFFBF"
        }
      }))
      
      node_table$font.color <- unname(sapply(node_table$id, function(x) {
        if(x %in% node_colors$node_id) {
          node_colors[which(node_colors$node_id == x),"text_col"]
        } else {
          "#000000"
        }
      }))
      
      if(log_norm) {
        node_table$title <- paste(node_table$title,
                                  paste("Log exp FC", pathway_exp_values[node_table$id]),
                                  paste("Log PSF", pathway_psf_values[node_table$id]),
                                  sep = "<br>")
      } else {
        node_table$title <- paste(node_table$title,
                                  paste("Exp FC", pathway_exp_values[node_table$id]),
                                  paste("PSF", pathway_psf_values[node_table$id]),
                                  sep = "<br>")
      }
      
    } else {
      
      if(mapping_data) {
        node_table$title <- paste(node_table$title,
                                  paste("Exp FC", pathway_exp_values[node_table$id]),
                                  paste("PSF", pathway_psf_values[node_table$id]),
                                  sep = "<br>")
      }
      
    }
    
    if(!is.null(highlight_nodes)) {
      
      node_table$color.border <- unname(sapply(node_table$id, function(x) {
        if(x %in% highlight_nodes) {
          if(length(highlight_color) > 1) {
            highlight_color_vector[x]
          } else {
            highlight_color
          }
        } else {
          "#BFFFBF"
        }
      }))
      
    }
    
    node_table$image <- "unselected"
    
    if(!is.null(node_colors)) {
      
      magick::autoviewer_disable() ### to avoid legend plotting before netowrk rendering
      legend_img <- magick::image_device(width = 480, height = 480)
      plot.new()
      
      color_legend_maker(x = 0.05, y = 0, leg = 0.9, cols = c(pal1(10), pal2(10)), title = col_legend_title, lims = range(pathway_node_values), digits=3, prompt=FALSE,
                         lwd=4, outline=TRUE, subtitle = "", fsize = 1.3)
      
      temp_legend <- tempfile()
      
      magick::image_write(magick::image_trim(legend_img, fuzz = 0), path = temp_legend)
      
      legend_path <- paste('data:image/png;base64', RCurl::base64Encode(readBin(temp_legend, 'raw', file.info(temp_legend)[1, 'size']), 'txt'), sep = ',')
      
      legend_data_frame <- data.frame(
        id = as.character(max(as.integer(node_table$id)) + 1),
        label = "Color legend", 
        widthConstraint.maximum = NA, 
        font.color = "#000000", 
        size = 30, 
        font.size = 32, 
        borderWidth = 0, 
        color.border = "",
        color.background = "",
        shape = "image", 
        title = "",
        x = max(graphical_data$node_coords$x_end),
        y = min(graphical_data$node_coords$y_end) - 10,
        image = legend_path
      )
      
      node_table <- rbind(node_table, legend_data_frame)
    }
    
    
    saving_pathway_name <- gsub("(", "", pathway$attrs$title, fixed = T)
    saving_pathway_name <- gsub(")", "", saving_pathway_name, fixed = T)
    saving_pathway_name <- gsub("-", "_", saving_pathway_name)
    saving_pathway_name <- gsub(" ", "_", saving_pathway_name)
    saving_pathway_name <- gsub("___", "_", saving_pathway_name)
    
    visNetwork(nodes = node_table, edges = edge_table, width = "100%", height = "800px") %>% 
      visIgraphLayout(layout = "layout_nicely") %>% 
      visExport(name = saving_pathway_name) %>%
      visInteraction(navigationButtons = TRUE, multiselect = T)
    
  }
  
}


pathway, plot_type = "visnet",
color_nodes = NULL, sample_id = "mean", log_norm = T, 
highlight_nodes = NULL, highlight_color = "red", plot_sink_values = F,
y_adj_text = 0, y_adj_sink = 0, use_old_images = F



kegg_node_mapper <- function(group_graphics, kegg_pathway_graphics, pathway_name, pathway_image, highlight.genes = NULL, 
                             edge_mapping = FALSE, advanced_edge_ampping = FALSE, show_changes = FALSE, 
                             edge_in_mode = FALSE, highlight_color = "green", opacity = 0, color.genes = NULL, 
                             color_bar_psf_mode = F, col_legend_title, color_bar_lims = NULL, draw_color_bar = F
) {
  
  img <- image_draw(pathway_image)
  
  if(is.null(color.genes)) {
    for(i in 1:nrow(kegg_pathway_graphics$node_coords)) {
      
      if(kegg_pathway_graphics$node_coords[i,"changed"]) {
        lty_type <- "dotted"
      }
      
      if(kegg_pathway_graphics$node_coords[i,"sink"]) {
        border_color <- "blue"
        lty_type <- "dashed"
      } else {
        border_color <- "red"
        lty_type <- "solid"
      }
      
      if(show_changes) {
        if(kegg_pathway_graphics$node_coords[i,"changed"]) {
          border_color <- "purple"
        }
        if(!kegg_pathway_graphics$node_coords[i,"exist"]) {
          border_color <- "red"
          lty_type <- "dashed"
        }
      } else {
        if(!kegg_pathway_graphics$node_coords[i,"exist"]) {
          border_color <- NA
        }
      }
      
      if(kegg_pathway_graphics$node_coords[i,"node_class"] == "bio_event") {
        
        plotrix::textbox(x = c(kegg_pathway_graphics$node_coords[i,"x_start"], kegg_pathway_graphics$node_coords[i,"x_end"]),
                         y = c(kegg_pathway_graphics$node_coords[i,"y_start"], kegg_pathway_graphics$node_coords[i,"y_end"]), 
                         textlist = kegg_pathway_graphics$node_coords[i,"gr_name"], 
                         fill = "#32CD32", lty = lty_type, lwd = 2, leading = 0.3,
                         justify = 'c', border = border_color, adj = c(0,0.25))
        
      } else {
        rect( kegg_pathway_graphics$node_coords[i,"x_start"], 
              kegg_pathway_graphics$node_coords[i,"y_start"], 
              kegg_pathway_graphics$node_coords[i,"x_end"], 
              kegg_pathway_graphics$node_coords[i,"y_end"], 
              border = border_color, lty = lty_type, lwd=2)
      }
      
      # text( kegg.pathway$pathway.info[[i]]$graphics$x, height - kegg.pathway$pathway.info[[i]]$graphics$y, kegg.pathway$pathway.info[[i]]$graphics$label, 
      #       cex = ifelse( nchar(kegg.pathway$pathway.info[[i]]$graphics$label)<=6,.3/2,.24/2), col = "black", family="sans" )
    }
  } else {
    
    ### node exp coloring
    
    node_graphics <- cbind(kegg_pathway_graphics$node_coords, t(
      sapply(kegg_pathway_graphics$node_coords$sink, function(x) {
        
        if(x) {
          c("blue", "dashed")
        } else {
          c("red", "solid")
        }
        
      })
    ))
    
    colnames(node_graphics)[15:16] <- c("border_color", "lty_type")
    
    node_graphics$border_color <- as.character(node_graphics$border_color)
    
    node_graphics$lty_type <- as.character(node_graphics$lty_type)
    
    node_graphics$x_center <- node_graphics$x_start + (node_graphics$x_end - node_graphics$x_start)/2
    
    
    rownames(color.genes) <- color.genes$node_id
    
    rownames(node_graphics) <- node_graphics$node_id
    
    rect( node_graphics$x_start, 
          node_graphics$y_start, 
          node_graphics$x_end, 
          node_graphics$y_end, 
          border = node_graphics$border_color, lty = node_graphics$lty_type, lwd=2)
    
    
    if(any(node_graphics$node_id %in% color.genes$node_id)) {
      coloring_set <- node_graphics[color.genes$node_id,]
      rect( coloring_set$x_start,
            coloring_set$y_start,
            coloring_set$x_end,
            coloring_set$y_end,
            border = coloring_set$border_color, lty = coloring_set$lty_type, lwd=2,
            col = adjustcolor( color.genes$col, alpha.f = 1)
      )
      
      text(x = coloring_set$x_center,
           y = coloring_set$y_start,
           labels = coloring_set$gr_name,
           col = color.genes$text_col, adj = c(0,0.2) + c(0.48, 1))
      
    }
    
  }
  
  if(any(kegg_pathway_graphics$node_coords[,"node_id"] %in% highlight.genes[,"node_id"])) {
    
    highlight_set <- kegg_pathway_graphics$node_coords[which(kegg_pathway_graphics$node_coords[,"node_id"] %in% highlight.genes$node_id),]
    
    rect( highlight_set$x_start, 
          highlight_set$y_start, 
          highlight_set$x_end, 
          highlight_set$y_end, 
          border = highlight_color, lty = "solid", lwd=2, col = adjustcolor( "#a3297a", alpha.f = opacity))
    
    
  }
  
  ### scale color bar
  
  if(draw_color_bar) {
    if(color_bar_psf_mode) {
      color_legend_maker(x = magick::image_info(img)$width - 230, y = 50, leg = 200, cols = c(pal1(10), pal2(10)), title = col_legend_title, lims = color_bar_lims, digits=3, prompt=FALSE,
                         lwd=4, outline=TRUE, subtitle = "", fsize = 1.3)
    } else {
      color_legend_maker(x = magick::image_info(img)$width - 230, y = 50, leg = 200, cols = exp_pal(20), title = col_legend_title, lims = exp_color_all$exp_lims, digits=3, prompt=FALSE,
                         lwd=4, outline=TRUE, subtitle = "", fsize = 1.3)
    }
  }
  
  
  ## color grop nodes
  if(length(group_graphics) > 0 ) {
    lapply(group_graphics, function(z) {
      rect( z$kegg.gr.x-z$kegg.gr.width*0.5, 
            z$kegg.gr.y+z$kegg.gr.height*0.5, 
            z$kegg.gr.x+z$kegg.gr.width*0.5, 
            z$kegg.gr.y-z$kegg.gr.height*0.5, 
            border = "yellow", lty = "dashed", lwd=2)
    })
  }
  
  if(edge_mapping) {
    if(advanced_edge_ampping) {
      
      from <- kegg_pathway_graphics$node_coords[which(kegg_pathway_graphics$node_coords$node_id == highlight.genes[1,"node_id"]),]
      to <- kegg_pathway_graphics$node_coords[which(kegg_pathway_graphics$node_coords$node_id == highlight.genes[2,"node_id"]),]
      
      shape::Arrows(
        x0 = (from$x_start + from$x_end)/2,
        x1 = (to$x_start + to$x_end)/2,
        y0 = (from$y_start + from$y_end)/2,
        y1 = (to$y_start + to$y_end)/2,
        col = "red", 
        lwd=2, arr.length = 0.2, arr.type = "simple"
      )
      
    } else {
      
      if(!is.null(kegg_pathway_graphics$edge_coords)) {
        
        if(!show_changes) {
          kegg_pathway_graphics$edge_coords <- kegg_pathway_graphics$edge_coords[which(kegg_pathway_graphics$edge_coords$lty == "solid"),]
          # kegg_pathway_graphics$edge_coords$lty <- "solid"
          kegg_pathway_graphics$edge_coords$col[which(kegg_pathway_graphics$edge_coords$col == "purple")] <- "red"
          kegg_pathway_graphics$edge_coords$col[which(kegg_pathway_graphics$edge_coords$col == "yellow")] <- "blue"
        }
        
        if(nrow(highlight.genes) == 0) {
          shape::Arrows(
            x0 = kegg_pathway_graphics$edge_coords$x0, 
            x1 = kegg_pathway_graphics$edge_coords$x1, 
            y0 = kegg_pathway_graphics$edge_coords$y0, 
            y1 = kegg_pathway_graphics$edge_coords$y1,
            col = kegg_pathway_graphics$edge_coords$col,
            lty = kegg_pathway_graphics$edge_coords$lty,
            lwd=2, arr.length = 0.2, arr.type = "simple"
          )
        } else {
          if(edge_in_mode) {
            index <- unique(which(kegg_pathway_graphics$edge_coords$from %in% highlight.genes[,"node_id"] & kegg_pathway_graphics$edge_coords$to %in% highlight.genes[,"node_id"]))
          } else {
            index <- unique(c(which(kegg_pathway_graphics$edge_coords$from %in% highlight.genes[,"node_id"]), 
                              which(kegg_pathway_graphics$edge_coords$to %in% highlight.genes[,"node_id"])))
          }
          
          shape::Arrows(
            x0 = kegg_pathway_graphics$edge_coords$x0[index], 
            x1 = kegg_pathway_graphics$edge_coords$x1[index], 
            y0 = kegg_pathway_graphics$edge_coords$y0[index], 
            y1 = kegg_pathway_graphics$edge_coords$y1[index],
            col = kegg_pathway_graphics$edge_coords$col[index], 
            lty = kegg_pathway_graphics$edge_coords$lty[index],
            lwd=2, arr.length = 0.2, arr.type = "simple"
          )
        }
      }
    }
  }
  
  
  
  dev.off()
  
  return(img)
  
}