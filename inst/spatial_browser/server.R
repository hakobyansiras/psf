library(data.table)
library(DT)
library(miniUI)
library(shiny)
library(Seurat)
library(ggplot2)
library(psf)
library(magick)
library(shinyjs)
library(visNetwork)
library(plotly)

load("melanoma_spatial_demo_data.RData")

psf_subset <- function(psf_collection, pathway_name, sample_list) {
  
  pathway <- psf_collection[[pathway_name]]
  
  pathway$psf_activities <- pathway$psf_activities[,sample_list, drop = F]
  pathway$exp_fc <- pathway$exp_fc[,sample_list, drop = F]
  
  return(pathway)
}

plot_kegg_pathway <- function(graphnel_df, group_graphics, pathway_image, 
                              node_colors = NULL, custom_edge_mapping = FALSE, edge_mapping = FALSE, edge_in_mode = FALSE, highlight_nodes = NULL, highlight_color = "red",
                              col_legend_title = NULL, color_bar_lims, y_adj_sink = 3, node_fill_opacity = 0, present_node_modifications = FALSE, removed_nodes = NULL, y_adj_text = 3, custom_color_scale = NULL
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
    if(is.null(custom_color_scale)) {
      color_legend_maker(x = magick::image_info(img)$width - 230, y = 50, leg = 200, cols = c(psf:::pal1(10), psf:::pal2(10)), title = col_legend_title, lims = color_bar_lims, digits=3, prompt=FALSE,
                         lwd=4, outline=TRUE, subtitle = "", fsize = 1.3)
    } else {
      color_legend_maker(x = magick::image_info(img)$width - 230, y = 50, leg = 200, cols = c(custom_color_scale[1:50], custom_color_scale[51:100]), title = col_legend_title, lims = color_bar_lims, digits=3, prompt=FALSE,
                         lwd=4, outline=TRUE, subtitle = "", fsize = 1.3)
    }
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


#### node color coding function ####
node_color_generator <- function(pathway, sample_id = "mean", log_norm = TRUE, color_nodes = "psf_activities", custom_color_scale = NULL) {
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
      hover_text <- paste(paste("Log exp FC", round(log(pathway_exp_values + 0.00001), digits = 5)),
                          paste("Log PSF", round(log(pathway_psf_values + 0.00001), digits = 5)),
                          sep = "<br>")
      node_values <- cbind(round(log(pathway_exp_values + 0.00001), digits = 5), round(log(pathway_psf_values + 0.00001), digits = 5))
    } else {
      hover_text <- paste(paste("Exp FC", round(pathway_exp_values, digits = 5)),
                          paste("PSF", round(pathway_psf_values, digits = 5)),
                          sep = "<br>")
      node_values <- cbind(round(pathway_exp_values, digits = 5), round(pathway_psf_values, digits = 5))
    }
    
    names(hover_text) <- rownames(pathway$exp_fc)
    
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
    
    if(is.null(custom_color_scale)) {
      node_colors <- color_code(values = pathway_node_values, pal1 = pal1, pal2 = pal2, log_scale = log_norm)
      
      node_colors <- data.frame(node_id = names(pathway_node_values)[c(which(pathway_node_values <= 0), which(pathway_node_values > 0))], 
                                col = node_colors,
                                text_col = unname(sapply(node_colors, function(x) {c( "black", "white")[  1+(sum( col2rgb(x) *c(299, 587,114))/1000 < 123) ]})),
                                stringsAsFactors = F
      )
    } else {
      node_colors <- map_colors(values = pathway_node_values, colors = custom_color_scale)
      
      node_colors <- data.frame(node_id = names(pathway_node_values), 
                                col = node_colors,
                                text_col = unname(sapply(node_colors, function(x) {c( "black", "white")[  1+(sum( col2rgb(x) *c(299, 587,114))/1000 < 123) ]})),
                                stringsAsFactors = F
      )
    }
    
    
    
    node_colors <- node_colors[which(node_colors$node_id %in% graph::nodes(pathway$graph)[which(unlist(graph::nodeData(pathway$graph, attr = "type")) != "map")]),]
    
  } else {
    stop("Please provide pathway with evalueated activity")
  }
  
  if(color_nodes == "psf_activities") {
    col_legend_title = ifelse(log_norm, "Log PSF value", "PSF value")
  } else {
    col_legend_title = ifelse(log_norm, "Log FC value", "FC value")
  }
  
  return(list(node_colors = node_colors, col_legend_title = col_legend_title, color_bar_lims = range(pathway_node_values), hover_text = hover_text, node_values = node_values))
}

map_colors <- function(values, colors) {
  middle_index <- length(colors) / 2
  negative_colors <- rev(colors[1:middle_index])
  positive_colors <- colors[(middle_index + 1):length(colors)]
  max_value <- max(abs(values))
  mapped_colors <- sapply(values, function(x) {
    if (x < 0) {
      # Map negative values to negative_colors
      index <- max(1, ceiling((abs(x) / max_value) * middle_index))
      return(negative_colors[index])
    } else if (x > 0) {
      # Map positive values to positive_colors
      index <- max(1, ceiling((x / max_value) * middle_index))
      return(positive_colors[index])
    } else {
      # Assign the middle color for zero values
      return(colors[middle_index])
    }
  })
  return(mapped_colors)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#### visnet data extractor ####
vis_extract <- function(vis_table, node_colors = NULL, higlight_node = NULL, custom_color_scale = NULL) {
  node_table <- vis_table$node_table[,c("node_id", "label", "vis_width", "font.color", "size", "font.size", "borderWidth", "color.border", "color.background", "shape", "title", "x", "y")]
  
  colnames(node_table) <- c("id", "label", "widthConstraint.maximum", "font.color", "size", "font.size", "borderWidth", "color.border", "color.background", "shape", "title", "x", "y")
  
  if(!is.null(node_colors)) {
    
    magick::autoviewer_disable() ### to avoid legend plotting before netowrk rendering
    legend_img <- magick::image_device(width = 480, height = 480)
    plot.new()
    
    if(is.null(custom_color_scale)) {
      color_legend_maker(x = 0.05, y = 0, leg = 0.9, cols = c(psf:::pal1(10), psf:::pal2(10)), title = node_colors$col_legend_title, lims = node_colors$color_bar_lims, digits=3, prompt=FALSE,
                               lwd=4, outline=TRUE, subtitle = "", fsize = 1.3)
    } else {
      color_legend_maker(x = 0.05, y = 0, leg = 0.9, cols = c(custom_color_scale[1:50], custom_color_scale[51:100]), title = node_colors$col_legend_title, lims = node_colors$color_bar_lims, digits=3, prompt=FALSE,
                               lwd=4, outline=TRUE, subtitle = "", fsize = 1.3)
    }
    
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
      x = max(node_table$x),
      y = min(node_table$y) - 10,
      image = legend_path
    )
    
    node_table$image <- "unselected"
    node_table <- rbind(node_table, legend_data_frame)
    
    node_table$color.background <- unname(sapply(node_table$id, function(x) {
      if(x %in% node_colors$node_colors$node_id) {
        node_colors$node_colors[which(node_colors$node_colors$node_id == x),"col"]
      } else {
        "#BFFFBF"
      }
    }))
    
    node_table$color.border <- node_table$color.background
    if(!is.null(higlight_node)) {
      node_table$color.border[which(node_table$id %in% higlight_node)] <- "#48f542"
    }
    
    
    node_table$font.color <- unname(sapply(node_table$id, function(x) {
      if(x %in% node_colors$node_colors$node_id) {
        node_colors$node_colors[which(node_colors$node_colors$node_id == x),"text_col"]
      } else {
        "#000000"
      }
    }))
    
    node_table$title <- paste(node_table$title, node_colors$hover_text[node_table$id], sep = "<br>")
  }
  
  
  edge_table <- vis_table$edge_table[,c("id", "from", "to", "color", "arrows.to.enabled", "arrows.to.type", "label", "dashes")]
  return(list(node_table = node_table, edge_table = edge_table))
}

sink_name_to_id <- setNames(object = rownames(kegg_sink_to_process), nm = paste0(kegg_sink_to_process$Pathway_name, "->", kegg_sink_to_process$Sink))

# Your existing code for data preparation
# Extract the spatial image (stored as an array)
spatial_image <- as.raster(psf_saptial_results$spatial_image)

# Extract the low-resolution coordinates using GetTissueCoordinates
coords <- psf_saptial_results$coords # GetTissueCoordinates(object = psf_spatial$spatial_psf_obj, scale = "lowres")

coords$imagerow <- nrow(spatial_image) - coords$imagerow

# Create a data frame for plotting
spot_groups <- psf_saptial_results$meta.data$seurat_clusters # psf_spatial$spatial_psf_obj@meta.data$seurat_clusters
# Define cluster colors
unique_clusters <- sort(unique(spot_groups))
cluster_colors <- setNames(gg_color_hue(length(unique_clusters)), unique_clusters)

plotting_data <- data.frame(coords, spot_groups = spot_groups)

# All points are highlighted initially
plotting_data$Highlight <- "Highlighted"
plotting_data$ColorCode <- cluster_colors[plotting_data$spot_groups]
plotting_data$Opacity <- 1
plotting_data$text_info <- paste(rownames(plotting_data), '\n', 'Cluster:', plotting_data$spot_groups)

x_min_image <- 0
x_max_image <- ncol(spatial_image)
y_min_image <- 0
y_max_image <- nrow(spatial_image)


# Convert spot_groups to character for consistency
spot_groups <- as.character(spot_groups)

server <- function(input, output, session) {
  # Reactive values
  v <- reactiveValues(
    feature = NULL,
    plotly_data = plotting_data,
    plotly_title = "PSF clusters",
    value_type = "psf",
    spatial_plotly = NULL,
    graphnel_df = NULL,
    pathway_image = NULL,
    pathway_name = NULL,
    node_selection_error = "",
    plotly_selection = NULL,
    selected_feature = "",
    highlight_feature_id = NULL,
    vis_samples = NULL,
    feature_metdata = NULL,
    selected_node_name = setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("node_id", "component_id_s", "label"))
  )
  
  updateSelectInput(session, "cluster", choices = c("all", as.character(unique_clusters)), selected = "all")
  updateSelectizeInput(session, "selected_pathway", choices = names(psf_saptial_results$spatial_psf_collection), server = FALSE)
  
  #### rendering spatial plotly ####
  observe({
    # Create the Plotly plot
    v$spatial_plotly <- plot_ly(
      data = v$plotly_data,
      x = ~imagecol,
      y = ~imagerow,
      customdata = rownames(v$plotly_data),
      type = 'scatter',
      mode = 'markers',
      text = ~text_info,
      hoverinfo = 'text',
      marker = list(
        size = 6,
        color = ~ColorCode,
        opacity = ~Opacity
      )
    ) %>%
      layout(
        title = list(text = v$plotly_title, x = 0.2),
        images = list(
          list(
            source = raster2uri(spatial_image),
            xref = "x",
            yref = "y",
            x = x_min_image,
            y = y_max_image,
            sizex = (x_max_image - x_min_image),
            sizey = (y_max_image - y_min_image),
            sizing = "stretch",
            opacity = 1,
            layer = "below"
          )
        ),
        xaxis = list(
          title = "",
          range = c(x_min_image, x_max_image),
          scaleanchor = "y",
          showgrid = FALSE,
          zeroline = FALSE,
          showticklabels = FALSE
        ),
        yaxis = list(
          title = "",
          range = c(y_min_image, y_max_image),
          scaleanchor = "x",
          showgrid = FALSE,
          zeroline = FALSE,
          showticklabels = FALSE
        ),
        showlegend = FALSE,  # Hide legend to prevent clutter
        dragmode = "pan"
      ) %>%
      config(scrollZoom = TRUE)
  })
  
  
  # Observe the selected cluster and update the plot accordingly
  observeEvent(input$cluster, {
    
    v$selected_feature <- ""
    v$feature_metdata <- NULL
    v$plotly_selection <- NULL
    v$highlight_feature_id <- NULL
    v$plotly_title <- paste0("PSF cluster ", input$cluster)
    
    if(input$cluster != "all") {
      updateSelectizeInput(inputId = "selected_feature", 
                           choices = c("", paste0(psf_saptial_results$psf_ident_markers[[input$cluster]]$Pathway_name, "->", psf_saptial_results$psf_ident_markers[[input$cluster]]$Sink)), 
                           selected = "", server = TRUE)
    } else {
      updateSelectizeInput(inputId = "selected_feature", 
                           choices = NULL, server = TRUE)
    }
    
    if(input$cluster == "all") {
      plotting_data$ColorCode <- cluster_colors[plotting_data$spot_groups]
      plotting_data$Opacity <- 1
    } else {
      # Add a column to indicate whether to highlight the point
      plotting_data$Highlight <- ifelse(
        plotting_data$spot_groups == input$cluster, "Highlighted", "Normal"
      )
      
      # Set color and opacity based on whether the point is highlighted
      plotting_data$ColorCode <- ifelse(
        plotting_data$Highlight == "Highlighted",
        cluster_colors[plotting_data$spot_groups],
        "grey"
      )
      
      plotting_data$Opacity <- ifelse(
        plotting_data$Highlight == "Highlighted",
        1,
        0.3
      )
    }
    
    plotting_data$text_info <- paste(rownames(plotting_data), '\n', 'Cluster:', plotting_data$spot_groups)
    
    v$plotly_data <- plotting_data
    
  })
  
  observeEvent(input$selected_feature, {
    v$selected_feature <- input$selected_feature
  })
  
  # Initial spatial plot rendering
  output$spatial_plot <- renderPlotly({
    v$spatial_plotly %>%
      event_register()
  })
  
  #### plotly selection events ####
  observeEvent(event_data("plotly_selected"), {
    v$plotly_selection <- event_data("plotly_selected")$customdata
  })
  
  observeEvent(event_data("plotly_deselect"), {
    v$plotly_selection <- event_data("plotly_selected")$customdata
  })
  
  #### pathway plotting variables observer ####
  observeEvent(c(v$selected_feature, input$selected_pathway, input$app_mode, v$plotly_selection, input$cluster), {
    
    if(!is.null(v$plotly_selection)) {
      v$vis_samples <- v$plotly_selection
    } else {
      if(input$cluster == "all") {
        v$vis_samples <- rownames(psf_saptial_results$meta.data)
      } else {
        v$vis_samples <- rownames(psf_saptial_results$meta.data)[which(as.character(psf_saptial_results$meta.data$seurat_clusters) == input$cluster)]
      }
    }
    
    if(input$app_mode == "Top features") {
      if(v$selected_feature != "") {
        v$highlight_feature_id <- unlist(strsplit(sink_name_to_id[v$selected_feature], split = "; "))[1]
        v$pathway_name <- unlist(strsplit(sink_name_to_id[v$selected_feature], split = "; "))[2]
        v$plotly_title <- paste0(v$selected_feature, " PSF")
      }
      
    } else {
      if(!is.null(input$selected_pathway)) {
        v$pathway_name <- input$selected_pathway
      }
    }
    
  })
  
  observeEvent(v$selected_feature, {
    if(input$cluster != "all" & v$selected_feature != "") {
      metdata_text <- psf_saptial_results$psf_ident_markers[[input$cluster]][gsub("_", "-", sink_name_to_id[v$selected_feature]),c("avg_log2FC", "p_val_adj", "Process", "Dowstream_pathway")]
      
      v$feature_metdata <- paste0(
        paste("<font color=\"#000000\"><b>",
              "Avg log2FC: ", round(metdata_text$avg_log2FC, digits = 3),
              "; P_adj: ", metdata_text$p_val_adj, "<br>",
              "</b></font>"), 
        paste("<font color=\"#1f992f\"><b>",
              "Process: ", metdata_text$Process,  "<br>",
              "Downstream pathway:", metdata_text$Dowstream_pathway,
              "</b></font>")
      )
    } else {
      v$feature_metdata <- NULL
    }
  })
  
  #### Selected cluster feature information #### 
  output$feature_stat <- renderText({
    v$feature_metdata
  })
  
  
  #### colored featureplot data generator ####
  observeEvent(v$highlight_feature_id, {
    feature_values <- log(psf_saptial_results$spatial_psf_collection[[v$pathway_name]][["psf_activities"]][v$highlight_feature_id,])
    
    plotting_data <- data.frame(coords, spot_groups = spot_groups)
    plotting_data$ColorCode <- map_colors(values = feature_values, colors = Seurat:::SpatialColors(n = 100)) # Seurat:::SpatialColors(n = 100)[cut(feature_values, breaks = 100, labels = FALSE)]
    plotting_data$Opacity <- 1
    plotting_data$text_info <- paste(rownames(plotting_data), '\n', 'Cluster:', plotting_data$spot_groups, '\n', "Log PSF", round(feature_values, digits = 5))
    v$plotly_data <- plotting_data
  })
  
  #### pathway coloring and plotting observer ####
  observeEvent(c(v$pathway_name, v$vis_samples, v$highlight_feature_id, v$value_type, input$node_coloring_type), {
    if(!is.null(v$pathway_name)) {
      if(v$pathway_name != "") {
        clust_pathway <- psf_subset(psf_collection = psf_saptial_results$spatial_psf_collection, 
                                    pathway_name = v$pathway_name, 
                                    sample_list = v$vis_samples)
        
        if(!is.null(clust_pathway)) {
          v$node_colors <- node_color_generator(pathway = clust_pathway, log_norm = T, sample_id = "mean", 
                                                color_nodes = c("psf_activities", "fold_change")[as.integer(input$node_coloring_type)], 
                                                custom_color_scale = Seurat:::SpatialColors(n = 100))
        }
        
        v$graphnel_df <- graphnel_to_df(clust_pathway, extended = TRUE)
        
        v$graphnel_df$node_table[rownames(v$node_colors$node_values), "expression"] <- v$node_colors$node_values[,1]
        v$graphnel_df$node_table[rownames(v$node_colors$node_values), "signal"] <- v$node_colors$node_values[,2]
        
        
        img_path <- system.file("extdata", "old_imgs", paste0(gsub("path:", "", clust_pathway$attrs$name), ".png"), package = "psf")
        
        v$pathway_image <- magick::image_read(img_path)
        
        v$image_file <- plot_kegg_pathway(graphnel_df = v$graphnel_df, group_graphics = NULL, pathway_image = v$pathway_image,
                                          node_colors = v$node_colors$node_colors, edge_mapping = FALSE, edge_in_mode = FALSE, 
                                          highlight_nodes = v$highlight_feature_id, highlight_color = "green", 
                                          col_legend_title = v$node_colors$col_legend_title, color_bar_lims = v$node_colors$color_bar_lims, 
                                          present_node_modifications = FALSE, removed_nodes = NULL, custom_color_scale = Seurat:::SpatialColors(n = 100)) %>% 
          image_write(tempfile(fileext='png'), format = 'png')
      }
    }
  })
  
  
  #### image click events ####
  observeEvent(input$image_click, {
    
    if(!is.null(input$multiple_choice)) {
      if(input$multiple_choice == "yes") {
        v$selected_node_name <- unique(rbind(v$selected_node_name, v$graphnel_df$node_table[which(input$image_click[[1]] <= v$graphnel_df$node_table$x_end & input$image_click[[1]] >= v$graphnel_df$node_table$x_start &
                                                                                                    input$image_click[[2]] <= v$graphnel_df$node_table$y_end & input$image_click[[2]] >= v$graphnel_df$node_table$y_start), c("node_id", "component_id_s", "label")]))
      } else {
        v$selected_node_name <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("node_id", "component_id_s", "label"))
        v$selected_node_name <- v$graphnel_df$node_table[which(input$image_click[[1]] <= v$graphnel_df$node_table$x_end & input$image_click[[1]] >= v$graphnel_df$node_table$x_start &
                                                                 input$image_click[[2]] <= v$graphnel_df$node_table$y_end & input$image_click[[2]] >= v$graphnel_df$node_table$y_start), c("node_id", "component_id_s", "label")]
      }
    } else {
      v$selected_node_name <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("node_id", "component_id_s", "label"))
      v$selected_node_name <- v$graphnel_df$node_table[which(input$image_click[[1]] <= v$graphnel_df$node_table$x_end & input$image_click[[1]] >= v$graphnel_df$node_table$x_start &
                                                               input$image_click[[2]] <= v$graphnel_df$node_table$y_end & input$image_click[[2]] >= v$graphnel_df$node_table$y_start), c("node_id", "component_id_s", "label")]
    }
    
    
    v$image_file <- plot_kegg_pathway(graphnel_df = v$graphnel_df, group_graphics = NULL, pathway_image = v$pathway_image,
                                      node_colors = v$node_colors$node_colors, edge_mapping = FALSE, edge_in_mode = FALSE, highlight_nodes = v$selected_node_name$node_id, highlight_color = "red",
                                      col_legend_title = v$node_colors$col_legend_title, color_bar_lims = v$node_colors$color_bar_lims, present_node_modifications = FALSE, removed_nodes = NULL, custom_color_scale = Seurat:::SpatialColors(n = 100)) %>%
      image_write(tempfile(fileext='png'), format = 'png')
    
  })
  
  #### image brush events ####
  observeEvent(input$image_brush, {
    if(!is.null(input$multiple_choice)) {
      if(input$multiple_choice == "yes") {
        
        v$selected_node_name <- unique(rbind(v$selected_node_name, v$graphnel_df$node_table[which(input$image_brush[1:4]$xmax >= v$graphnel_df$node_table$x_end & input$image_brush[1:4]$xmin <= v$graphnel_df$node_table$x_start &
                                                                                                    input$image_brush[1:4]$ymax >= v$graphnel_df$node_table$y_end & input$image_brush[1:4]$ymin <= v$graphnel_df$node_table$y_start), c("node_id", "component_id_s", "label")]))
        
      } else {
        v$selected_node_name <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("node_id", "component_id_s", "label"))
        v$selected_node_name <- v$graphnel_df$node_table[which(input$image_brush[1:4]$xmax >= v$graphnel_df$node_table$x_end & input$image_brush[1:4]$xmin <= v$graphnel_df$node_table$x_start &
                                                                 input$image_brush[1:4]$ymax >= v$graphnel_df$node_table$y_end & input$image_brush[1:4]$ymin <= v$graphnel_df$node_table$y_start), c("node_id", "component_id_s", "label")]
      }
    } else {
      v$selected_node_name <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("node_id", "component_id_s", "label"))
      v$selected_node_name <- v$graphnel_df$node_table[which(input$image_brush[1:4]$xmax >= v$graphnel_df$node_table$x_end & input$image_brush[1:4]$xmin <= v$graphnel_df$node_table$x_start &
                                                               input$image_brush[1:4]$ymax >= v$graphnel_df$node_table$y_end & input$image_brush[1:4]$ymin <= v$graphnel_df$node_table$y_start), c("node_id", "component_id_s", "label")]
    }
    
    v$image_file <- plot_kegg_pathway(graphnel_df = v$graphnel_df, group_graphics = NULL, pathway_image = v$pathway_image,
                                      node_colors = v$node_colors$node_colors, edge_mapping = FALSE, edge_in_mode = FALSE, highlight_nodes = v$selected_node_name$node_id, highlight_color = "red",
                                      col_legend_title = v$node_colors$col_legend_title, color_bar_lims = v$node_colors$color_bar_lims, present_node_modifications = FALSE, removed_nodes = NULL, custom_color_scale = Seurat:::SpatialColors(n = 100)) %>%
      image_write(tempfile(fileext='png'), format = 'png')
    
    
  })
  
  observeEvent(input$clicked_node$nodes, {
    v$selected_node_name <- v$graphnel_df$node_table[which(v$graphnel_df$node_table$node_id == unlist(input$clicked_node$nodes)),c("node_id", "component_id_s", "label")]
  })
  
  
  #### plotly visualization by node selection ####
  observeEvent(c(v$selected_node_name, input$node_coloring_type), {
    if(nrow(v$selected_node_name) == 1) {
      
      plotting_data <- data.frame(coords, spot_groups = spot_groups)
      
      feature_values <- log(psf_saptial_results$spatial_psf_collection[[v$pathway_name]][[c("psf_activities", "exp_fc")[as.integer(input$node_coloring_type)]]][v$selected_node_name$node_id,])
      
      plotting_data$ColorCode <-  map_colors(values = feature_values, colors = Seurat:::SpatialColors(n = 100)) # Seurat:::SpatialColors(n = 100)[cut(feature_values, breaks = 100, labels = FALSE)]
      plotting_data$Opacity <- 1
      plotting_data$text_info <- paste(rownames(plotting_data), '\n', 'Cluster:', plotting_data$spot_groups, '\n', "Log PSF", round(feature_values, digits = 5))
      
      v$plotly_data <- plotting_data
      
      v$plotly_title <- paste0(v$selected_node_name$label, c(" Log PSF", " Log Exp FC")[as.integer(input$node_coloring_type)])
    } else {
      if(nrow(v$selected_node_name) > 1) {
        v$node_selection_error <- "Only one node can be visualized"
        delay(3000, v$node_selection_error <- "")
      }
    }
  })
  
  output$error_text <- renderText({
    paste("<font color=\"#ff0000\"><b>", v$node_selection_error, "</b></font>")
  })
  
  #### image hover events ####
  output$hover_text <- renderUI({
    h5(textOutput("hover_data"))
  })
  
  output$hover_data <- renderText({
    if(length(which(input$image_hover[[1]] <= v$graphnel_df$node_table$x_end & input$image_hover[[1]] >= v$graphnel_df$node_table$x_start &
                    input$image_hover[[2]] <= v$graphnel_df$node_table$y_end & input$image_hover[[2]] >= v$graphnel_df$node_table$y_start)) > 0 ) {
      
      gr_name <- v$graphnel_df$node_table[which(input$image_hover[[1]] <= v$graphnel_df$node_table$x_end & input$image_hover[[1]] >= v$graphnel_df$node_table$x_start &
                                                  input$image_hover[[2]] <= v$graphnel_df$node_table$y_end & input$image_hover[[2]] >= v$graphnel_df$node_table$y_start),"title"]
      
      if(!is.null(v$node_colors)) {
        node_id <- v$graphnel_df$node_table[which(input$image_hover[[1]] <= v$graphnel_df$node_table$x_end & input$image_hover[[1]] >= v$graphnel_df$node_table$x_start &
                                                    input$image_hover[[2]] <= v$graphnel_df$node_table$y_end & input$image_hover[[2]] >= v$graphnel_df$node_table$y_start),"node_id"]
        gr_name <- paste0(gr_name, "\n", v$node_colors$hover_text[node_id])
      }
      
      ### fix next line issue
      gsub("<br>", "\n", gr_name)
      
    }
  })
  
  
  output$pathway_image <- renderImage({
    return(list(src = ifelse(is.null(v$graphnel_df), "", v$image_file), contentType = "image/png", alt = "pathway_image"))
  }, deleteFile = TRUE)
  
  
  #### visnet rendering ####
  output$visnet <- renderVisNetwork({
    
    if(!is.null(v$graphnel_df)) {
      visNetwork(nodes = vis_extract(v$graphnel_df, node_colors = v$node_colors, higlight_node = v$highlight_feature_id, custom_color_scale = Seurat:::SpatialColors(n = 100))$node_table, 
                 edges = vis_extract(v$graphnel_df)$edge_table, width = "100%", height = "800px") %>% 
        visIgraphLayout(layout = "layout_nicely") %>% 
        visInteraction(navigationButtons = TRUE, multiselect = F, keyboard = F, selectConnectedEdges = T) %>%
        visOptions(highlightNearest = T) %>%
        visExport() %>%
        visEvents(click = "function(nodes) {
                      console.info('click')
                      console.info(nodes)
                      Shiny.onInputChange('clicked_node', {nodes : nodes.nodes, edges : nodes.edges});
                      ;}"
        )
    }
    
  })
  
}
