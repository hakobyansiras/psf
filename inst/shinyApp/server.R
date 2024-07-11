library(magick)
library(shiny)
library(shape)
library(psf)
library(DT)
library(plotly)
library(data.table)
library(shinyjs)
library(visNetwork)
library(shinyjqui)
library(ggplot2)
library(igraph)
library(shinyhelper)
### library(plotrix)

load("whole_data_unit.RData")
protein_drug_interactions <- fread("Protein_drug_interactions.tsv")
options(shiny.maxRequestSize=50*1024^2)

#############################################
######### Pathway curation block ############
#############################################

#### pathway kgml and image downloader ####
kegg_data_downloader <- function(pathway_name, collection_name) {
  
  pathway_id <- pathway_codes_new[pathway_name]
  
  kgml <- paste0("http://rest.kegg.jp/get/", pathway_id, "/kgml")
  
  image <- paste0("http://rest.kegg.jp/get/", pathway_id, "/image")
  
  kgml_path <- tempfile()
  
  download.file(kgml, destfile = kgml_path, method = "auto")
  
  # image_path <- tempfile()
  image_path <- paste0("collection_dir/", collection_name, "/pathway_images/", pathway_id, ".png")
  
  download.file(image, destfile = image_path, method = "auto")
  
  img <- magick::image_read(image_path)
  
  graph <- psf::generate.kegg.collection.from.kgml(kgml_path)[[1]]
  
  return(list(pathway = graph, image = img))
  
}

#### kegg pathway plotter ####
plot_kegg_pathway <- function(graphnel_df, group_graphics, pathway_image, 
                              node_colors = NULL, custom_edge_mapping = FALSE, edge_mapping = FALSE, edge_in_mode = FALSE, highlight_nodes = NULL, highlight_color = "red",
                              col_legend_title = NULL, color_bar_lims, y_adj_sink = 3, node_fill_opacity = 0, present_node_modifications = FALSE, removed_nodes = NULL, y_adj_text = 3
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

#### function for graphical data generation ####
graphical_data_generator <- function(pathway, include_changes = FALSE) {
  
  entrez_id <- unname(sapply(pathway$graph@nodes, function(y) {
    ifelse(is.null(unlist(graph::nodeData(pathway$graph, y, attr = "genes"))),
           as.character(unlist(graph::nodeData(pathway$graph, y, attr = "label"))), 
           paste(as.character(unlist(graph::nodeData(pathway$graph, y, attr = "genes"))), collapse = ","))
  }))
  
  hover_name <- unname(sapply(pathway$graph@nodes, function(y) {
    ifelse(is.null(unlist(graph::nodeData(pathway$graph, y, attr = "genes"))),
           paste0(as.character(unlist(graph::nodeData(pathway$graph, y, attr = "label"))), " (",
                  unlist(strsplit(kegg_compounds_to_full_name[as.character(unlist(graph::nodeData(pathway$graph, y, attr = "label"))),], split = ";"))[1],
                  ")"
           ),
           paste0(
             paste0(as.character(unlist(graph::nodeData(pathway$graph, y, attr = "genes"))), 
                    paste0("(", entrez_to_symbol[as.character(unlist(graph::nodeData(pathway$graph, y, attr = "genes"))),], ")")),
             collapse = ", "
           )
    )
  }))
  
  kegg_coords_new <- data.frame(
    pathway = rep(pathway$attrs$name, length(graph::nodeData(pathway$graph, attr = "kegg.gr.x"))),
    x_start = as.integer(unlist(graph::nodeData(pathway$graph, attr = "kegg.gr.x"))) - as.integer(unlist(graph::nodeData(pathway$graph, attr = "kegg.gr.width")))*0.5, 
    y_start = as.integer(unlist(graph::nodeData(pathway$graph, attr = "kegg.gr.y"))) - as.integer(unlist(graph::nodeData(pathway$graph, attr = "kegg.gr.height")))*0.5,
    x_end  = as.integer(unlist(graph::nodeData(pathway$graph, attr = "kegg.gr.x"))) + as.integer(unlist(graph::nodeData(pathway$graph, attr = "kegg.gr.width")))*0.5,
    y_end =  as.integer(unlist(graph::nodeData(pathway$graph, attr = "kegg.gr.y"))) + as.integer(unlist(graph::nodeData(pathway$graph, attr = "kegg.gr.height")))*0.5,
    node_name = as.character(unlist(graph::nodeData(pathway$graph, attr = "kegg.name"))), 
    gr_name = as.character(unlist(graph::nodeData(pathway$graph, attr = "label"))),
    node_id = as.character(unlist(graph::nodeData(pathway$graph, attr = "kegg.id"))),
    entrez_id = entrez_id, hover_name = hover_name,
    sink = (unlist(graph::nodeData(pathway$graph, attr = "kegg.id")) %in% pathway$sink.nodes), 
    exist = (unlist(graph::nodeData(pathway$graph, attr = "existence")) == "exist"),
    changed = (unlist(graph::nodeData(pathway$graph, attr = "data_source")) != "kegg"),
    node_class = as.character(unlist(graph::nodeData(pathway$graph, attr = "type"))), stringsAsFactors = F 
  )
  
  kegg_arrows_type <- c("simple", "simple", "simple", "simple", "T", "simple", "simple", "simple", "T", "simple", "", "simple", "simple", "T", "simple", "simple")
  names(kegg_arrows_type) <- c("activation", "binding/association", "compound", "dephosphorylation", "dissociation", "expression", "glycosylation", "indirect effect", "inhibition", "missing interaction", "n/a", "phosphorylation", "reaction", "repression", "state change", "ubiquitination" )
  
  line_col <- c("red", "red", "red", "red", "blue","red", "red", "red", "blue", "red", "red", "red", "red", "blue", "red", "red")
  names(line_col) <- c("activation", "binding/association", "compound", "dephosphorylation", "dissociation", "expression", "glycosylation", "indirect effect", "inhibition", "missing interaction", "n/a", "phosphorylation", "reaction", "repression", "state change", "ubiquitination" )
  
  # from <- as.character(graph::edgeMatrix(pathway$graph)[1,])
  # to <- as.character(graph::edgeMatrix(pathway$graph)[2,])
  # 
  # if(sum(is.na(c(match(from, kegg_coords_new$node_id), match(to, kegg_coords_new$node_id)))) > 0) {
  #   splitted_interactions <- strsplit(graph::edgeNames(pathway$graph), split = "~")
  #   from <- as.character(sapply(splitted_interactions, "[[", 1))
  #   to <- as.character(sapply(splitted_interactions, "[[", 2))
  # }
  
  if(length(names(pathway$graph@edgeData@data)) > 0) {
    splitted_interactions <- strsplit(names(pathway$graph@edgeData@data), split = "|", fixed = T)
    from <- as.character(sapply(splitted_interactions, "[[", 1))
    to <- as.character(sapply(splitted_interactions, "[[", 2))
  } else {
    return(list(node_coords = kegg_coords_new, edge_coords = NULL))
  }
  
  from_index <- match(from, kegg_coords_new$node_id)
  to_index <- match(to, kegg_coords_new$node_id)
  
  col <- line_col[unname(unlist(graph::edgeData(pathway$graph, from = from, to = to, attr = "subtype1")))]
  col[which(unname(unlist(graph::edgeData(pathway$graph, from = from, to = to, attr = "data_source"))) != "kegg" & col == "red")] <- "purple"
  col[which(unname(unlist(graph::edgeData(pathway$graph, from = from, to = to, attr = "data_source"))) != "kegg" & col == "blue")] <- "yellow"
  
  lty_type <- c("solid", "solid", "dotted")
  names(lty_type) <- c("exist", "added", "removed")
  
  arrow_data <- data.frame(
    from = from, to = to,
    x0 = (kegg_coords_new$x_start[from_index] + kegg_coords_new$x_end[from_index])/2,
    x1 = (kegg_coords_new$x_start[to_index] + kegg_coords_new$x_end[to_index])/2,
    y0 = (kegg_coords_new$y_start[from_index] + kegg_coords_new$y_end[from_index])/2,
    y1 = (kegg_coords_new$y_start[to_index] + kegg_coords_new$y_end[to_index])/2,
    arr.type = kegg_arrows_type[unname(unlist(graph::edgeData(pathway$graph, from = from, to = to,attr = "subtype1")))],
    col = col, lty = lty_type[unname(unlist(graph::edgeData(pathway$graph, from = from, to = to,attr = "existence")))],
    stringsAsFactors = F
  )
  
  if(!include_changes) {
    arrow_data <- arrow_data[which(arrow_data$lty != "dotted"),]
  }
  
  return(list(node_coords = kegg_coords_new, edge_coords = arrow_data))
}

#### kegg mapper function ####
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

#### database search function ####
node_pairs_generator <- function(selected_node_name, direction = TRUE, direction_state = TRUE, pathway_realted = FALSE, pathway_data) {
  
  if(pathway_realted) {
    node_a <- unlist(strsplit(selected_node_name[1,"component_id_s"], split = ","))
    
    pathway_nodes <- do.call(rbind, apply(pathway_data, 1, function(x) {
      cbind(unlist(strsplit(x[3], split = ",")), 
            rep(x[1], length(unlist(strsplit(x[3], split = ",")))))
    }))
    
    node_names <- setNames(object = pathway_nodes[,2], nm = pathway_nodes[,1])
    
    if(direction) {
      
      if(direction_state) {
        pathway_filtered_database <- all_interactions[entrez_a %in% pathway_nodes[,1] & entrez_b %in% node_a]
      } else {
        pathway_filtered_database <- all_interactions[entrez_a %in% node_a & entrez_b %in% pathway_nodes[,1]]
      }
      
    } else {
      
      pathway_filtered_database <- all_interactions[entrez_a %in% pathway_nodes[,1] & entrez_b %in% node_a | entrez_a %in% node_a & entrez_b %in% pathway_nodes[,1]]
      
    }
    
    pathway_filtered_database <- cbind(pathway_filtered_database, node_names[pathway_filtered_database$entrez_a], node_names[pathway_filtered_database$entrez_b])
    
    colnames(pathway_filtered_database)[9:10] <- c("node_a_id", "node_b_id")
    
    return(pathway_filtered_database)
    
  } else {
    node_a <- unlist(strsplit(selected_node_name[1,"component_id_s"], split = ","))
    node_b <- unlist(strsplit(selected_node_name[2,"component_id_s"], split = ","))
    
    if(direction) {
      if(direction_state) {
        filtered_interactons <- all_interactions[entrez_a %in% node_a & entrez_b %in% node_b]
        node_a_id <- rep(selected_node_name[1,"node_id"], nrow(filtered_interactons))
        node_b_id <- rep(selected_node_name[2,"node_id"], nrow(filtered_interactons))
        return(cbind(filtered_interactons, node_a_id, node_b_id))
      } else {
        filtered_interactons <- all_interactions[entrez_a %in% node_b & entrez_b %in% node_a]
        node_a_id <- rep(selected_node_name[1,"node_id"], nrow(filtered_interactons))
        node_b_id <- rep(selected_node_name[2,"node_id"], nrow(filtered_interactons))
        return(cbind(filtered_interactons, node_b_id, node_a_id))
      }
    } else {
      node_a_id <- selected_node_name[1,"node_id"]
      node_b_id <- selected_node_name[2,"node_id"]
      filtered_interactons_ab <- all_interactions[entrez_a %in% node_a & entrez_b %in% node_b]
      filtered_interactons_ba <- all_interactions[entrez_a %in% node_b & entrez_b %in% node_a]
      return(
        rbind(
          cbind(filtered_interactons_ab, rep(node_a_id, nrow(filtered_interactons_ab)), rep(node_b_id, nrow(filtered_interactons_ab))),
          cbind(filtered_interactons_ba, rep(node_b_id, nrow(filtered_interactons_ba)), rep(node_a_id, nrow(filtered_interactons_ba)))
        )
      )
    }
  }
  
}

#### visnet data extractor ####
vis_extract <- function(vis_table, node_colors = NULL) {
  node_table <- vis_table$node_table[,c("node_id", "label", "vis_width", "font.color", "size", "font.size", "borderWidth", "color.border", "color.background", "shape", "title", "x", "y")]
  
  colnames(node_table) <- c("id", "label", "widthConstraint.maximum", "font.color", "size", "font.size", "borderWidth", "color.border", "color.background", "shape", "title", "x", "y")
  
  if(!is.null(node_colors)) {
    
    magick::autoviewer_disable() ### to avoid legend plotting before netowrk rendering
    legend_img <- magick::image_device(width = 480, height = 480)
    plot.new()
    
    color_legend_maker(x = 0.05, y = 0, leg = 0.9, cols = c(pal1(10), pal2(10)), title = node_colors$col_legend_title, lims = node_colors$color_bar_lims, digits=3, prompt=FALSE,
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

visnet_creator <- function(graphical_data, node_colors = NULL, col_legend_title = "", color_bar_lims = "") {
  graphical_data$edge_coords <- graphical_data$edge_coords[which(graphical_data$edge_coords$lty == "solid"),]
  graphical_data$node_coords <- graphical_data$node_coords[which(graphical_data$node_coords$exist),]
  
  node_shapes <- unname(sapply(graphical_data$node_coords$node_name, function(x) {
    if(grepl("cpd",x)) {
      "dot"
    } else {
      "box"
    }
  }))
  
  node_widths <- unname(sapply(1:nrow(graphical_data$node_coords), function(x) {
    if(grepl("path", graphical_data$node_coords[x, "node_name"])) {
      (graphical_data$node_coords[x,"x_end"] - graphical_data$node_coords[x,"x_start"])*2
    } else {
      NA
    }
  }))
  
  if(!is.null(node_colors)) {
    color <- unname(sapply(graphical_data$node_coords$node_id, function(x) {
      if(x %in% node_colors$node_id) {
        node_colors[which(node_colors$node_id == x),"col"]
      } else {
        "#BFFFBF"
      }
    }))
    
    font_color <- unname(sapply(graphical_data$node_coords$node_id, function(x) {
      if(x %in% node_colors$node_id) {
        node_colors[which(node_colors$node_id == x),"text_col"]
      } else {
        "#000000"
      }
    }))
    
    border_color <- unname(sapply(graphical_data$node_coords$sink, function(x) {
      if(x) {
        "#0099cc"
      } else {
        "#BFFFBF"
      }
    }))
    
  } else {
    color <- unname(sapply(graphical_data$node_coords$sink, function(x) {
      if(x) {
        "#0099cc"
      } else {
        "#BFFFBF"
      }
    }))
    
    border_color <- rep("#BFFFBF", nrow(graphical_data$node_coords))
    
    font_color <- rep("#000000", nrow(graphical_data$node_coords))
    
  }
  
  size <- unname(sapply(graphical_data$node_coords$node_name, function(x) {
    if(grepl("cpd",x)) {
      10
    } else {
      25
    }
  }))
  
  nodes <- data.frame(id = graphical_data$node_coords$node_id,
                      image = rep("unselected", nrow(graphical_data$node_coords)),
                      label = graphical_data$node_coords$gr_name,
                      shape = node_shapes, 
                      color = color, 
                      color.border = border_color, 
                      title = graphical_data$node_coords$hover_name,
                      borderWidth = 2,
                      font.size = rep(22, nrow(graphical_data$node_coords)), 
                      size = size,
                      font.color = font_color,
                      widthConstraint.maximum = node_widths,
                      x = (graphical_data$node_coords$x_start + graphical_data$node_coords$x_end)/2,
                      y = (graphical_data$node_coords$y_start + graphical_data$node_coords$y_end)/2,
                      row.names = graphical_data$node_coords$node_id, stringsAsFactors = F
  )
  
  if(!is.null(node_colors)) {
    
    legend_img <- magick::image_device(width = 480, height = 480)
    plot.new()
    
    color_legend_maker(x = 0.05, y = 0, leg = 0.9, cols = c(pal1(10), pal2(10)), title = col_legend_title, lims = color_bar_lims, digits=3, prompt=FALSE,
                       lwd=4, outline=TRUE, subtitle = "", fsize = 1.3)
    
    temp_legend <- tempfile()
    
    magick::image_write(magick::image_trim(legend_img, fuzz = 0), path = temp_legend)
    
    legend_path <- paste('data:image/png;base64', RCurl::base64Encode(readBin(temp_legend, 'raw', file.info(temp_legend)[1, 'size']), 'txt'), sep = ',')
    
    legend_data_frame <- data.frame(
      id = as.character(max(as.integer(graphical_data$node_coords$node_id)) + 1),
      image = legend_path,
      label = "Color legend", 
      shape = "image",
      color = "", 
      color.border = "", 
      title = "",
      borderWidth = 0,
      font.size = 32, size = 30, font.color = "#000000",
      x = max(graphical_data$node_coords$x_end),
      y = min(graphical_data$node_coords$y_end) - 10
    )
    
    nodes <- rbind(nodes, legend_data_frame)
    
  }
  
  arrows_type <- c("arrow","bar")
  names(arrows_type) <- c("simple", "T")
  
  edges <- data.frame(from = graphical_data$edge_coords$from, to = graphical_data$edge_coords$to,
                      color = graphical_data$edge_coords$col,
                      arrows.to.enabled = rep(TRUE, length(graphical_data$edge_coords$col)),
                      arrows.to.type = arrows_type[graphical_data$edge_coords$arr.type]
  )
  
  return(list(nodes = nodes, edges = edges))
  
}

database_table_design <- function(interaction_table) {
  if(nrow(interaction_table) > 0) {
    interaction_table <- as.data.frame(t(apply(interaction_table, 1, function(x) {
      if(x[5] != "") {
        ids <- unlist(strsplit(x[5], split = ";"))[1:2]
        x[5] <- paste(paste0("<a href=\"", "https://www.ncbi.nlm.nih.gov/pubmed/", ids, "\" target=\"_blank\">", ids, "</a>"), collapse = ",")
      }
      x
    })))
  }
  return(interaction_table)
}

network_table_design <- function(table) {
  if("from" %in% colnames(table)) {
    return(table[,c("from", "to", "type", "subtype1", "subtype2", "weight", "existence", "change_info", "data_source")])
  } else {
    return(table[,c("node_id", "label", "component_id_s", "type", "psf_function", "expression", "signal", "existence", "change_info", "data_source")])
  }
}

#############################################
######## Gene exp mapping and PSF ###########
#############################################

#### psf calculation and coloring function ####
psf_signal_calculator_and_coloring <- function(entrez_fc, pathway, pathway_name, update_mod = FALSE, old_graph = NULL, no_color_mode = F) {
  
  if(update_mod) {
    psf_graph <- psf.from.env.entrez.fc(entrez.fc = entrez_fc,
                                 kegg.collection = old_graph[[1]][pathway_name], calculate.significance = F, sum = FALSE, map_exp_data = FALSE)
    
  } else {
    psf_graph <- psf.from.env.entrez.fc(entrez.fc = entrez_fc,
                                        kegg.collection = setNames(object = list(pathway), nm = pathway_name), 
                                        calculate.significance = F, sum = FALSE)
  }
  
  exp_values <- unlist(graph::nodeData(psf_graph[[1]][[pathway_name]]$graph, attr = "expression"))[which(unlist(graph::nodeData(psf_graph[[1]][[pathway_name]]$graph, attr = "type")) == "gene")]
  
  exp_values <- log(exp_values[order(exp_values)] + 0.00001)
  
  if(no_color_mode) {
    exp_colors <- NULL
  } else {
    exp_colors <- color_code(values = exp_values, pal1 = pal1, pal2 = pal2)
    
    # exp_colors <- c(pal1(10)[cut(exp_values[which(exp_values <= 0)],10)], 
    #                 pal2(10)[cut(exp_values[which(exp_values > 0)],10)])
    
    # exp_colors <- c(pal1(length(unname(exp_values[which(exp_values <= 1)])))[rank(exp_values[which(exp_values <= 1)])],
    #                 pal2(length(unname(exp_values[which(exp_values > 1)])))[rank(exp_values[which(exp_values > 1)])]
    # )
    
    # exp_colors <- pal(length(unname(exp_values)))[rank(exp_values)]
    
    exp_colors = data.frame(node_id = names(exp_values)[c(which(exp_values <= 0), which(exp_values > 0))], 
                            col = exp_colors,
                            text_col = unname(sapply(exp_colors, function(x) {c( "black", "white")[  1+(sum( col2rgb(x) *c(299, 587,114))/1000 < 123) ]})),
                            stringsAsFactors = F
    )
  }
  
  signal_values <- unlist(graph::nodeData(psf_graph[[1]][[pathway_name]]$graph, attr = "signal"))[which(unlist(graph::nodeData(psf_graph[[1]][[pathway_name]]$graph, attr = "type")) == "gene")]
  
  signal_values <- log(signal_values[order(signal_values)] + 0.00001)
  
  psf_colors <- color_code(values = signal_values, pal1 = pal1, pal2 = pal2)
  
  # psf_colors <- c(pal1(10)[cut(signal_values[which(signal_values <= 0)],10)], 
  #                 pal2(10)[cut(signal_values[which(signal_values > 0)],10)])
  
  # psf_colors <- c(pal1(length(unname(signal_values[which(signal_values <= 1)])))[rank(signal_values[which(signal_values <= 1)])],
  #                 pal2(length(unname(signal_values[which(signal_values > 1)])))[rank(signal_values[which(signal_values > 1)])]
  # )
  
  # psf_colors <- pal(length(unname(signal_values)))[rank(signal_values)]
  
  psf_colors = data.frame(node_id = names(signal_values)[c(which(signal_values <= 0), which(signal_values > 0))], 
                          col = psf_colors,
                          text_col = unname(sapply(psf_colors, function(x) {c( "black", "white")[  1+(sum( col2rgb(x) *c(299, 587,114))/1000 < 123) ]})),
                          stringsAsFactors = F
  )
  
  sink_names <- unlist(graph::nodeData(pathway$graph, pathway$sink.nodes, attr = "label"))
  
  sink_order_by_coord <- order(as.integer(unname(unlist(graph::nodeData(pathway$graph, pathway$sink.nodes, attr = "kegg.gr.y")))), decreasing = T)
  
  sink_names <- setNames(nm = names(sink_names), object = paste0(sink_names, "_", 1:length(sink_names)))
  
  signal_data <- Reduce(rbind, 
                        lapply(psf_graph, function(x) {
                          
                          sink_signals <- round(log(x[[pathway_name]]$signal.at.sink[order(x[[pathway_name]]$signal.at.sink)]), digits = 5)
                          
                          sink_colors <- color_code(values = sink_signals, pal1 = pal1, pal2 = pal2)  
                          
                          sink_signals <- sink_signals[c(which(sink_signals <= 0), which(sink_signals > 0))]
                          
                          data.frame(sink_name = sink_names[names(sink_signals)], 
                                     signal = sink_signals, 
                                     dot_color = sink_colors,
                                     stringsAsFactors = F)
                          
                        })
  )
  
  signal_data$sink_name <- factor(signal_data$sink_name, levels = sink_names[sink_order_by_coord])
  
  return(list(psf_colors = psf_colors, exp_colors = exp_colors, 
              exp_values = exp_values, signal_values = signal_values, 
              sink_signals = signal_data, psf_graph = psf_graph))
  
}

#### node color coding function ####
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
  
  return(list(node_colors = node_colors, col_legend_title = col_legend_title, color_bar_lims = range(pathway_node_values), hover_text = hover_text, node_values = node_values))
}

# #### check user data ####
# if(file.exists("collection_dir") & length(dir("collection_dir/pathways/") > 0)) {
#   pathway <- readRDS(dir("collection_dir/pathways", full.names = T)[1])
#   pathway_name <- gsub(".RDS", "", dir("collection_dir/pathways")[1])
#   graphnel_df <- graphnel_to_df(pathway, extended = TRUE)
#   
#   if(file.exists(paste0("collection_dir/pathway_images/", pathway_codes_new[pathway_name], ".png"))) {
#     image <- magick::image_read(paste0("collection_dir/pathway_images/", pathway_codes_new[pathway_name], ".png"))
#     
#     rendering_image <- plot_kegg_pathway(graphnel_df = graphnel_df, group_graphics = pathway$group_nodes, pathway_image = image,
#                                          node_colors = NULL, edge_mapping = FALSE, edge_in_mode = FALSE, highlight_nodes = NULL, highlight_color = "red",
#                                          col_legend_title = NULL, color_bar_lims, present_node_modifications = FALSE, removed_nodes = pathway$removed_nodes) %>%
#       image_write(tempfile(fileext='png'), format = 'png')
#   } else {
#     image <- NULL
#     rendering_image <- NULL
#   }
# } else {
#   
#   if(!file.exists("collection_dir")) {
#     dir.create("collection_dir")
#     dir.create("collection_dir/pathways")
#     dir.create("collection_dir/pathway_images")
#   }
#   
#   pathway_name <- "Chemokine_signaling_pathway"
#   pathway_data <- kegg_data_downloader(pathway_name)
#   pathway <- pathway_data$pathway
#   image <- pathway_data$image
#   
#   graphnel_df <- graphnel_to_df(pathway, extended = TRUE)
#   
#   rendering_image <- plot_kegg_pathway(graphnel_df = graphnel_df, group_graphics = pathway$group_nodes, pathway_image = image,
#                                        node_colors = NULL, edge_mapping = FALSE, edge_in_mode = FALSE, highlight_nodes = NULL, highlight_color = "red",
#                                        col_legend_title = NULL, color_bar_lims, present_node_modifications = FALSE, removed_nodes = pathway$removed_nodes) %>%
#     image_write(tempfile(fileext='png'), format = 'png')
# }


shinyServer(function(input, output, session) {
  
  #### reactive values ####
  v <- reactiveValues(image_file = NULL,
                      pathway_name = NULL,
                      collection_name = NULL,
                      pathway = NULL, pathway_image = NULL,
                      graphnel_df = NULL,
                      undo_stack = list(),
                      redo_stack = NULL,
                      pathway_load_alert = 0,
                      selected_node_name = setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("node_id", "component_id_s", "label")),
                      from_nodes = NULL, database_out = NULL, searching_state = FALSE, 
                      database_search_node = NULL, selected_interaction = NULL, collection_load_error = "", edge_adding_error = "", edge_adding_error_visnet = "", node_data_adding_error_visnet = "",
                      event_node_type = scan(file = "event_names.txt", what = "character", sep = "\t"), pathway_change_alert = 0, undo_alert = 0,
                      edge_edit_stage = FALSE, previous_row_selection = 0, allow_edge_draw = TRUE,
                      psf_and_colors = NULL, draw_color_bar = FALSE, col_legend_title = "", color_bar_lims = NULL, color_bar_psf_mode = FALSE, 
                      node_colors = NULL, # pathway_exp_colored = FALSE, pathway_psf_colored = FALSE, 
                      mapping_value_type = NULL, exp_uploaded = NULL, entrez_fc = NULL,
                      allow_graph_param_update = FALSE, sink_values_plot = NULL,
                      edge_network_table_filtered = NULL, node_network_table_filtered = NULL, 
                      visnet_coord_ranges = NULL, node_edit_mode = FALSE, loaded_collection = NULL, drug_table = NULL, visnet_highlighted_nodes = NULL                   
  )
  
  observeEvent(input$add_load_collection, {
    
    if(length(dir("collection_dir/")) > 0) {
      show('collection_selector')
    } else {
      hide('collection_selector')
    }
    
    showModal(modalDialog(
      title = NULL,
      size = "m",
      footer = tagList(column(12, align="center", 
                              useShinyjs(),
                              div(id = "collecton_source_radio", radioButtons("collection_source_selector", label = NULL, choices = list("New collection" = 1, "Local collection" = 2, "Collection zip" = 3), selected = 1, inline = T)),
                              div(id = "new_collection_name_input", textInput("collection_name", label = "Collection name*")),
                              hidden(div(id = "collection_selector", style = "vertical-align:top; width: 210px;", selectizeInput("selected_collection", label = "Select collection", choices = c(dir("collection_dir/"), ""), selected = ""))),
                              hidden(div(id = "collection_loader", style="vertical-align:top; width: 210px;", fileInput("load_collection", label = "Load collection", accept = ".zip"))),
                              htmlOutput('collection_load_error'),
                              actionButton("load_selected_collecton", label = "Add/Load"), modalButton("Cancel")
      )
      ),
      easyClose = TRUE
    ))
    
  })
  
  observeEvent(input$collection_source_selector, {
    if(input$collection_source_selector == 1) {
      show('new_collection_name_input')
      hide('collection_selector')
      # updateSelectizeInput(session, "selected_collection", selected = "", server = TRUE)
      hide('collection_loader')
    }
    if(input$collection_source_selector == 2) {
      hide('new_collection_name_input')
      # updateTextInput(session, "collection_name", value = "")
      show('collection_selector')
      hide('collection_loader')
    } 
    if(input$collection_source_selector == 3) {
      hide('new_collection_name_input')
      # updateTextInput(session, "collection_name", value = "")
      hide('collection_selector')
      # updateSelectizeInput(session, "selected_collection", selected = "", server = TRUE)
      show('collection_loader')
    }
  })
  
  observeEvent(input$load_selected_collecton, {
    
    v$image_file <- NULL
    v$pathway_name <- NULL
    v$pathway <- NULL
    v$graphnel_df <- NULL
    v$pathway_image <- NULL
    undo_stack <- list()
    redo_stack <- NULL
    
    if(input$collection_source_selector == 1) {
      if(input$collection_name != "") {
        v$collection_name <- input$collection_name
        
        if(!file.exists("collection_dir")) {
          dir.create("collection_dir/")
        }
        
        dir.create(paste0("collection_dir/", input$collection_name))
        dir.create(paste0("collection_dir/", input$collection_name, "/pathways"))
        dir.create(paste0("collection_dir/", input$collection_name, "/pathway_images"))
        
        updateSelectizeInput(session, "selected_pathway", choices = "", label = "Select pathway", server = TRUE)
        show('pathway_loading_widgets')
        hide('curation_tools')
        removeModal()
      } else {
        v$collection_load_error <- "Please speciy name of the new collection"
        delay(3000, v$collection_load_error <- "")
      }
    }
    if(input$collection_source_selector == 2) {
      if(input$selected_collection != "") {
        v$collection_name <- input$selected_collection
        
        collection_pathways <- gsub(".RDS", "", dir(paste0("collection_dir/", v$collection_name, "/pathways")))
        
        if(any(collection_pathways %in% kegg_human_pathway_list$Name)) {
          showTab(session, select = F, inputId = "net_vis_panels", target = "1")
        }
        
        updateSelectizeInput(session, "selected_pathway", choices = collection_pathways, selected = collection_pathways[1], label = "Select pathway", server = TRUE)
        show('pathway_loading_widgets')
        removeModal()
        
        if(length(collection_pathways) == 0) {
          hide('curation_tools')
        }
        
      } else {
        v$collection_load_error <- "Please select collection before loading"
        delay(3000, v$collection_load_error <- "")
      }
    } 
    if(input$collection_source_selector == 3) {
      if(length(input$load_collection$datapath) > 0) {
        #### collection loader ####
        if(grepl("zip", input$load_collection$datapath)) {
          v$collection_name <- gsub("/", "", unzip(zipfile = input$load_collection$datapath, list = TRUE)$Name[1])
          unzip(zipfile = input$load_collection$datapath, exdir = "collection_dir")
          
          collection_pathways <- gsub(".RDS", "", dir(paste0("collection_dir/", v$collection_name, "/pathways")))
          
          if(any(collection_pathways %in% kegg_human_pathway_list$Name)) {
            showTab(session, select = F, inputId = "net_vis_panels", target = "1")
          }
          
          updateSelectizeInput(session, "selected_pathway", choices = collection_pathways, selected = collection_pathways[1], label = "Select pathway", server = TRUE)
          show('pathway_loading_widgets')
        }
        removeModal()
        
        if(length(collection_pathways) == 0) {
          hide('curation_tools')
        }
        
      } else {
        v$collection_load_error <- "Please upload collection before loading"
        delay(3000, v$collection_load_error <- "")
      }
    }
    
  })
  
  output$collection_load_error <- renderText({
    paste("<font color=\"#ff0000\"><b>", v$collection_load_error, "</b></font>")
  })
  
  
  #### New pathway popup ####
  observeEvent(input$add_pathway, {
    
    showModal(modalDialog(
      title = NULL,
      size = "m",
      footer = tagList(column(12, align="center", 
                              useShinyjs(),
                              div(style="vertical-align:top; width: 210px;", selectizeInput("new_pathway_source", label = "Pathway source", choices = c("KEGG", "DF", "Editor"), selected = "KEGG")),
                              div(id = "new_kegg_pathway_selector",
                                  div(style = "display: inline-block; vertical-align:top; width: 210px;",selectizeInput("new_kegg_pathway_organism", label = "Select organism", choices = "hsa", selected = "hsa")),
                                  div(style = "display: inline-block; vertical-align:top; width: 210px;",selectizeInput("new_kegg_pathway", label = "Select new pathway", choices = kegg_human_pathway_list$Name, selected = ""))
                              ),
                              hidden(div(id = "dfs_loader", 
                                         div(style = "display: inline-block; vertical-align:top; width: 210px;", fileInput("df_node_table", label = "Load node table", accept = ".tsv")),
                                         div(style = "display: inline-block; vertical-align:top; width: 210px;", fileInput("df_edge_table", label = "Load edge table", accept = ".tsv")))
                                     ),
                              hidden(div(id = "pathway_name_text_input", 
                                         textInput("new_pathway_name", label = "Pathway name*"))
                                     ),
                              actionButton("add_new_pathway", label = "Add"), modalButton("Cancel")
      )
      ),
      easyClose = TRUE
    ))
    
  })
  
  
  observeEvent(input$new_pathway_source, {
    if(input$new_pathway_source == "KEGG") {
      show('new_kegg_pathway_selector')
      hide('dfs_loader')
      hide('pathway_name_text_input')
    }
    if(input$new_pathway_source == "DF") {
      hide('new_kegg_pathway_selector')
      show('dfs_loader')
      show('pathway_name_text_input')
    }
    if(input$new_pathway_source == "Editor") {
      hide('dfs_loader')
      hide('new_kegg_pathway_selector')
      show('pathway_name_text_input')
    }
  })
  
  #### add new pathway to collection ####
  observeEvent(input$add_new_pathway, {
    
    if(input$new_pathway_source == "KEGG") {
      pathway_data <- kegg_data_downloader(input$new_kegg_pathway, collection_name = v$collection_name)
      saveRDS(pathway_data$pathway, file = paste0("collection_dir/", v$collection_name, "/pathways/", input$new_kegg_pathway, ".RDS"))
      updateSelectizeInput(session, "selected_pathway", choices = gsub(".RDS", "", dir(paste0("collection_dir/", v$collection_name, "/pathways"))), selected = input$new_kegg_pathway, label = "Select pathway", server = TRUE)
      removeModal()
      
    }
    if(input$new_pathway_source == "DF") {
      hideTab(session, select = F, inputId = "net_vis_panels", target = "1")
      
    }
    if(input$new_pathway_source == "Editor") {
      hideTab(session, select = F, inputId = "net_vis_panels", target = "1")
      
    }
    
  })
  
  #### pathway selection by selectize input ####
  observeEvent(input$load_pathway, {
    if(input$selected_pathway != "") {
      v$pathway_name <- input$selected_pathway
      
      v$image_file <- NULL
      
      v$pathway <- readRDS(paste0("collection_dir/", v$collection_name, "/pathways/", v$pathway_name, ".RDS"))
      v$pathway_image <- magick::image_read(paste0("collection_dir/", v$collection_name, "/pathway_images/", pathway_codes_new[v$pathway_name], ".png"))
      
      v$graphnel_df <- graphnel_to_df(v$pathway, extended = TRUE)
      
      hide('vis_buttons')
      hide('pi_and_topology_analysis')
      
      v$undo_stack <- list(v$graphnel_df)
      v$redo_stack <- NULL
      
      v$psf_and_colors <- NULL
      
      # v$pathway_exp_colored <- FALSE
      # v$pathway_psf_colored <- FALSE
      
      v$image_file <- plot_kegg_pathway(graphnel_df = v$graphnel_df, group_graphics = v$pathway$group_nodes, pathway_image = v$pathway_image,
                        node_colors = NULL, edge_mapping = FALSE, edge_in_mode = FALSE, highlight_nodes = NULL, highlight_color = "red",
                        col_legend_title = NULL, color_bar_lims, present_node_modifications = input$show_modified_nodes, removed_nodes = v$pathway$removed_nodes) %>%
        image_write(tempfile(fileext='png'), format = 'png')
      
      show('curation_tools')
      
    }
  })
  
  #### pathway name rendering ####
  output$collection_and_pathway_names <- renderText({
    paste("<font color=\"#1f992f\"><b>", "Collection: ", v$collection_name, "<br>", gsub("_", " ", v$pathway_name), "</b></font>")
  })
  
  #### swithcher for highlighting changed items ####
  ## will come back to this later
  # observeEvent(input$show_changes, {
  #   if(input$show_changes) {
  #     v$graphical_data <- graphical_data_generator(v$pathway, include_changes = TRUE)
  #   } else {
  #     v$graphical_data <- graphical_data_generator(v$pathway)
  #   }
  #   v$image_file <- kegg_node_mapper(group_graphics = v$pathway_data$group_graphics, kegg_pathway_graphics = v$graphical_data, pathway_name = v$pathway_name, pathway_image = v$pathway_image, show_changes = input$show_changes) %>% 
  #     image_write(tempfile(fileext='png'), format = 'png')
  # })
  
  observeEvent(input$show_modified_nodes, {
    if(!is.null(v$graphnel_df)) {
      v$image_file <- plot_kegg_pathway(graphnel_df = v$graphnel_df, group_graphics = v$pathway$group_nodes, pathway_image = v$pathway_image,
                                        node_colors = v$node_colors$node_colors, edge_mapping = FALSE, highlight_nodes = v$selected_node_name$node_id, highlight_color = "red",
                                        col_legend_title = v$node_colors$col_legend_title, color_bar_lims = v$node_colors$color_bar_lims, present_node_modifications = input$show_modified_nodes, removed_nodes = v$pathway$removed_nodes) %>%
        image_write(tempfile(fileext='png'), format = 'png')
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
    
    hide('exp_change_panel')
    
    if(!is.null(v$node_colors)) {
      v$image_file <- plot_kegg_pathway(graphnel_df = v$graphnel_df, group_graphics = v$pathway$group_nodes, pathway_image = v$pathway_image,
                                        node_colors = v$node_colors$node_colors, edge_mapping = FALSE, edge_in_mode = FALSE, highlight_nodes = v$selected_node_name$node_id, highlight_color = "red",
                                        col_legend_title = v$node_colors$col_legend_title, color_bar_lims = v$node_colors$color_bar_lims, present_node_modifications = input$show_modified_nodes, removed_nodes = v$pathway$removed_nodes) %>%
        image_write(tempfile(fileext='png'), format = 'png')
      
    } else {
      
      if(!is.null(v$from_nodes)) {
        
        v$image_file <- plot_kegg_pathway(graphnel_df = v$graphnel_df, group_graphics = v$pathway$group_nodes, pathway_image = v$pathway_image,
                                          node_colors = v$node_colors$node_colors, edge_mapping = FALSE, edge_in_mode = FALSE, highlight_nodes = c(v$from_nodes, v$selected_node_name$node_id), highlight_color = c(rep("green", length(v$from_nodes)), rep("red", nrow(v$selected_node_name))),
                                          col_legend_title = v$node_colors$col_legend_title, color_bar_lims = v$node_colors$color_bar_lims, present_node_modifications = input$show_modified_nodes, removed_nodes = v$pathway$removed_nodes) %>%
          image_write(tempfile(fileext='png'), format = 'png')
        
        
      } else {
        v$image_file <- plot_kegg_pathway(graphnel_df = v$graphnel_df, group_graphics = v$pathway$group_nodes, pathway_image = v$pathway_image,
                                          node_colors = v$node_colors$node_colors, edge_mapping = FALSE, edge_in_mode = FALSE, highlight_nodes = v$selected_node_name$node_id, highlight_color = "red",
                                          col_legend_title = v$node_colors$col_legend_title, color_bar_lims = v$node_colors$color_bar_lims, present_node_modifications = input$show_modified_nodes, removed_nodes = v$pathway$removed_nodes) %>%
          image_write(tempfile(fileext='png'), format = 'png')
      }
      
    }
    
    
  })
  
  #### image brush events ####
  observeEvent(input$image_brush, {
    if(!input$event_node_mode) {
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
      
      hide('exp_change_panel')
      
      if(!is.null(v$node_colors)) {
        
        v$image_file <- plot_kegg_pathway(graphnel_df = v$graphnel_df, group_graphics = v$pathway$group_nodes, pathway_image = v$pathway_image,
                                          node_colors = v$node_colors$node_colors, edge_mapping = FALSE, edge_in_mode = FALSE, highlight_nodes = v$selected_node_name$node_id, highlight_color = "red",
                                          col_legend_title = v$node_colors$col_legend_title, color_bar_lims = v$node_colors$color_bar_lims, present_node_modifications = input$show_modified_nodes, removed_nodes = v$pathway$removed_nodes) %>%
          image_write(tempfile(fileext='png'), format = 'png')
        
      } else {
        
        if(!is.null(v$from_nodes)) {
          
          v$image_file <- plot_kegg_pathway(graphnel_df = v$graphnel_df, group_graphics = v$pathway$group_nodes, pathway_image = v$pathway_image,
                                            node_colors = v$node_colors$node_colors, edge_mapping = FALSE, edge_in_mode = FALSE, highlight_nodes = c(v$from_nodes, v$selected_node_name$node_id), highlight_color = c(rep("green", length(v$from_nodes)), rep("red", nrow(v$selected_node_name))),
                                            col_legend_title = v$node_colors$col_legend_title, color_bar_lims = v$node_colors$color_bar_lims, present_node_modifications = input$show_modified_nodes, removed_nodes = v$pathway$removed_nodes) %>%
            image_write(tempfile(fileext='png'), format = 'png')
          
          
        } else {
          v$image_file <- plot_kegg_pathway(graphnel_df = v$graphnel_df, group_graphics = v$pathway$group_nodes, pathway_image = v$pathway_image,
                                            node_colors = v$node_colors$node_colors, edge_mapping = FALSE, edge_in_mode = FALSE, highlight_nodes = v$selected_node_name$node_id, highlight_color = "red",
                                            col_legend_title = v$node_colors$col_legend_title, color_bar_lims = v$node_colors$color_bar_lims, present_node_modifications = input$show_modified_nodes, removed_nodes = v$pathway$removed_nodes) %>%
            image_write(tempfile(fileext='png'), format = 'png')
        }
        
      }
      
    }
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
  
  
  #### edge drawing ####
  observeEvent(input$draw_edge, {
    
    if(v$allow_edge_draw & !is.null(v$graphnel_df)) {
      
      if(!is.null(v$node_colors)) {
        
        v$image_file <- plot_kegg_pathway(graphnel_df = v$graphnel_df, group_graphics = v$pathway$group_nodes, pathway_image = v$pathway_image,
                                          node_colors = v$node_colors$node_colors, edge_mapping = TRUE, edge_in_mode = input$ingoing_edge, highlight_nodes = v$selected_node_name$node_id, highlight_color = "red",
                                          col_legend_title = v$node_colors$col_legend_title, color_bar_lims = v$node_colors$color_bar_lims, present_node_modifications = input$show_modified_nodes, removed_nodes = v$pathway$removed_nodes) %>%
          image_write(tempfile(fileext='png'), format = 'png')
        
      } else {
        v$image_file <- plot_kegg_pathway(graphnel_df = v$graphnel_df, group_graphics = v$pathway$group_nodes, pathway_image = v$pathway_image,
                                          node_colors = v$node_colors$node_colors, edge_mapping = TRUE, edge_in_mode = input$ingoing_edge, highlight_nodes = v$selected_node_name$node_id, highlight_color = "red",
                                          col_legend_title = v$node_colors$col_legend_title, color_bar_lims = v$node_colors$color_bar_lims, present_node_modifications = input$show_modified_nodes, removed_nodes = v$pathway$removed_nodes) %>%
          image_write(tempfile(fileext='png'), format = 'png')
      }
    }
    
  })
  
  #### disconnected node(s) checker ####
  observeEvent(input$check_disconnected_nodes, {
    disconnected_nodes <- setdiff(v$graphnel_df$node_table$node_id, unique(c(v$graphnel_df$edge_table$from, v$graphnel_df$edge_table$to)))
    
    visNetworkProxy("visnet") %>%
      visSelectNodes(id = disconnected_nodes, highlightEdges = F)
    
    v$image_file <- plot_kegg_pathway(graphnel_df = v$graphnel_df, group_graphics = v$pathway$group_nodes, pathway_image = v$pathway_image,
                                      node_colors = v$node_colors$node_colors, edge_mapping = FALSE, edge_in_mode = FALSE, highlight_nodes = disconnected_nodes, highlight_color = "red",
                                      col_legend_title = v$node_colors$col_legend_title, color_bar_lims = v$node_colors$color_bar_lims, node_fill_opacity = 0.5, present_node_modifications = input$show_modified_nodes, removed_nodes = v$pathway$removed_nodes) %>%
      image_write(tempfile(fileext='png'), format = 'png')
  })
  
  #### duplicated node(s) checker ####
  observeEvent(input$check_duplicated_nodes, {
    duplicated_labels <- v$graphnel_df$node_table$label[which(duplicated(v$graphnel_df$node_table$label))]
    
    duplicated_nodes <- v$graphnel_df$node_table$node_id[which(v$graphnel_df$node_table$label %in% duplicated_labels)]
    
    visNetworkProxy("visnet") %>%
      visSelectNodes(id = duplicated_nodes, highlightEdges = F)
    
    v$image_file <- plot_kegg_pathway(graphnel_df = v$graphnel_df, group_graphics = v$pathway$group_nodes, pathway_image = v$pathway_image,
                                      node_colors = v$node_colors$node_colors, edge_mapping = FALSE, edge_in_mode = FALSE, highlight_nodes = duplicated_nodes, highlight_color = "red",
                                      col_legend_title = v$node_colors$col_legend_title, color_bar_lims = v$node_colors$color_bar_lims, node_fill_opacity = 0.5, present_node_modifications = input$show_modified_nodes, removed_nodes = v$pathway$removed_nodes) %>%
      image_write(tempfile(fileext='png'), format = 'png')
    
  })
  
  #### clear image marks ####
  ## previous observer activator c(input$app_mode, input$clear)
  observeEvent(input$clear, {
    if(!is.null(v$graphnel_df)) {
      removeUI(selector = '#edge_search')
      hide('edge_search_dialog')
      hide('add_edge_button')
      hide('transition_maker_button')
      hide('edge_delete_button')
      hide('dir_switch')
      hide('delete_nodes')
      hide('edit_edge')
      show('add_new_node')
      hide('edge_panel')
      hide('edge_creating_panel')
      hide("event_node_adding_panel")
      hide('edge_attr_editing_panel_vis')
      hide('exp_change_panel')
      v$allow_edge_draw <- TRUE
      updateCheckboxInput(session, "event_node_mode", value = FALSE)
      v$selected_node_name <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("node_id", "component_id_s", "label"))
      v$from_nodes <- NULL
      
      v$node_colors <- NULL
      v$graphnel_df$node_table$expression <- 1
      v$graphnel_df$node_table$signal <- 1
      
      v$image_file <- plot_kegg_pathway(graphnel_df = v$graphnel_df, group_graphics = v$pathway$group_nodes, pathway_image = v$pathway_image,
                                        node_colors = v$node_colors$node_colors, edge_mapping = FALSE, edge_in_mode = FALSE, highlight_nodes = NULL, highlight_color = "red",
                                        col_legend_title = v$node_colors$col_legend_title, color_bar_lims = v$node_colors$color_bar_lims, present_node_modifications = input$show_modified_nodes, removed_nodes = v$pathway$removed_nodes) %>%
        image_write(tempfile(fileext='png'), format = 'png')
      
      ### continue from here
      # v$visnet_list <- visnet_creator(v$graphical_data)
    }
    
  }, ignoreInit = TRUE)
  
  #### right click dialog ####
  observeEvent(input$right_click, {
    
    if(input$app_mode == "Curation") {
      hide('change_exp_value')
      hide('add_drug_button')
      hide('edge_creating_panel')
      hide('change_edge_weight_button')
      hide('color_scale_fixer')
      if(nrow(v$selected_node_name) == 0) {
        hide('edge_search_dialog')
        hide('add_edge_button')
        hide('event_association')
        hide('transition_maker_button')
        hide('edge_delete_button')
        hide('dir_switch')
        hide('delete_nodes')
        hide('edit_edge')
        hide('connect_to_button')
      } else {
        if(nrow(v$selected_node_name) == 1) {
          show('edge_search_dialog')
          show('connect_to_button')
          show('transition_maker_button')
          show('event_association')
          show('delete_nodes')
          hide('add_new_node')
          hide('add_edge_button')
          hide('dir_switch')
          hide('edge_delete_button')
          hide('edit_edge')
        } else {
          if(nrow(v$selected_node_name) == 2) {
            if(nrow(as.data.table(v$graphnel_df$edge_table)[from == v$selected_node_name[1,"node_id"] & to == v$selected_node_name[2,"node_id"] | from == v$selected_node_name[2,"node_id"] & to == v$selected_node_name[1,"node_id"]]) == 0) {
              show('add_edge_button')
              hide('dir_switch')
              hide('edge_delete_button')
              hide('edit_edge')
            } else {
              if(nrow(as.data.table(v$graphnel_df$edge_table)[from == v$selected_node_name[1,"node_id"] & to == v$selected_node_name[2,"node_id"] | from == v$selected_node_name[2,"node_id"] & to == v$selected_node_name[1,"node_id"]]) == 1) {
                hide('add_edge_button')
                show('dir_switch')
                show('edge_delete_button')
                show('edit_edge')
              } else {
                hide('add_edge_button')
                hide('dir_switch')
                hide('edit_edge')
                show('edge_delete_button')
              }
            }
            
            show('connect_to_button')
            show('edge_search_dialog')
            show('transition_maker_button')
            show('event_association')
            show('delete_nodes')
            hide('add_new_node')
          } else {
            hide('edge_search_dialog')
            hide('add_edge_button')
            hide('edge_delete_button')
            hide('dir_switch')
            hide('add_new_node')
            hide('edit_edge')
            show('connect_to_button')
            show('transition_maker_button')
            show('event_association')
            show('delete_nodes')
          }
        }
      }
      
    } else {
      
      hide('connect_to_button')
      hide('edge_search_dialog')
      hide('add_edge_button')
      hide('event_association')
      hide('transition_maker_button')
      hide('edge_delete_button')
      hide('dir_switch')
      hide('delete_nodes')
      hide('edit_edge')
      hide('color_scale_fixer')
      hide('change_edge_weight_button')
      hide('change_exp_value')
      hide('add_drug_button')
      
      if(v$allow_graph_param_update) {
        if(nrow(v$selected_node_name) == 1) {
          show('change_exp_value')
          show('add_drug_button')
        }
        
        if(nrow(v$selected_node_name) == 2) {
          
          if(paste(v$selected_node_name[c(1,2),"node_id"], collapse = ";") %in% paste(v$graphnel_df$edge_table$from, v$graphnel_df$edge_table$to, sep = ";") |
             paste(v$selected_node_name[c(2,1),"node_id"], collapse = ";") %in% paste(v$graphnel_df$edge_table$from, v$graphnel_df$edge_table$to, sep = ";")
          ) {
            show('change_edge_weight_button')
          }
        }
      }
      
    }
    
    if(!is.null(v$graphnel_df)) {
      show('dialog_content')
    }
    
  })
  
  observeEvent(input$left_click, {
    hide('dialog_content')
  })
  
  #### interaction table output ####
  output$edge_table <- DT::renderDataTable({datatable(
    database_table_design(v$database_out),
    filter = 'bottom',
    class   = 'cell-border compact hover',
    fillContainer = T,
    escape = F,
    selection = list(mode = 'single'),
    options = list(dom = 'rtp',
                   pageLength = 100,
                   searchHighlight = TRUE,
                   columnDefs = list(list(
                     targets = 6,
                     render   = JS(
                       "function(data, type, row, meta) {",
                       "return type === 'display' && data.length > 10 ?",                                    "'<span title=\"' + data + '\">' +
                                data.substr(0, 20) + '...</span>' : data;", "}")
                   )),
                   initComplete = JS("function(settings, json) {",
                                     "$(this.api().table().header()).css({'background-color': '#3474B7', 'color': '#fff'});",
                                     "}"),
                   scrollY = 250,scrollX = TRUE
    ),
    style = 'bootstrap', editable = FALSE)
  })
  
  #### pair search in pathway databases ####
  observeEvent(input$edge_search_dialog, {
    
    updateCheckboxInput(session, inputId = "direction_state", value = FALSE)
    v$searching_state <- TRUE
    v$database_search_node <- v$selected_node_name
    hide('add_selected_edge')
    show('edge_panel')
    
    if(nrow(v$database_search_node) == 1) {
      
      v$database_out <- node_pairs_generator(v$database_search_node, pathway_realted = TRUE, pathway_data = v$graphnel_df$node_table, direction = input$direction_state, direction_state = input$direction_side)
      
    } else {
      
      v$database_out <- node_pairs_generator(v$database_search_node, direction = input$direction_state, direction_state = input$direction_side)
      
    }
    
  })
  
  observe({
    if(!is.null(v$database_search_node)) {
      if(!all(v$database_search_node[,"node_id"] %in% v$selected_node_name[,"node_id"])) {
        v$searching_state <- FALSE
        hide('edge_panel')
      }
    }
  })
  
  observe({
    if(!is.null(v$database_edit_node)) {
      if(!all(v$database_edit_node[,"node_id"] %in% v$selected_node_name[,"node_id"])) {
        hide('edge_creating_panel')
      }
    }
  })
  
  observeEvent(input$close_table, {
    v$searching_state <- FALSE
    hide('edge_panel')
  })
  
  observeEvent(input$selected_pathway, {
    v$searching_state <- FALSE
    hide('edge_panel')
    hide('edge_creating_panel')
  })
  
  #### database data updater based on direction switcher ####
  observe({
    if(v$searching_state) {
      if(nrow(v$database_search_node) == 1) {
        v$database_out <- node_pairs_generator(v$database_search_node, pathway_realted = TRUE, pathway_data = v$graphnel_df$node_table, direction = input$direction_state, direction_state = input$direction_side)
      } else {
        v$database_out <- node_pairs_generator(v$database_search_node, direction = input$direction_state, direction_state = input$direction_side)
      }
    }
  })
  
  #### database table row selection ####
  observeEvent(input$edge_table_rows_selected, {
    if(!is.null(input$edge_table_rows_selected)) {
      if(v$searching_state) {
        
        show('add_selected_edge')
        if(input$edge_table_rows_selected != v$previous_row_selection) {
          highlight_genes <- rbind(v$graphnel_df$node_table[which(v$graphnel_df$node_table$node_id %in% v$database_out[input$edge_table_rows_selected, 9]), c("node_id", "component_id_s")],
                                   v$graphnel_df$node_table[which(v$graphnel_df$node_table$node_id %in% v$database_out[input$edge_table_rows_selected, 10]), c("node_id", "component_id_s")])
          
          v$image_file <- plot_kegg_pathway(graphnel_df = v$graphnel_df, group_graphics = v$pathway$group_nodes, pathway_image = v$pathway_image,
                                            node_colors = v$node_colors$node_colors, custom_edge_mapping = TRUE, edge_mapping = TRUE, edge_in_mode = TRUE, highlight_nodes = highlight_genes$node_id, highlight_color = "red",
                                            col_legend_title = v$node_colors$col_legend_title, color_bar_lims = v$node_colors$color_bar_lims, present_node_modifications = input$show_modified_nodes, removed_nodes = v$pathway$removed_nodes) %>%
            image_write(tempfile(fileext='png'), format = 'png')
          
          v$previous_row_selection <- input$edge_table_rows_selected
        }
        
      } else {
        hide('add_selected_edge')
      }
    } else {
      hide('add_selected_edge')
    }
  }, ignoreNULL = F)
  
  #### adding new edge from database ####
  observeEvent(input$add_selected_edge,{
    v$edge_edit_stage <- FALSE
    updateActionButton(session, "create_edge", label = "Create edge")
    hide('edge_panel')
    show('edge_creating_panel')
    
    updateSelectizeInput(session, "interaction_source", selected = "database")
    updateSelectizeInput(session, "subtype1", selected = "")
    updateSelectizeInput(session, "subtype2", selected = "")
    
    v$selected_interaction <- list(source = paste(v$database_out[input$edge_table_rows_selected, "INTERACTION_DATA_SOURCE"], "database"), selected_interaction = rbind(v$graphnel_df$node_table[which(v$graphnel_df$node_table$node_id %in% v$database_out[input$edge_table_rows_selected, 9]), c("label", "node_id")],
                                                                                                                                                                       v$graphnel_df$node_table[which(v$graphnel_df$node_table$node_id %in% v$database_out[input$edge_table_rows_selected, 10]), c("label", "node_id")]))
  })
  
  #### adding new edge from image ####
  observeEvent(input$add_edge, {
    v$edge_edit_stage <- FALSE
    updateActionButton(session, "create_edge", label = "Create edge")
    updateSelectizeInput(session, "interaction_source", selected = "image")
    updateSelectizeInput(session, "subtype1", selected = "")
    updateSelectizeInput(session, "subtype2", selected = "")
    show('edge_creating_panel')
    show('image_edge_direction_switch')
    
    v$selected_interaction <- list(source = "database", selected_interaction = v$selected_node_name[,c("label", "node_id")])
    
    v$image_file <- plot_kegg_pathway(graphnel_df = v$graphnel_df, group_graphics = v$pathway$group_nodes, pathway_image = v$pathway_image,
                                      node_colors = v$node_colors$node_colors, custom_edge_mapping = TRUE, edge_mapping = TRUE, edge_in_mode = TRUE, highlight_nodes = v$selected_node_name$node_id, highlight_color = "red",
                                      col_legend_title = v$node_colors$col_legend_title, color_bar_lims = v$node_colors$color_bar_lims, present_node_modifications = input$show_modified_nodes, removed_nodes = v$pathway$removed_nodes) %>%
      image_write(tempfile(fileext='png'), format = 'png')
    
  })
  
  #### many to many edge adder ####
  observeEvent(input$connect_to, {
    v$from_nodes <- v$selected_node_name$node_id
    
    v$image_file <- plot_kegg_pathway(graphnel_df = v$graphnel_df, group_graphics = v$pathway$group_nodes, pathway_image = v$pathway_image,
                                      node_colors = v$node_colors$node_colors, highlight_nodes = v$selected_node_name$node_id, highlight_color = "green",
                                      col_legend_title = v$node_colors$col_legend_title, color_bar_lims = v$node_colors$color_bar_lims, present_node_modifications = input$show_modified_nodes, removed_nodes = v$pathway$removed_nodes) %>%
      image_write(tempfile(fileext='png'), format = 'png')
    
    updateSelectizeInput(session, "interaction_source", selected = "image")
    show('edge_creating_panel')
    show('node_selecion_message')
  })
  
  #### edge direction switcher for image based edge adding ####
  observeEvent(input$image_edge_direction, {
    if(input$image_edge_direction) {
      v$selected_interaction <- list(source = "database", selected_interaction = v$selected_node_name[,c("label", "node_id")])
      v$image_file <- plot_kegg_pathway(graphnel_df = v$graphnel_df, group_graphics = v$pathway$group_nodes, pathway_image = v$pathway_image,
                                        node_colors = v$node_colors$node_colors, custom_edge_mapping = TRUE, edge_mapping = TRUE, edge_in_mode = TRUE, highlight_nodes = v$selected_node_name$node_id, highlight_color = "red",
                                        col_legend_title = v$node_colors$col_legend_title, color_bar_lims = v$node_colors$color_bar_lims, present_node_modifications = input$show_modified_nodes, removed_nodes = v$pathway$removed_nodes) %>%
        image_write(tempfile(fileext='png'), format = 'png')
    } else {
      v$selected_interaction <- list(source = "database", selected_interaction = v$selected_node_name[2:1,c("label", "node_id")])
      v$image_file <- plot_kegg_pathway(graphnel_df = v$graphnel_df, group_graphics = v$pathway$group_nodes, pathway_image = v$pathway_image,
                                        node_colors = v$node_colors$node_colors, custom_edge_mapping = TRUE, edge_mapping = TRUE, edge_in_mode = TRUE, highlight_nodes = v$selected_node_name$node_id[2:1], highlight_color = "red",
                                        col_legend_title = v$node_colors$col_legend_title, color_bar_lims = v$node_colors$color_bar_lims, present_node_modifications = input$show_modified_nodes, removed_nodes = v$pathway$removed_nodes) %>%
        image_write(tempfile(fileext='png'), format = 'png')
    }
  }, ignoreInit = T)
  
  #### editing edge attributes ####
  observeEvent(input$edit_edge, {
    v$database_edit_node <- v$selected_node_name
    v$edge_edit_stage <- TRUE
    selected_edge <- as.data.table(v$graphnel_df$edge_table)[from == v$selected_node_name[1,"node_id"] & to == v$selected_node_name[2,"node_id"] | from == v$selected_node_name[2,"node_id"] & to == v$selected_node_name[1,"node_id"]]
    
    updateSelectizeInput(session, "subtype1", selected = selected_edge$subtype1) # unname(unlist(graph::edgeData(v$pathway$grap, selected_edge$from, selected_edge$to, attr = "subtype1"))))
    updateSelectizeInput(session, "subtype2", selected = selected_edge$subtype2) # unname(unlist(graph::edgeData(v$pathway$grap, selected_edge$from, selected_edge$to, attr = "subtype2"))))
    updateSelectizeInput(session, "interaction_source", selected = "")
    updateActionButton(session, "create_edge", label = "Edit edge")
    show('edge_creating_panel')
    
    selected_interaction <- rbind(v$graphnel_df$node_table[which(v$graphnel_df$node_table$node_id == selected_edge$from), c("label", "node_id")],
                                  v$graphnel_df$node_table[which(v$graphnel_df$node_table$node_id == selected_edge$to), c("label", "node_id")])
    
    v$selected_interaction <- list(source = "", selected_interaction = selected_interaction)
    
  })
  
  output$edge_info <- renderText({
    paste("<b>", paste0(v$selected_interaction$selected_interaction[1,"label"], " to ", v$selected_interaction$selected_interaction[2,"label"]), "</b>")
  })
  
  observeEvent(input$close_panel, {
    hide('edge_creating_panel')
    hide('image_edge_direction_switch')
  })
  
  output$edge_adding_error <- renderText({
    paste("<font color=\"#ff0000\"><b>", v$edge_adding_error, "</b></font>")
  })
  
  #### edge assigning to graph and editing existed one ####
  observeEvent(input$create_edge, {
    
    if(input$subtype1 == "" | is.null(input$interaction_source)) {
      v$edge_adding_error <- "Fill all necessary fields"
      delay(3000, v$edge_adding_error <- "")
    } else {
      edge_id <- paste0(v$selected_interaction$selected_interaction[1,"node_id"], "|", v$selected_interaction$selected_interaction[2,"node_id"])
      if(v$edge_edit_stage) {
        
        v$graphnel_df$edge_table[edge_id, "subtype1"] <-  input$subtype1
        v$graphnel_df$edge_table[edge_id, "subtype2"] <-  input$subtype2
        v$graphnel_df$edge_table[edge_id, "color"] <- line_col[input$subtype1]
        v$graphnel_df$edge_table[edge_id, "arrows.to.type"] <- kegg_arrows_type[input$subtype1]
        v$graphnel_df$edge_table[edge_id, "dashes"] <- input$subtype2 == "indirect effect" | input$subtype1 == "indirect effect"
        v$graphnel_df$edge_table[edge_id, "change_info"] <-  "attr_change"
        v$graphnel_df$edge_table[edge_id, "data_source"] <-  paste(input$interaction_source, collapse = ", ")
        
        v$pathway_change_alert <- v$pathway_change_alert + 1
        hide('edge_creating_panel')
      } else {
        
        if(!is.null(v$from_nodes)) {
          
          new_edges <- Reduce(rbind, lapply(v$from_nodes, function(x) {
            
            data.frame(
              id = paste0(rep(x, length(v$selected_node_name$node_id)), "|", v$selected_node_name$node_id),
              from = rep(x, length(v$selected_node_name$node_id)),
              to = v$selected_node_name$node_id,
              color = line_col[input$subtype1],
              arrows.to.enabled = TRUE,
              arrows.to.type = kegg_arrows_type[input$subtype1],
              label = "",
              dashes = input$subtype2 == "indirect effect" | input$subtype1 == "indirect effect",
              type = "",
              subtype1 = input$subtype1,
              subtype2 = input$subtype2,
              state = "",
              weight = 1,
              existence = "exist",
              change_info = "added",
              data_source = "image",
              row.names = paste0(rep(x, length(v$selected_node_name$node_id)), "|", v$selected_node_name$node_id)
            )
            
          }))
          
          v$graphnel_df$edge_table[rownames(new_edges),] <- new_edges
          hide('node_selecion_message')
          v$from_nodes <- NULL
        } else {
          data_source <- input$interaction_source
          data_source[which(data_source == "database")] <- v$selected_interaction$source
          data_source <- paste(data_source, collapse = ", ")
          if(!(edge_id %in% rownames(v$graphnel_df$edge_table))) {
            
            new_edge <- data.frame(
              id = edge_id,
              from = v$selected_interaction$selected_interaction[1,"node_id"],
              to = v$selected_interaction$selected_interaction[2,"node_id"],
              color = line_col[input$subtype1],
              arrows.to.enabled = TRUE,
              arrows.to.type = kegg_arrows_type[input$subtype1],
              label = "",
              dashes = input$subtype2 == "indirect effect" | input$subtype1 == "indirect effect",
              type = "",
              subtype1 = input$subtype1,
              subtype2 = input$subtype2,
              state = "",
              weight = 1,
              existence = "exist",
              change_info = "added",
              data_source = data_source,
              row.names = edge_id
            )
            
            v$graphnel_df$edge_table <- rbind(v$graphnel_df$edge_table, new_edge)
          }
        }
        
        v$pathway_change_alert <- v$pathway_change_alert + 1
        hide('edge_creating_panel')
        hide('image_edge_direction_switch')
      }
    }
    
  })
  
  #### edge direction switch ####
  observeEvent(input$edge_direction_switch, {
    index <- unique(which(v$graphnel_df$edge_table$from %in% v$selected_node_name[,"node_id"] & v$graphnel_df$edge_table$to %in% v$selected_node_name[,"node_id"]))
    
    from <- v$graphnel_df$edge_table$from[index]
    to <- v$graphnel_df$edge_table$to[index]
    
    v$graphnel_df$edge_table$from[index] <- to
    v$graphnel_df$edge_table$to[index] <- from
    
    v$pathway_change_alert <- v$pathway_change_alert + 1
  })
  
  #### delete edge ####
  observeEvent(input$edge_delete, {
    index <- unique(which(v$graphnel_df$edge_table$from %in% v$selected_node_name[,"node_id"] & v$graphnel_df$edge_table$to %in% v$selected_node_name[,"node_id"]))
    
    v$graphnel_df$edge_table[index, "existence"] <- "removed"
    
    v$pathway_change_alert <- v$pathway_change_alert + 1
  })
  
  output$node_adding_error <- renderText({
    paste("<font color=\"#ff0000\"><b>", v$node_adding_error, "</b></font>")
  })
  
  #### event type selectize input updater
  observe({
    updateSelectizeInput(session, "event_node_name", choices = v$event_node_type, selected = NULL, server = T)
  })
  
  #### event association ####
  ## will come to this later
  # observeEvent(input$associate_with_event_node, {
  #   if("bio_event" %in% unlist(graph::nodeData(v$pathway$graph, v$selected_node_name[,"node_id"], attr = "type"))) {
  #     event_ids <- names(which(unlist(graph::nodeData(v$pathway$graph, v$selected_node_name[,"node_id"], attr = "type")) == "bio_event"))
  #     gene_ids <- names(which(unlist(graph::nodeData(v$pathway$graph, v$selected_node_name[,"node_id"], attr = "type")) != "bio_event"))
  #     
  #     if(length(event_ids) > 0 & length(event_ids) < 2 & length(gene_ids) > 0) {
  #       v$pathway$graph <- graph::addEdge(gene_ids, event_ids, v$pathway$graph)
  #       graph::edgeData(v$pathway$graph, gene_ids, event_ids, attr = "data_source") <- "image"
  #       graph::edgeData(v$pathway$graph, gene_ids, event_ids, attr = "subtype1") <- "activation"
  #       graph::edgeData(v$pathway$graph, gene_ids, event_ids, attr = "subtype2") <- "indirect effect"
  #       
  #       v$pathway_change_alert <- v$pathway_change_alert + 1
  #     }
  #   } else {
  #     v$allow_edge_draw <- FALSE
  #     show("event_node_adding_panel")
  #   }
  #   
  # })
  # 
  # #### adding new node ####
  # observeEvent(input$add_node, {
  #   
  #   if(is.null(input$image_brush[1:4])) {
  #     v$node_adding_error <- "Please specify node location"
  #     delay(3000, v$node_adding_error <- NULL)
  #   } else {
  #     new_node_num <- as.character(max(as.integer(v$pathway$graph@nodes)) + 1)
  #     v$pathway$graph <- graph::addNode(new_node_num, v$pathway$graph)
  #     graph::nodeData(v$pathway$graph, new_node_num, attr = "label") <- input$event_node_name
  #     graph::nodeData(v$pathway$graph, new_node_num, attr = "kegg.name") <- input$event_node_name
  #     graph::nodeData(v$pathway$graph, new_node_num, attr = "data_source") <- input$node_source
  #     graph::nodeData(v$pathway$graph, new_node_num, attr = "type") <- "bio_event"
  #     graph::nodeData(v$pathway$graph, new_node_num, attr = "kegg.id") <- new_node_num
  #     graph::nodeData(v$pathway$graph, new_node_num, attr = "kegg.gr.x") <- input$image_brush[[1]] + (input$image_brush[[2]] - input$image_brush[[1]])*0.5
  #     graph::nodeData(v$pathway$graph, new_node_num, attr = "kegg.gr.y") <- input$image_brush[[3]] + (input$image_brush[[4]] - input$image_brush[[3]])*0.5
  #     graph::nodeData(v$pathway$graph, new_node_num, attr = "kegg.gr.width") <- input$image_brush[[2]] - input$image_brush[[1]]
  #     graph::nodeData(v$pathway$graph, new_node_num, attr = "kegg.gr.height") <- input$image_brush[[4]] - input$image_brush[[3]]
  #     
  #     print(v$selected_node_name[,"node_id"])
  #     v$pathway$graph <- graph::addEdge(v$selected_node_name[,"node_id"], new_node_num, v$pathway$graph)
  #     graph::edgeData(v$pathway$graph, v$selected_node_name[,"node_id"], new_node_num, attr = "data_source") <- input$node_source
  #     graph::edgeData(v$pathway$graph, v$selected_node_name[,"node_id"], new_node_num, attr = "subtype1") <- "activation"
  #     graph::edgeData(v$pathway$graph, v$selected_node_name[,"node_id"], new_node_num, attr = "subtype2") <- "indirect effect"
  #     
  #     v$pathway_change_alert <- v$pathway_change_alert + 1
  #     hide('event_node_adding_panel')
  #     v$allow_edge_draw <- TRUE
  #     updateCheckboxInput(session, "event_node_mode", value = FALSE)
  #   }
  #   
  #   if(!(input$event_node_name %in% v$event_node_type)) {
  #     v$event_node_type <- unique(c(scan(file = "event_names.txt", what = "character", sep = "\t"), input$event_node_name))
  #     write(v$event_node_type, file = "event_names.txt", sep = "\t")
  #   }
  #   
  # })
  
  observeEvent(input$close_node_panel, {
    hide('event_node_adding_panel')
    v$allow_edge_draw <- TRUE
    updateCheckboxInput(session, "event_node_mode", value = FALSE)
  })
  
  #### make node as a transition node ####
  observeEvent(input$make_transition, {
    
    v$graphnel_df$node_table[v$selected_node_name[,"node_id"],"type"] <- "transition_gene"
    v$graphnel_df$node_table[v$selected_node_name[,"node_id"],"change_info"] <- "class_change"
    v$graphnel_df$node_table[v$selected_node_name[,"node_id"],"data_source"] <- "image"
    
    v$pathway_change_alert <- v$pathway_change_alert + 1
  })
  
  ## continue from here
  #### delete node ####
  observeEvent(input$node_delete, {
    
    index <- unique(which(v$graphnel_df$edge_table$from %in% v$selected_node_name[,"node_id"] | v$graphnel_df$edge_table$to %in% v$selected_node_name[,"node_id"]))
    
    v$graphnel_df$edge_table[index, "existence"] <- "removed"
    v$graphnel_df$node_table[v$selected_node_name[,"node_id"], "existence"] <- "removed"
    
    v$pathway_change_alert <- v$pathway_change_alert + 1
  })
  
  
  ## will come to this later
  #### undo redo observer ####
  observeEvent(v$pathway_change_alert, {
    if(v$pathway_change_alert != 0) {
      if(length(v$undo_stack) > 10) {
        v$undo_stack <- list(v$undo_stack, v$graphnel_df)[[-1]]
      } else {
        if(length(v$undo_stack) == 1) {
          v$undo_stack <- list(v$undo_stack[[1]], v$graphnel_df)
        } else {
          v$undo_stack <- c(v$undo_stack, list(v$graphnel_df))
        }
      }
      v$redo_stack <- NULL
    }
    # View(v$undo_stack)
  }, suspended = F)
  
  observeEvent(input$undo, {
    if(length(v$undo_stack) > 1) {
      v$redo_stack <- c(v$redo_stack, list(v$undo_stack[[length(v$undo_stack)]]))
      v$graphnel_df <- v$undo_stack[[length(v$undo_stack) - 1]]
      v$undo_stack <- v$undo_stack[-length(v$undo_stack)]
      v$undo_alert <- v$undo_alert + 1
    }
  })
  
  observeEvent(input$redo, {
    if(!is.null(v$redo_stack)) {
      if(length(v$redo_stack) > 0) {
        v$graphnel_df <- v$redo_stack[[length(v$redo_stack)]]
        v$undo_stack <- c(v$undo_stack, list(v$graphnel_df))
        v$redo_stack <- v$redo_stack[-length(v$redo_stack)]
        v$undo_alert <- v$undo_alert + 1
      }
    }
  })
  
  #### observer for saving pathway graph after any changes ####
  observeEvent(v$pathway_change_alert, {
    if(v$pathway_change_alert != 0) {
      
      removed_nodes <- v$graphnel_df$node_table[which(v$graphnel_df$node_table$existence == "removed"),]
      removed_edges <- v$graphnel_df$edge_table[which(v$graphnel_df$edge_table$existence == "removed"),]
      
      v$graphnel_df$edge_table <- v$graphnel_df$edge_table[which(v$graphnel_df$edge_table$existence != "removed"),]
      v$graphnel_df$node_table <- v$graphnel_df$node_table[which(v$graphnel_df$node_table$existence != "removed"),]
      
      v$graphnel_df$node_table$sink <- FALSE
      v$graphnel_df$node_table[unique(v$graphnel_df$edge_table$to[which(!(v$graphnel_df$edge_table$to %in% v$graphnel_df$edge_table$from))]), "sink"] <- TRUE
      v$graphnel_df$node_table$color.background <- ifelse(v$graphnel_df$node_table$sink, "#0099cc", "#BFFFBF")
      
      v$pathway[c("sink.nodes", "order", "graph")] <- df_to_graphnel(v$graphnel_df$node_table, v$graphnel_df$edge_table)[c("sink.nodes", "order", "graph")]
      
      v$pathway$graphnel_df <- v$graphnel_df
      
      v$pathway$removed_nodes <- rbind(v$pathway$removed_nodes, removed_nodes)
      v$pathway$removed_edges <- rbind(v$pathway$removed_edges, removed_edges)
      
      v$selected_node_name <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("node_id", "component_id_s", "label"))
      
      v$image_file <- plot_kegg_pathway(graphnel_df = v$graphnel_df, group_graphics = v$pathway$group_nodes, pathway_image = v$pathway_image,
                                        node_colors = v$node_colors$node_colors, edge_mapping = FALSE, highlight_nodes = v$selected_node_name$node_id, highlight_color = "red",
                                        col_legend_title = v$node_colors$col_legend_title, color_bar_lims = v$node_colors$color_bar_lims, present_node_modifications = input$show_modified_nodes, removed_nodes = v$pathway$removed_nodes) %>%
        image_write(tempfile(fileext='png'), format = 'png')
      
      saveRDS(v$pathway, file = paste0("collection_dir/", v$collection_name, "/pathways/", v$pathway_name, ".RDS"))
    }
  }, suspended = F)
  
  #### observer for saving pathway graph after undo or redo operations ####
  observeEvent(v$undo_alert, {
    if(v$undo_alert != 0) {
      
      removed_nodes <- v$graphnel_df$node_table[which(v$graphnel_df$node_table$existence == "removed"),]
      removed_edges <- v$graphnel_df$edge_table[which(v$graphnel_df$edge_table$existence == "removed"),]
      
      v$graphnel_df$edge_table <- v$graphnel_df$edge_table[which(v$graphnel_df$edge_table$existence != "removed"),]
      v$graphnel_df$node_table <- v$graphnel_df$node_table[which(v$graphnel_df$node_table$existence != "removed"),]
      
      v$graphnel_df$node_table$sink <- FALSE
      v$graphnel_df$node_table[unique(v$graphnel_df$edge_table$to[which(!(v$graphnel_df$edge_table$to %in% v$graphnel_df$edge_table$from))]), "sink"] <- TRUE
      v$graphnel_df$node_table$color.background <- ifelse(v$graphnel_df$node_table$sink, "#0099cc", "#BFFFBF")
      
      
      v$pathway[c("sink.nodes", "order", "graph")] <- df_to_graphnel(v$graphnel_df$node_table, v$graphnel_df$edge_table)[c("sink.nodes", "order", "graph")] # v$pathway <- df_to_graphnel(v$graphnel_df$node_table, v$graphnel_df$edge_table)
      
      v$pathway$graphnel_df <- v$graphnel_df
      
      v$pathway$removed_nodes <- rbind(v$pathway$removed_nodes, removed_nodes)
      v$pathway$removed_edges <- rbind(v$pathway$removed_edges, removed_edges)
      
      v$selected_node_name <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("node_id", "component_id_s", "label"))
      
      v$image_file <- plot_kegg_pathway(graphnel_df = v$graphnel_df, group_graphics = v$pathway$group_nodes, pathway_image = v$pathway_image,
                                        node_colors = v$node_colors$node_colors, edge_mapping = FALSE, highlight_nodes = v$selected_node_name$node_id, highlight_color = "red",
                                        col_legend_title = v$node_colors$col_legend_title, color_bar_lims = v$node_colors$color_bar_lims, present_node_modifications = input$show_modified_nodes, removed_nodes = v$pathway$removed_nodes) %>%
        image_write(tempfile(fileext='png'), format = 'png')
      
      
      saveRDS(v$pathway, file = paste0("collection_dir/", v$collection_name, "/pathways/", v$pathway_name, ".RDS"))
    }
  }, suspended = F)
  
  #### download handler for edited all pathways ####
  output$download_collection <- downloadHandler(
    filename = function(){
      paste0(v$collection_name, ".zip")
    },
    content = function(file) {
      work_dir <- getwd()
      setwd("collection_dir/")
      zip(zipfile = file, files = v$collection_name)
      setwd(work_dir)
    }
  )
  
  ##########################################
  ######## KEGG visualization ##############
  ##########################################
  
  #### expression upload ####
  observeEvent(input$file, {
    
    if(isTRUE(tryCatch( { as.matrix(suppressWarnings(fread(input$file$datapath, stringsAsFactors = F))) }
                 , error = function(e) {TRUE}))) {
        
      output$fc_table_load_error <- renderText({paste("<font color=\"#ff0000\"><b>", "Something is wrong with a file", "</b></font>")})
      hide('vis_buttons')
      hide('pi_and_topology_analysis')
      v$allow_graph_param_update <- FALSE
      
    } else {
      output$fc_table_load_error <- NULL
      v$exp_uploaded <- as.matrix(suppressWarnings(fread(input$file$datapath, stringsAsFactors = F)))
      
      rownames(v$exp_uploaded) <- v$exp_uploaded[,1]
      
      v$entrez_fc <- v$exp_uploaded[,-1, drop = F]
      
      if(isTRUE(tryCatch( { run_psf(entrez.fc = v$entrez_fc, kegg.collection = setNames(object = list(v$pathway), nm = v$pathway_name), calculate.significance = F, ncores = 4) }
                          , error = function(e1) {TRUE}))) {
        
        output$fc_table_load_error <- renderText({paste("<font color=\"#ff0000\"><b>", "Something is wrong with a file", "</b></font>")})
        hide('vis_buttons')
        hide('pi_and_topology_analysis')
        v$allow_graph_param_update <- FALSE
        
      } else {
        output$fc_table_load_error <- NULL
        
        v$pathway_psf <- run_psf(entrez.fc = v$entrez_fc, kegg.collection = setNames(object = list(v$pathway), nm = v$pathway_name), calculate.significance = F, ncores = 4)
        
        show('sample_selection_type_input')
        
        v$node_colors <- node_color_generator(pathway = v$pathway_psf[[v$pathway_name]], log_norm = input$log_norm_checkbox, sample_id = "mean", color_nodes = c("psf_activities", "fold_change")[as.integer(input$node_coloring_type)])
        
        v$graphnel_df$node_table[rownames(v$node_colors$node_values), "expression"] <- v$node_colors$node_values[,1]
        v$graphnel_df$node_table[rownames(v$node_colors$node_values), "signal"] <- v$node_colors$node_values[,2]
        
        v$image_file <- plot_kegg_pathway(graphnel_df = v$graphnel_df, group_graphics = v$pathway$group_nodes, pathway_image = v$pathway_image,
                                          node_colors = v$node_colors$node_colors, edge_mapping = FALSE, edge_in_mode = FALSE, highlight_nodes = NULL, highlight_color = "red",
                                          col_legend_title = v$node_colors$col_legend_title, color_bar_lims = v$node_colors$color_bar_lims, present_node_modifications = input$show_modified_nodes, removed_nodes = v$pathway$removed_nodes) %>%
          image_write(tempfile(fileext='png'), format = 'png')
        
        v$allow_graph_param_update <- TRUE
        
        updateSelectizeInput(session, "node_influence", 
                             choices = v$graphical_data$node_coords$gr_name[which(v$graphical_data$node_coords$node_class == "gene")], 
                             selected = v$graphical_data$node_coords$gr_name[which(v$graphical_data$node_coords$node_class == "gene")][1], server = TRUE)
        
        show('vis_buttons')
        show('pi_and_topology_analysis')
        
      }
      
    }
    
  })
  
  #### psf calculation and sample data mapping ####
  observeEvent(input$sample_selection_type, {
    
    if(input$sample_selection_type == 1) {
      
      if(v$allow_graph_param_update) {
        v$node_colors <- node_color_generator(pathway = v$pathway_psf[[v$pathway_name]], log_norm = input$log_norm_checkbox, sample_id = "mean", color_nodes = c("psf_activities", "fold_change")[as.integer(input$node_coloring_type)])
        
        v$graphnel_df$node_table[rownames(v$node_colors$node_values), "expression"] <- v$node_colors$node_values[,1]
        v$graphnel_df$node_table[rownames(v$node_colors$node_values), "signal"] <- v$node_colors$node_values[,2]
        
        v$image_file <- plot_kegg_pathway(graphnel_df = v$graphnel_df, group_graphics = v$pathway$group_nodes, pathway_image = v$pathway_image,
                                          node_colors = v$node_colors$node_colors, edge_mapping = FALSE, edge_in_mode = FALSE, highlight_nodes = NULL, highlight_color = "red",
                                          col_legend_title = v$node_colors$col_legend_title, color_bar_lims = v$node_colors$color_bar_lims, present_node_modifications = input$show_modified_nodes, removed_nodes = v$pathway$removed_nodes) %>%
          image_write(tempfile(fileext='png'), format = 'png')
        
        hide('sample_selection_input')
      }
      
    } else {
      if(input$sample_selection_type == 2) {
        updateSelectizeInput(session, "selected_sample", label = "Select sample", choices = colnames(v$entrez_fc), selected = colnames(v$entrez_fc)[1], server = T)
        show('sample_selection_input')
      }
      ### will work on this later
      if(input$sample_selection_type == 3) {
        updateSelectizeInput(session, "selected_sample", label = "Select group", choices = "", selected = NULL, server = T)
        show('sample_selection_input')
      }
    }
    
  })
  
  
  observeEvent(c(input$selected_sample, input$node_coloring_type, input$log_norm_checkbox), {
    if(v$allow_graph_param_update) {
      if(input$selected_sample %in% colnames(v$entrez_fc)) {
        v$node_colors <- node_color_generator(pathway = v$pathway_psf[[v$pathway_name]], log_norm = input$log_norm_checkbox, sample_id = input$selected_sample, color_nodes = c("psf_activities", "fold_change")[as.integer(input$node_coloring_type)])
      } else {
        v$node_colors <- node_color_generator(pathway = v$pathway_psf[[v$pathway_name]], log_norm = input$log_norm_checkbox, sample_id = "mean", color_nodes = c("psf_activities", "fold_change")[as.integer(input$node_coloring_type)])
      }  
      
      v$graphnel_df$node_table[rownames(v$node_colors$node_values), "expression"] <- v$node_colors$node_values[,1]
      v$graphnel_df$node_table[rownames(v$node_colors$node_values), "signal"] <- v$node_colors$node_values[,2]
        
      v$image_file <- plot_kegg_pathway(graphnel_df = v$graphnel_df, group_graphics = v$pathway$group_nodes, pathway_image = v$pathway_image,
                                        node_colors = v$node_colors$node_colors, edge_mapping = FALSE, edge_in_mode = FALSE, highlight_nodes = NULL, highlight_color = "red",
                                        col_legend_title = v$node_colors$col_legend_title, color_bar_lims = v$node_colors$color_bar_lims, present_node_modifications = input$show_modified_nodes, removed_nodes = v$pathway$removed_nodes) %>%
        image_write(tempfile(fileext='png'), format = 'png')
    }
    
  })
  
  #### node expression change observer ####
  observeEvent(input$change_exp, {
    print(v$selected_node_name[,"node_id"])
    updateNumericInput(session, "exp_slider", label = paste0(v$selected_node_name[,"label"], " FC value"),
                       value = v$pathway_psf[[v$pathway_name]]$exp_fc[v$selected_node_name[,"node_id"], input$selected_sample], 
                       min = 0, max = 10, step = 0.001)
    
    show('exp_change_panel')
    
  })
  
  observeEvent(input$close_slider, {
    v$searching_state <- FALSE
    hide('exp_change_panel')
  })
  
  
  #### PSF updater ####
  observeEvent(input$change_exp_submit, {
    
    hide('exp_change_panel')
    
    ### cleaning old changes
    # v$psf_and_colors <- psf_signal_calculator_and_coloring(entrez_fc = v$entrez_fc, pathway = v$pathway,
    #                                                        pathway_name = v$pathway_name, update_mod = FALSE)
    
    ### editing new changed FC value
    graph::nodeData(v$psf_and_colors$psf_graph[[1]][[v$pathway_name]]$graph, v$selected_node_name[,"node_id"],attr = "expression") <- input$exp_slider
    
    v$psf_and_colors <- psf_signal_calculator_and_coloring(entrez_fc = v$entrez_fc, pathway = v$pathway,
                                                           pathway_name = v$pathway_name, update_mod = TRUE, old_graph = v$psf_and_colors$psf_graph)
    
    # v$pathway_exp_colored = TRUE
    # 
    # v$pathway_psf_colored = TRUE
    
    v$pathway_psf_colored = TRUE
    
    v$color_bar_lims <- range(v$psf_and_colors$signal_values)
    
    v$mapping_value_type <- "PSF"
    
    v$color_bar_psf_mode <- TRUE
    
    v$selected_node_name <- NULL
    
    v$image_file <- kegg_node_mapper(group_graphics = v$pathway_data$group_graphics, kegg_pathway_graphics = v$graphical_data, pathway_name = v$pathway_name, pathway_image = v$pathway_image, highlight.genes = v$selected_node_name, color.genes = v$psf_and_colors$psf_colors, color_bar_psf_mode = v$color_bar_psf_mode, col_legend_title = v$col_legend_title, color_bar_lims = v$color_bar_lims, draw_color_bar = v$draw_color_bar) %>% 
      image_write(tempfile(fileext='png'), format = 'png')
    
    v$visnet_list <- visnet_creator(v$graphical_data, node_colors = v$psf_and_colors$psf_colors, col_legend_title = v$col_legend_title, color_bar_lims = v$color_bar_lims)
    
    ### temporarely commented ggplot of sink values
    # sink_gg_plot <- ggplot(v$psf_and_colors$sink_signals, aes(x=sink_name, y=signal, fill=signal)) +
    #   geom_boxplot(color="black", lwd=0.2, outlier.shape=NA) +
    #   geom_point(color = v$psf_and_colors$sink_signals$dot_color) +
    #   geom_jitter(color = v$psf_and_colors$sink_signals$dot_color, width = 0.2) +
    #   guides(fill=FALSE) +
    #   coord_flip() +
    #   theme_bw()
    # 
    # v$sink_values_plot <- ggplotly(sink_gg_plot, width = 400, height = magick::image_info(v$pathway_image)$height)
    
  })
  
  observeEvent(input$change_edge_weight, {
    
    selected_edge <- as.data.table(v$graphical_data$edge_coords)[from == v$selected_node_name[1,"node_id"] & to == v$selected_node_name[2,"node_id"] | from == v$selected_node_name[2,"node_id"] & to == v$selected_node_name[1,"node_id"]]
    
    selected_interaction <- rbind(v$graphical_data$node_coords[which(v$graphical_data$node_coords$node_id == selected_edge$from), c("gr_name", "node_id")],
                                  v$graphical_data$node_coords[which(v$graphical_data$node_coords$node_id == selected_edge$to), c("gr_name", "node_id")])
    
    impact <- unname(unlist(graph::edgeData(v$psf_and_colors$psf_graph[[1]][[v$pathway_name]]$graph, selected_interaction[1,"node_id"], selected_interaction[2,"node_id"], attr = "impact")))
    
    if(impact == 1) {
      edge_type <- "activation"
    } else {
      edge_type <- "inhibition"
    }
    
    updateSelectizeInput(session, "vis_subtype", selected = edge_type)
    
    updateNumericInput(session, "vis_edge_weight", label = paste0(selected_interaction[1,"gr_name"], " to ", selected_interaction[2,"gr_name"], " weight"),
                       value = unlist(unname(graph::edgeData(v$psf_and_colors$psf_graph[[1]][[v$pathway_name]]$graph, selected_interaction[1,"node_id"], selected_interaction[2,"node_id"], attr = "weight"))), 
                       min = 0, max = 1, step = 0.001)
    
    show('edge_attr_editing_panel_vis')
    
  })
  
  observeEvent(input$close_edge_attr_vis_panel, {
    hide('edge_attr_editing_panel_vis')
  })
  
  observeEvent(input$edit_vis_edge, {
    
    selected_edge <- as.data.table(v$graphical_data$edge_coords)[from == v$selected_node_name[1,"node_id"] & to == v$selected_node_name[2,"node_id"] | from == v$selected_node_name[2,"node_id"] & to == v$selected_node_name[1,"node_id"]]
    
    selected_interaction <- rbind(v$graphical_data$node_coords[which(v$graphical_data$node_coords$node_id == selected_edge$from), c("gr_name", "node_id")],
                                  v$graphical_data$node_coords[which(v$graphical_data$node_coords$node_id == selected_edge$to), c("gr_name", "node_id")])
    
    graph::edgeData(v$psf_and_colors$psf_graph[[1]][[v$pathway_name]]$graph, selected_interaction[1,"node_id"], selected_interaction[2,"node_id"], attr = "weight") <- input$vis_edge_weight
    
    graph::edgeData(v$psf_and_colors$psf_graph[[1]][[v$pathway_name]]$graph, selected_interaction[1,"node_id"], selected_interaction[2,"node_id"], attr = "subtype2") <- input$vis_subtype
    
    v$psf_and_colors <- psf_signal_calculator_and_coloring(entrez_fc = v$entrez_fc, pathway = v$pathway,
                                                           pathway_name = v$pathway_name, update_mod = TRUE, old_graph = v$psf_and_colors$psf_graph)
    
    
    v$pathway_psf_colored = TRUE
    
    v$color_bar_lims <- range(v$psf_and_colors$signal_values)
    
    v$mapping_value_type <- "PSF"
    
    v$color_bar_psf_mode <- TRUE
    
    v$selected_node_name <- NULL
    
    v$image_file <- kegg_node_mapper(group_graphics = v$pathway_data$group_graphics, kegg_pathway_graphics = v$graphical_data, pathway_name = v$pathway_name, pathway_image = v$pathway_image, highlight.genes = v$selected_node_name, color.genes = v$psf_and_colors$psf_colors, color_bar_psf_mode = v$color_bar_psf_mode, col_legend_title = v$col_legend_title, color_bar_lims = v$color_bar_lims, draw_color_bar = v$draw_color_bar) %>% 
      image_write(tempfile(fileext='png'), format = 'png')
    
    v$visnet_list <- visnet_creator(v$graphical_data, node_colors = v$psf_and_colors$psf_colors, col_legend_title = v$col_legend_title, color_bar_lims = v$color_bar_lims)
    
    hide('edge_attr_editing_panel_vis')
    
  })
  
  #### reset PSF calculations with default expression values ####
  # observeEvent(input$reset_psf, {
  #   v$psf_and_colors <- psf_signal_calculator_and_coloring(entrez_fc = v$entrez_fc, pathway = v$pathway,
  #                                                          pathway_name = v$pathway_name, update_mod = FALSE)
  #   
  #   v$pathway_psf_colored = TRUE
  #   
  #   v$color_bar_lims <- range(v$psf_and_colors$signal_values)
  #   
  #   v$mapping_value_type <- "PSF"
  #   
  #   v$color_bar_psf_mode <- TRUE
  #   
  #   v$image_file <- kegg_node_mapper(group_graphics = group_graphics[[v$pathway_name]], kegg_pathway_graphics = v$graphical_data, pathway_name = v$pathway_name, pathway_image = v$pathway_image, highlight.genes = v$selected_node_name, color.genes = v$psf_and_colors$psf_colors, color_bar_psf_mode = v$color_bar_psf_mode, col_legend_title = v$col_legend_title, color_bar_lims = v$color_bar_lims, draw_color_bar = v$draw_color_bar) %>% 
  #     image_write(tempfile(fileext='png'), format = 'png')
  #   
  #   v$visnet_list <- visnet_creator(v$graphical_data, node_colors = v$psf_and_colors$psf_colors)
  #   
  # })
  
  
  ########################################
  ##### image and network rendering ######
  ########################################
  
  #### pathway image rendering #### 
  output$pathway_image <- renderImage({
    return(list(src = ifelse(is.null(v$graphnel_df), "", v$image_file), contentType = "image/png", alt = "pathway_image"))
  }, deleteFile = TRUE)
  
  #### sink data plot rendering ####
  # output$sink_plot <- renderPlotly({
  #   v$sink_values_plot
  # })

  
  #### visnet rendering ####
  output$visnet <- renderVisNetwork({
    
    if(!is.null(v$graphnel_df)) {
      visNetwork(nodes = vis_extract(v$graphnel_df, node_colors = v$node_colors)$node_table, edges = vis_extract(v$graphnel_df)$edge_table, width = "100%", height = "800px") %>% 
        visIgraphLayout(layout = "layout_nicely") %>% 
        visInteraction(navigationButtons = TRUE, multiselect = T, keyboard = F, selectConnectedEdges = T) %>%
        visOptions(manipulation = list(enabled = TRUE, 
                                       addNode = htmlwidgets::JS("function(data, callback) {
                                                                callback(data);
                                                                console.info('add_node')
                                                                Shiny.onInputChange('add_node', {node_data : data});
                                                                ;}"),
                                       editNode = htmlwidgets::JS("function(data, callback) {
                                                                callback(data);
                                                                console.info('edit_node')
                                                                Shiny.onInputChange('edit_node', {edit_node_data : data});
                                                                ;}"),
                                       # editEdge = htmlwidgets::JS("function(data, callback) {
                                       #                            callback(data);
                                       #                            console.info('vis_edit_edge')
                                       #                            Shiny.onInputChange('vis_edit_edge', {edit_edge_data : data});
                                       #                            ;}"),
                                       editEdge = TRUE,
                                       deleteNode = TRUE, initiallyActive = TRUE
        ), highlightNearest = T) %>%
        visExport() %>%
        visEvents(click = "function(nodes) {
                      console.info('click')
                      console.info(nodes)
                      Shiny.onInputChange('clicked_node', {nodes : nodes.nodes, edges : nodes.edges});
                      ;}"
                  # , beforeDrawing = "function(ctx) {
                  #       const image = document.getElementById('scream');
                  #       const canvasWidth = ctx.canvas.width;
                  #       const canvasHeight = ctx.canvas.height;
                  #       
                  #       var x = -image.width/2;
                  #       var y = -image.height/2;
                  # 
                  #       ctx.drawImage(image, -1014, -1014);}"
        )
    }
    
  })
  
  #### visnet editor ###
  observeEvent(input$add_node, {

    updateTextInput(session, inputId = "node_label", value = "")
    updateTextInput(session, inputId = "node_genes", value = "")
    updateSelectizeInput(session, inputId = "node_type", selected = "")
    updateSelectizeInput(session, inputId = "node_function", selected = "mean")
    updateNumericInput(session, inputId = "node_fc", value = 1)
    
    show('node_attrs')

  })
  
  observeEvent(input$edit_node, {
    
    v$node_edit_mode <- TRUE
    updateTextInput(session, inputId = "node_label", value = v$graphnel_df$node_table[unlist(input$clicked_node$nodes), "label"])
    updateTextInput(session, inputId = "node_genes", value = v$graphnel_df$node_table[unlist(input$clicked_node$nodes), "component_id_s"])
    updateSelectizeInput(session, inputId = "node_type", selected = v$graphnel_df$node_table[unlist(input$clicked_node$nodes), "type"])
    updateSelectizeInput(session, inputId = "node_function", selected = v$graphnel_df$node_table[unlist(input$clicked_node$nodes), "psf_function"])
    updateNumericInput(session, inputId = "node_fc", value = v$graphnel_df$node_table[unlist(input$clicked_node$nodes), "expression"])
    updateNumericInput(session, inputId = "node_width", value = v$graphnel_df$node_table[unlist(input$clicked_node$nodes), "node_width"])
    updateNumericInput(session, inputId = "node_height", value = v$graphnel_df$node_table[unlist(input$clicked_node$nodes), "node_height"])
    
    show('node_attrs')
    
  })
  
  
  observeEvent(input$node_type, {
    node_default_dims <- list(gene = c(46, 17), compound = c(8, 8), map = c(150, 35), linker = c(46, 17), event = c(150, 35))
    if(input$node_type != "") {
      updateNumericInput(session, inputId = "node_width", value = node_default_dims[[input$node_type]][1])
      updateNumericInput(session, inputId = "node_height", value = node_default_dims[[input$node_type]][2])
    }
  })
  
  
  observe({
    invalidateLater(3000)
    visNetworkProxy("visnet") %>%
      visGetPositions()
    if(!is.null(input$visnet_positions)) {
      v$visnet_coord_ranges <- list(x = range(as.data.frame(do.call(rbind, input$visnet_positions))[,"x"]),
                                    y = range(as.data.frame(do.call(rbind, input$visnet_positions))[,"y"]))
      v$visnet_coords <- list(x = as.data.frame(do.call(rbind, input$visnet_positions))[,"x"],
                              y = as.data.frame(do.call(rbind, input$visnet_positions))[,"y"])
    }
  })
  
  observeEvent(input$close_node_attr_panel, {
    hide('node_attrs')
  })
  
  output$node_data_adding_error <- renderText({
    paste("<font color=\"#ff0000\"><b>", v$node_data_adding_error_visnet, "</b></font>")
  })
  
  observeEvent(input$submit_node_attrs, {
    
    if(input$node_label == "" | input$node_genes == "" | input$node_type == "") {
      v$node_data_adding_error_visnet <- "Fill all necessary fields"
      delay(3000, v$node_data_adding_error_visnet <- "")
    } else {
      if(v$node_edit_mode) {
        v$graphnel_df$node_table[unlist(input$clicked_node$nodes), "label"] <- input$node_label
        v$graphnel_df$node_table[unlist(input$clicked_node$nodes), "component_id_s"] <- input$node_genes
        v$graphnel_df$node_table[unlist(input$clicked_node$nodes), "type"] <- input$node_type
        v$graphnel_df$node_table[unlist(input$clicked_node$nodes), "psf_function"] <- input$node_function
        v$graphnel_df$node_table[unlist(input$clicked_node$nodes), "expression"] <- input$node_fc
        v$graphnel_df$node_table[unlist(input$clicked_node$nodes), "node_width"] <- input$node_width
        v$graphnel_df$node_table[unlist(input$clicked_node$nodes), "node_height"] <- input$node_height
        
        v$pathway_change_alert <- v$pathway_change_alert + 1
        v$node_edit_mode <- FALSE
        hide('node_attrs')
      } else {
        ### when visnet editor is used update visnet data frame directly without rebuilding node and edge tables
        node_ids <- as.integer(v$graphnel_df$node_table$node_id)
        
        ### continue here
        x2 = (input$add_node$node_data$x - v$visnet_coord_ranges$x[1]) * (max(v$graphnel_df$node_table$x) - min(v$graphnel_df$node_table$x)) / (v$visnet_coord_ranges$x[2] - v$visnet_coord_ranges$x[1]) + min(v$graphnel_df$node_table$x)
        
        y2 = (input$add_node$node_data$y - v$visnet_coord_ranges$y[1]) * (max(v$graphnel_df$node_table$y) - min(v$graphnel_df$node_table$y)) / (v$visnet_coord_ranges$y[2] - v$visnet_coord_ranges$y[1]) + min(v$graphnel_df$node_table$y)
        
        new_node <- data.frame(
          node_id = max(node_ids) + 1,
          label = input$node_label,
          component_id_s = input$node_genes,
          type = input$node_type,
          psf_function = input$node_function,
          expression = input$node_fc,
          signal = 1,
          sink = F,
          x_start = x2 - input$node_width*0.5,
          y_start = y2 - input$node_height*0.5,
          x_end  = x2 + input$node_width*0.5,
          y_end =  y2 + input$node_height*0.5,
          node_width = input$node_width,
          node_height = input$node_height,
          title = paste0(input$node_genes, "(", input$node_label, ")", "<br>", "Node id ", max(node_ids) + 1),
          shape = ifelse(input$node_type == "compound", "dot", "box"),
          color.background = "#BFFFBF",
          color.border = "#BFFFBF",
          borderWidth = 2,
          font.size = 22,
          size = ifelse(input$node_type == "compound", 10, 25),
          font.color = "#000000",
          vis_width = ifelse(input$node_type %in% c("map", "event"), input$node_width*2, NA),
          x = x2, # input$add_node$node_data$x,
          y = y2, # input$add_node$node_data$y,
          existence = "added",
          change_info = "no_change",
          data_source = "",
          row.names = max(node_ids) + 1,
          stringsAsFactors = F
        )
        
        v$graphnel_df$node_table <- rbind(v$graphnel_df$node_table,
                                          new_node
        )
        
        v$pathway_change_alert <- v$pathway_change_alert + 1
        
        # visNetworkProxy("visnet") %>%
        #   visUpdateNodes(vis_extract(v$graphnel_df)$node_table)
        
        hide('node_attrs')
      }
    }
  })
  
  
  observeEvent(input$visnet_graphChange, {
    
    if(input$visnet_graphChange$cmd == "addEdge") {
      
      updateSelectizeInput(session, inputId = "edge_type", selected = "PPrel")
      updateSelectizeInput(session, inputId = "subtype1_visnet", selected = "")
      updateSelectizeInput(session, inputId = "subtype2_visnet", selected = "")
      updateSelectizeInput(session, inputId = "interaction_source_visnetwork", selected = "image") 
      updateNumericInput(session, inputId = "edge_wight", value = 1)
      
      show('edge_attrs')
      
    }
    ### will work on it later
    if(input$visnet_graphChange$cmd == "editEdge") {
      
      updateSelectizeInput(session, inputId = "edge_type", selected = v$graphnel_df$edge_table[input$visnet_graphChange$id, "type"])
      updateSelectizeInput(session, inputId = "subtype1_visnet", selected = v$graphnel_df$edge_table[input$visnet_graphChange$id, "subtype1"])
      updateSelectizeInput(session, inputId = "subtype2_visnet", selected = v$graphnel_df$edge_table[input$visnet_graphChange$id, "subtype2"])
      updateSelectizeInput(session, inputId = "interaction_source_visnetwork", selected = v$graphnel_df$edge_table[input$visnet_graphChange$id, "data_source"])
      updateNumericInput(session, inputId = "edge_wight", value = v$graphnel_df$edge_table[input$visnet_graphChange$id, "weight"])
      updateActionButton(session, inputId = "submit_edge_attrs", label = "Edit edge")

      show('edge_attrs')
    }
    
    if(input$visnet_graphChange$cmd == "deleteElements") {
      
      if(length(unlist(input$visnet_graphChange$nodes)) > 0){
        v$graphnel_df$node_table[unlist(input$visnet_graphChange$nodes), "existence"] <- "removed"
      }
      if(length(unlist(input$visnet_graphChange$edges)) > 0) {
        v$graphnel_df$edge_table[unlist(input$visnet_graphChange$edges), "existence"] <- "removed"
      }
      v$pathway_change_alert <- v$pathway_change_alert + 1
    }
  })
  
  output$edge_data_adding_error <- renderText({
    paste("<font color=\"#ff0000\"><b>", v$edge_adding_error_visnet, "</b></font>")
  })
  
  #### visnet edit/add edge ####
  observeEvent(input$submit_edge_attrs, {
    if(input$subtype1_visnet == "" | is.null(input$interaction_source_visnetwork)) {
      v$edge_adding_error_visnet <- "Fill all necessary fields"
      delay(3000, v$edge_adding_error_visnet <- "")
    } else {
      new_edge_id <- paste0(input$visnet_graphChange$from, "|", input$visnet_graphChange$to)
      
      data_source <- input$interaction_source_visnetwork
      # data_source[which(data_source == "database")] <- v$selected_interaction$source
      data_source <- paste(data_source, collapse = ", ")
      
      if(input$visnet_graphChange$cmd == "editEdge") {
        v$graphnel_df$edge_table[input$visnet_graphChange$id, "from"] <- input$visnet_graphChange$from
        v$graphnel_df$edge_table[input$visnet_graphChange$id, "to"] <- input$visnet_graphChange$to
        v$graphnel_df$edge_table[input$visnet_graphChange$id, "subtype1"] <-  input$subtype1_visnet
        v$graphnel_df$edge_table[input$visnet_graphChange$id, "subtype2"] <-  input$subtype2_visnet
        v$graphnel_df$edge_table[input$visnet_graphChange$id, "color"] <- line_col[input$subtype1_visnet]
        v$graphnel_df$edge_table[input$visnet_graphChange$id, "arrows.to.type"] <- kegg_arrows_type[input$subtype1_visnet]
        v$graphnel_df$edge_table[input$visnet_graphChange$id, "dashes"] <- input$subtype2_visnet == "indirect effect" | input$subtype1_visnet == "indirect effect"
        v$graphnel_df$edge_table[input$visnet_graphChange$id, "change_info"] <-  "attr_change"
        v$graphnel_df$edge_table[input$visnet_graphChange$id, "data_source"] <-  paste(input$interaction_source, collapse = ", ")
        v$graphnel_df$edge_table[input$visnet_graphChange$id, "weight"] <- input$edge_wight
      } else {
        if(!(new_edge_id %in% rownames(v$graphnel_df$edge_table))) {
          
          new_edge <- data.frame(
            id = new_edge_id,
            from = input$visnet_graphChange$from,
            to = input$visnet_graphChange$to,
            color = line_col[input$subtype1_visnet],
            arrows.to.enabled = TRUE,
            arrows.to.type = kegg_arrows_type[input$subtype1_visnet],
            label = "",
            dashes = input$subtype1_visnet == "indirect effect" | input$subtype2_visnet == "indirect effect",
            type = "",
            subtype1 = input$subtype1_visnet,
            subtype2 = input$subtype2_visnet,
            state = "",
            weight = input$edge_wight,
            existence = "exist",
            change_info = "added",
            data_source = data_source,
            row.names = new_edge_id
          )
          
          v$graphnel_df$edge_table <- rbind(v$graphnel_df$edge_table, new_edge)
        }
      }
      
      v$pathway_change_alert <- v$pathway_change_alert + 1
      hide("edge_attrs")
      updateActionButton(session, inputId = "submit_edge_attrs", label = "Save")
    }
  })
  
  observeEvent(input$close_edge_attr_panel, {
    hide("edge_attrs")
  })
  
  
  #### undernetwork tables ####
  observe({
    output$network_node_table <- renderDataTable({
      datatable(
        v$graphnel_df$node_table,          
        escape=FALSE,
        # filter = 'top',
        # class = "cell-border",
        editable = TRUE,
        class   = 'cell-border compact hover',
        selection = list(mode = "none", target = 'row'),
        options = list(searchHighlight = TRUE,
                       order = list(0, 'asc'),
                       pageLength = 100,
                       autoWidth  = T,
                       columnDefs = list(list(
                         targets  = c(c(1:9), 15),
                         render   = JS(
                           "function(data, type, row, meta) {",
                           "return type === 'display' && data.length > 10 ?",                                    "'<span title=\"' + data + '\">' +
                                data.substr(0, 10) + '...</span>' : data;", "}"))),
                       initComplete = JS("function(settings, json) {",
                                         "$(this.api().table().header()).css({'background-color': '#3474B7', 'color': '#fff'});",
                                         "}"),
                       scrollY = 200,scrollX = TRUE
        ),
        style = 'bootstrap'
      )
    })
    
    network_node_proxy <- dataTableProxy('network_node_table', deferUntilFlush = FALSE)
    
    # observeEvent(input$network_node_table_cell_edit, {
    #   if(!is.null(v$node_network_table_filtered)) {
    #     node_ind <- which(v$node_table$name == v$node_network_table_filtered$name[input$network_node_table_cell_edit$row])
    #     new_network_node_row <- v$node_table[node_ind,]
    #   } else {
    #     node_ind <- which(v$node_table$name == v$node_network_table$name[input$network_node_table_cell_edit$row])
    #     new_network_node_row <- v$node_table[node_ind,]
    #   }
    #   
    #   new_network_node_row[10] <- input$network_node_table_cell_edit$value
    #   new_network_node_row[8] <- USER$User
    #   new_network_node_row <- as.character(new_network_node_row[-2])
    #   
    #   if(input$network_node_table_cell_edit$col == 10) {
    #     # new_network_node_row[8] <- as.character(Sys.time())
    #     node_curating_function(v$node_table, v$edge_table, new_network_node_row, TRUE, node_ind)
    #   }
    #   
    # })
    
    output$network_edge_table <- renderDataTable({
      datatable(
        v$graphnel_df$edge_table, 
        escape=FALSE,
        # filter = 'top',
        class = "cell-border",
        selection = list(mode = "none", target = 'row'),
        editable = TRUE,
        options = list(searchHighlight = TRUE,
                       # order = list(13, 'desc'),
                       pageLength = 100,
                       autoWidth  = T,
                       # columnDefs = list(list(
                       #   targets  = c(1:6),
                       #   render   = JS(
                       #     "function(data, type, row, meta) {",
                       #     "return type === 'display' && data.length > 10 ?",                                    "'<span title=\"' + data + '\">' +
                       #          data.substr(0, 10) + '...</span>' : data;", "}"))),
                       initComplete = JS("function(settings, json) {",
                                         "$(this.api().table().header()).css({'background-color': '#3474B7', 'color': '#fff'});",
                                         "}"),
                       scrollY = 200,scrollX = TRUE
        ),
        style = 'bootstrap'
      )
    })
    
    network_edge_proxy <- dataTableProxy('network_edge_table', deferUntilFlush = FALSE)
    
    # observeEvent(input$network_edge_table_cell_edit, {
    #   if(!is.null(v$edge_network_table_filtered)) {
    #     edge_ind <- which(v$edge_table$source == v$edge_network_table_filtered$source[input$network_edge_table_cell_edit$row] & v$edge_table$target == v$edge_network_table_filtered$target[input$network_edge_table_cell_edit$row])
    #     new_network_edge_row <- v$edge_table[edge_ind,]
    #   } else {
    #     edge_ind <- which(v$edge_table$source == v$edge_network_table$source[input$network_edge_table_cell_edit$row] & v$edge_table$target == v$edge_network_table$target[input$network_edge_table_cell_edit$row])
    #     new_network_edge_row <- v$edge_table[edge_ind,]
    #   }
    #   
    #   new_network_edge_row[14] <- input$network_edge_table_cell_edit$value
    #   new_network_edge_row[12] <- USER$User
    #   new_network_edge_row <- as.character(new_network_edge_row)
    #   
    #   if(input$network_edge_table_cell_edit$col == 14) {
    #     # new_network_edge_row[13] <- as.character(Sys.time())
    #     # edge_curating_function(v$edge_table, v$node_table, edge_row, edit_mode, row_ind, user)
    #     edge_curating_function(v$edge_table, v$node_table, new_network_edge_row, TRUE, edge_ind)
    #   }
    #   
    # })
    
    observe({
      if(!is.null(unlist(input$clicked_node$edges))) {
        # visNetworkProxy("visnet") %>%
        #   visGetEdges(input = "edge_list")
        
        # edge_matrix <- unlist(sapply(unlist(input$clicked_node$edges), function(x) {
        #   c(input$edge_list[[x]]$from, input$edge_list[[x]]$to)
        # }))
        
        # edge_matrix[1,] <- as.character(v$visnet_list$edges[edge_matrix[1,],1])
        # edge_matrix[2,] <- as.character(v$visnet_list$edges[edge_matrix[2,],2])
        
        
        v$edge_network_table_filtered <- v$graphnel_df$edge_table[unlist(input$clicked_node$edges),]
        replaceData(network_edge_proxy, v$edge_network_table_filtered, resetPaging = FALSE)
        
        # if(!is.null(edge_matrix)) {
        #   v$edge_network_table_filtered <- v$edge_network_table[which(v$edge_network_table$from %in% edge_matrix[1,] & v$edge_network_table$to %in% edge_matrix[2,]),]
        #   replaceData(network_edge_proxy, v$edge_network_table_filtered, resetPaging = FALSE)
        # }
      } else {
        replaceData(network_edge_proxy, v$graphnel_df$edge_table, resetPaging = FALSE)
        v$edge_network_table_filtered <- NULL
      } 
    })
    
    observe({
      if(!is.null(unlist(input$clicked_node$nodes))) {
        v$node_network_table_filtered <- v$graphnel_df$node_table[unlist(input$clicked_node$nodes),] # v$node_network_table[which(v$node_network_table$node_id %in% selected_nods),]
        replaceData(network_node_proxy, v$node_network_table_filtered, resetPaging = FALSE)
      } else {
        replaceData(network_node_proxy, v$graphnel_df$node_table, resetPaging = FALSE)
        v$node_network_table_filtered <- NULL
      }
    })
    
  })
  
  
  #### partial influence node detection by network node click ####
  observeEvent(input$get_partial_influence, {
    output$partial_inf_err <- NULL
    if(!is.null(unlist(input$clicked_node$nodes)) & input$influence_type != "None") {
      
      if(input$sample_selection_type == 1) {
        sample_id = "all"
      } 
      if(input$sample_selection_type == 2) {
        sample_id = input$selected_sample
      }
      # for group samples, will finish in the future
      # if(input$sample_selection_type == 3) {
      #   sample_id = input$selected_sample
      # }
      
      influence_nodes <- run_pi(pathway = v$pathway_psf[[v$pathway_name]], unlist(input$clicked_node$nodes), influence_direction = input$influence_type, sample_id = sample_id, node_combinations = 1, get_influence_matrix = FALSE, ncores = 4)
      
      if(length(sample_id) > 1) {
        influence_nodes <- influence_nodes$mean_influence_params$node_id
      } else {
        influence_nodes <- influence_nodes$ordered_influence_nodes[,1]
      }
      
      if(input$influence_node_num == "all") {
        influence_nodes <- influence_nodes
      } else {
        if(length(influence_nodes) > as.numeric(input$influence_node_num)) {
          influence_nodes <- influence_nodes[1:as.numeric(input$influence_node_num)]
        }
      }
      
      vis_node_table <- vis_extract(v$graphnel_df, node_colors = v$node_colors)$node_table
      
      vis_node_table[unlist(input$clicked_node$nodes),"font.color"] <- "#ff0000"
      # new_node_data[names(influence_nodes),"size"] <- 30
      vis_node_table[influence_nodes,"label"] <- paste(1:length(influence_nodes), "\n", vis_node_table[influence_nodes,"label"])
      
      ### update coordinates to be in visnet scale
      vis_node_table[,"x"] <- unlist(v$visnet_coords$x)
      vis_node_table[,"y"] <- unlist(v$visnet_coords$y)
      
      visNetworkProxy("visnet") %>%
        visUpdateNodes(vis_node_table)
      
    } else {
      output$partial_inf_err <- renderText({paste("<font color=\"#ff0000\"><b>", "Select network node and influence type", "</b></font>")})
    }
    
  })
  
  observeEvent(input$show_drug_interactions, {
    
    if(!is.null(unlist(input$clicked_node$nodes))) {
      selected_gene_nodes <- unlist(input$clicked_node$nodes)[which(v$graphnel_df$node_table[unlist(input$clicked_node$nodes), "type"] == "gene")]
      entrez_ids <- unique(sapply(v$graphnel_df$node_table[selected_gene_nodes,"component_id_s"], function(x) {unlist(strsplit(x, split = ","))}))
      v$drug_table <- protein_drug_interactions[entrez_id %in% entrez_ids]
      
    } else {
      entrez_ids <- unique(sapply(v$graphnel_df$node_table[which(v$graphnel_df$node_table$type == "gene"),"component_id_s"], function(x) {unlist(strsplit(x, split = ","))}))
      v$drug_table <- protein_drug_interactions[entrez_id %in% entrez_ids]
    }
    
    show('drug_table')
    
  })
  
  observeEvent(input$close_drag_table, {
    v$searching_state <- FALSE
    hide('drug_table')
  })
  
  #### drug table output ####
  output$drug_data_table <- DT::renderDataTable({datatable(
    v$drug_table[,c("gene_name", "interaction_claim_source", "interaction_types", "drug_name", "drug_concept_id", "interaction_group_score", "PMIDs")],
    filter = 'bottom',
    class   = 'cell-border compact hover',
    fillContainer = T,
    escape = F,
    selection = list(mode = 'single'),
    options = list(dom = 'rtp',
                   pageLength = 100,
                   searchHighlight = TRUE,
                   initComplete = JS("function(settings, json) {",
                                     "$(this.api().table().header()).css({'background-color': '#3474B7', 'color': '#fff'});",
                                     "}"),
                   scrollY = 250,scrollX = TRUE
    ),
    style = 'bootstrap', editable = FALSE)
  })
  
  #### visnet search ####
  observeEvent(input$search_go, {
    id <- v$graphnel_df$node_table$node_id[grep(input$network_search, v$graphnel_df$node_table$title, ignore.case = T)]
    if(length(id) > 0) {
      visNetworkProxy("visnet") %>%
        visFocus(id = id[1], scale = 2)
    }
  })
  
})
