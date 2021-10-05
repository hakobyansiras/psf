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
### library(plotrix)

load("whole_data_unit.RData")

#############################################
######### Pathway curation block ############
#############################################

#### pathway kgml and image downloader ####
kegg_data_downloader <- function(pathway_name) {
  
  pathway_id <- pathway_codes_new[pathway_name]
  
  kgml <- paste0("http://rest.kegg.jp/get/", pathway_id, "/kgml")
  
  image <- paste0("http://rest.kegg.jp/get/", pathway_id, "/image")
  
  kgml_path <- tempfile()
  
  download.file(kgml, destfile = kgml_path, method = "auto")
  
  image_path <- tempfile()
  
  download.file(image, destfile = image_path, method = "auto")
  
  img <- magick::image_read(image_path)
  
  graph <- psf::generate.kegg.collection.from.kgml(kgml_path)[[1]]
  
  graph::edgeDataDefaults(graph$graph, attr = "existence") <- "exist"
  graph::nodeDataDefaults(graph$graph, attr = "existence") <- "exist"
  graph::edgeDataDefaults(graph$graph, attr = "change_info") <- "no_change"
  graph::nodeDataDefaults(graph$graph, attr = "change_info") <- "no_change"
  graph::edgeDataDefaults(graph$graph, attr = "data_source") <- "kegg"
  graph::nodeDataDefaults(graph$graph, attr = "data_source") <- "kegg"  
  
  graph::edgeDataDefaults(graph$graph, attr = "weight") <- 1
  graph::nodeDataDefaults(graph$graph, attr = "psf_function") <- "mean"
  
  return(list(pathway = graph, image = img))
  
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
      rect( z$graphics$x-z$graphics$width*0.5, 
            z$graphics$y+z$graphics$height*0.5, 
            z$graphics$x+z$graphics$width*0.5, 
            z$graphics$y-z$graphics$height*0.5, 
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

#### disconnected node detector ####
disconnected_node_detector <- function (g) {
  edges = graph::edges(g)
  connected.nodes = vector()
  for (snode in names(edges)) {
    if (length(edges[[snode]]) > 0) {
      if (!(snode %in% connected.nodes)) 
        connected.nodes = c(connected.nodes, snode)
      for (tnode in edges[[snode]]) {
        if (!(tnode %in% connected.nodes)) 
          connected.nodes = c(connected.nodes, tnode)
      }
    }
  }
  all.nodes = graph::nodes(g)
  disconnected.nodes = setdiff(all.nodes, connected.nodes)
  return(disconnected.nodes)
}

#### database search function ####
node_pairs_generator <- function(selected_node_name, direction = TRUE, direction_state = TRUE, pathway_realted = FALSE, pathway_data) {
  
  if(pathway_realted) {
    node_a <- unlist(strsplit(selected_node_name[1,"entrez_id"], split = ","))
    
    pathway_nodes <- do.call(rbind, apply(pathway_data, 1, function(x) {
      cbind(unlist(strsplit(x[9], split = ",")), 
            rep(x[8], length(unlist(strsplit(x[9], split = ",")))))
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
    node_a <- unlist(strsplit(selected_node_name[1,"entrez_id"], split = ","))
    node_b <- unlist(strsplit(selected_node_name[2,"entrez_id"], split = ","))
    
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

#### visnet data generator ####
visnet_creator <- function(graphical_data, node_colors = NULL) {
  graphical_data$edge_coords <- graphical_data$edge_coords[which(graphical_data$edge_coords$lty == "solid"),]
  graphical_data$node_coords <- graphical_data$node_coords[which(graphical_data$node_coords$exist),]
  
  node_shapes <- unname(sapply(graphical_data$node_coords$node_name, function(x) {
    if(grepl("cpd",x)) {
      "dot"
    } else {
      "box"
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
    
  } else {
    color <- unname(sapply(graphical_data$node_coords$sink, function(x) {
      if(x) {
        "#0099cc"
      } else {
        "#BFFFBF"
      }
    }))
    
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
                      label = graphical_data$node_coords$gr_name,
                      shape = node_shapes, color = color, 
                      title = graphical_data$node_coords$hover_name,
                      font.size = rep(22, nrow(graphical_data$node_coords)), size = size,
                      font.color = font_color,
                      x = (graphical_data$node_coords$x_start + graphical_data$node_coords$x_end)/2,
                      y = (graphical_data$node_coords$y_start + graphical_data$node_coords$y_end)/2
  )
  
  arrows_type <- c("arrow","bar")
  names(arrows_type) <- c("simple", "T")
  
  edges <- data.frame(from = graphical_data$edge_coords$from, to = graphical_data$edge_coords$to,
                      color = graphical_data$edge_coords$col,
                      arrows.to.enabled = rep(TRUE, length(graphical_data$edge_coords$col)),
                      arrows.to.type = arrows_type[graphical_data$edge_coords$arr.type]
  )
  
  return(list(nodes = nodes, edges = edges))
  
}

#### pathway processing function ####
pathway_graph_processing <- function(pathway) {
  
  if(sum(graph::edgeData(pathway$graph, attr = "existence") == "removed") > 0) {
    removed_edges <- t(simplify2array(strsplit(names(which(graph::edgeData(pathway$graph, attr = "existence") == "removed")), split = "|", fixed = T )))
    
    pathway$graph <- graph::removeEdge(from = removed_edges[,1], to = removed_edges[,2], pathway$graph)
  }
  
  removed_nodes <- names(which(graph::nodeData(pathway$graph, attr = "existence") == "removed"))
  
  if(length(removed_nodes) > 0) {
    pathway$graph <- graph::removeNode(node = removed_nodes, pathway$graph)
  }
  
  return(pathway)
}

#### edge redirect function ####
edge_redirector <- function(from, to, g) {
  
  attrs <- AnnotationDbi::unlist2(graph::edgeData(g, from = from, to = to), recursive = F)
  attrs["change_info"] <- "reversed"
  
  g <- graph::removeEdge(from = from, to = to, graph = g)
  
  g <- graph::addEdge(from = to, to = from, graph = g)
  
  for(i in names(attrs)) {
    graph::edgeData(g, from = to, to = from, attr = i) <- attrs[[i]]
  }
  return(g)
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

#############################################
######## Gene exp mapping and PSF ###########
#############################################

#### color maker function ####
color_code <- function(values, pal1, pal2) {
  
  if(all(values > 0)) {
    calc_colors <- pal2(10)[cut(values[which(values > 0)],10)]
  } else {
    if(all(values <= 0)) {
      calc_colors <- pal1(10)[cut(values[which(values <= 0)],10)]
    } else {
      calc_colors <- c(pal1(10)[cut(values[which(values <= 0)],10)], 
                       pal2(10)[cut(values[which(values > 0)],10)])
    }
  }
  
  return(calc_colors)
}

#### psf calculation and coloring function ####
psf_signal_calculator_and_coloring <- function(entrez_fc, pathway, pathway_name, update_mod = FALSE, old_graph = NULL, no_color_mode = F) {
  
  if(update_mod) {
    psf_graph <- psf_flow_update(entrez.fc = entrez_fc,
                                 kegg.collection = old_graph[[1]][pathway_name], calculate.significance = F, sum = FALSE)
    
  } else {
    psf_graph <- psf.from.env.entrez.fc_new(entrez.fc = entrez_fc,
                                        kegg.collection = setNames(object = list(pathway), nm = pathway_name), 
                                        calculate.significance = F, sum = FALSE)
  }
  
  exp_values <- unlist(graph::nodeData(psf_graph[[1]][[pathway_name]]$graph, attr = "expression"))[which(unlist(graph::nodeData(psf_graph[[1]][[pathway_name]]$graph, attr = "type")) == "gene")]
  
  exp_values <- log(exp_values[order(exp_values)])
  
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
  
  signal_values <- log(signal_values[order(signal_values)])
  
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

chemokine_pathway_data <- kegg_data_downloader("Chemokine_signaling_pathway")

shinyServer(function(input, output, session) {
  
  #### reactive values ####
  v <- reactiveValues(image_file = kegg_node_mapper(group_graphics = group_graphics[["Chemokine_signaling_pathway"]], kegg_pathway_graphics = graphical_data_generator(chemokine_pathway_data$pathway), pathway_name = "Chemokine_signaling_pathway", pathway_image = chemokine_pathway_data$image, show_changes = F) %>% 
                        image_write(tempfile(fileext='png'), format = 'png'),
                      pathway_name = "Chemokine_signaling_pathway",
                      pathway = chemokine_pathway_data$pathway, pathway_image = chemokine_pathway_data$image, pathway_data = chemokine_pathway_data, graphical_data = graphical_data_generator(chemokine_pathway_data$pathway), 
                      selectize_state = TRUE, ordered_node_id = 1, ordered_nodes = NULL,
                      selected_node_name = setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("node_name", "gr_name", "node_id", "entrez_id")),
                      database_out = NULL, clickX = NULL, clicky = NULL, searching_state = FALSE, 
                      database_search_node = NULL, selected_interaction = NULL, edge_adding_error = "", edge_adding_error = "",
                      event_node_type = scan(file = "event_names.txt", what = "character", sep = "\t"), pathway_change_alert = 0, 
                      edge_edit_stage = FALSE, previous_row_selection = 0, allow_edge_draw = TRUE,
                      psf_and_colors = NULL, draw_color_bar = FALSE, col_legend_title = "", color_bar_lims = NULL, color_bar_psf_mode = FALSE, 
                      pathway_exp_colored = FALSE, pathway_psf_colored = FALSE, mapping_value_type = NULL, exp_uploaded = NULL, entrez_fc = NULL,
                      allow_graph_param_update = FALSE, sink_values_plot = NULL
  )
  
  #### selectize input update ####
  observe({
    updateSelectizeInput(session, "selected_pathway", 
                         choices = names(new_kegg_collection), 
                         selected = v$pathway_name, server = TRUE)
  })
  
  
  #### pathway selection by selectize input ####
  observeEvent(input$load_pathway, {
    if(input$selected_pathway != "") {
      v$pathway_name <- input$selected_pathway
      
      v$image_file <- NULL
        
      v$pathway_data <- kegg_data_downloader(input$selected_pathway)
        
      v$pathway <- v$pathway_data$pathway
        
      v$pathway_image <- v$pathway_data$image
        
      v$graphical_data <- graphical_data_generator(v$pathway)
      
      hide('vis_buttons')
      
      v$psf_and_colors <- NULL
      
      v$pathway_exp_colored <- FALSE
      v$pathway_psf_colored <- FALSE
      
      v$image_file <- kegg_node_mapper(group_graphics = group_graphics[[v$pathway_name]], kegg_pathway_graphics = v$graphical_data, pathway_name = v$pathway_name, pathway_image = v$pathway_image, show_changes = input$show_changes) %>% 
        image_write(tempfile(fileext='png'), format = 'png')
      
    }
  })
  
  #### load pathway from saved RData ####
  observeEvent(input$kegg_data, {
    
    # print(input$kegg_data$datapath)
    
    # unzip(input$kegg_data$datapath)
    
    # load(gsub(".zip", ".RData", input$kegg_data$datapath))
    # 
    # v$pathway <- saved_graph
    # 
    # v$pathway_name <- saved_pathway_name
    # 
    # v$pathway_image <- magick::image_read(gsub(".zip", ".png", input$kegg_data$datapath))
    # 
    # v$graphical_data <- graphical_data_generator(v$pathway)
    # 
    # if(tryCatch( { kegg_node_mapper(group_graphics = group_graphics[[v$pathway_name]], kegg_pathway_graphics = v$graphical_data, pathway_name = v$pathway_name, pathway_image = v$pathway_image, show_changes = input$show_changes) }
    #              , error = function(e) {TRUE})) {
    #   
    #   output$pathway_load_error <- renderText({paste("<font color=\"#ff0000\"><b>", "Something is wrong with a file", "</b></font>")})
    #   
    # } else {
    #   output$pathway_load_error <- NULL
    #   v$image_file <- kegg_node_mapper(group_graphics = group_graphics[[v$pathway_name]], kegg_pathway_graphics = v$graphical_data, pathway_name = v$pathway_name, pathway_image = v$pathway_image, show_changes = input$show_changes) %>% 
    #     image_write(tempfile(fileext='png'), format = 'png')
    # }
    
  })
  
  #### swithcher for highlighting changed items ####
  observeEvent(input$show_changes, {
    if(input$show_changes) {
      v$graphical_data <- graphical_data_generator(v$pathway, include_changes = TRUE)
    } else {
      v$graphical_data <- graphical_data_generator(v$pathway)
    }
    v$image_file <- kegg_node_mapper(group_graphics = group_graphics[[v$pathway_name]], kegg_pathway_graphics = v$graphical_data, pathway_name = v$pathway_name, pathway_image = v$pathway_image, show_changes = input$show_changes) %>% 
      image_write(tempfile(fileext='png'), format = 'png')
  })
  
  #### image click events ####
  observeEvent(input$image_click, {
    
    if(!is.null(input$multiple_choice)) {
      if(input$multiple_choice == "yes") {
        v$selected_node_name <- unique(rbind(v$selected_node_name, v$graphical_data$node_coords[which(input$image_click[[1]] <= v$graphical_data$node_coords$x_end & input$image_click[[1]] >= v$graphical_data$node_coords$x_start &
                                                                                                        input$image_click[[2]] <= v$graphical_data$node_coords$y_end & input$image_click[[2]] >= v$graphical_data$node_coords$y_start), c("node_name", "gr_name", "node_id", "entrez_id")]))
      } else {
        v$selected_node_name <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("node_name", "gr_name", "node_id", "entrez_id"))
        v$selected_node_name <- v$graphical_data$node_coords[which(input$image_click[[1]] <= v$graphical_data$node_coords$x_end & input$image_click[[1]] >= v$graphical_data$node_coords$x_start &
                                                                     input$image_click[[2]] <= v$graphical_data$node_coords$y_end & input$image_click[[2]] >= v$graphical_data$node_coords$y_start), c("node_name", "gr_name", "node_id", "entrez_id")]
      }
    } else {
      v$selected_node_name <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("node_name", "gr_name", "node_id", "entrez_id"))
      v$selected_node_name <- v$graphical_data$node_coords[which(input$image_click[[1]] <= v$graphical_data$node_coords$x_end & input$image_click[[1]] >= v$graphical_data$node_coords$x_start &
                                                                   input$image_click[[2]] <= v$graphical_data$node_coords$y_end & input$image_click[[2]] >= v$graphical_data$node_coords$y_start), c("node_name", "gr_name", "node_id", "entrez_id")]
    }
    
    # if(nrow(v$selected_node_name) < 3) {
    #   output$edge_search <- renderUI(actionButton("search_edge", label = "Search edge"))
    # } else {
    #   removeUI(selector = '#edge_search')
    # }
    
    hide('exp_change_panel')
    
    if(any(v$pathway_exp_colored, v$pathway_psf_colored)) {
      
      if(v$pathway_exp_colored) {
        colors <- v$psf_and_colors$exp_colors
      } else {
        colors <- v$psf_and_colors$psf_colors
      }
      
      v$image_file <- kegg_node_mapper(group_graphics = group_graphics[[v$pathway_name]], kegg_pathway_graphics = v$graphical_data, pathway_name = v$pathway_name, pathway_image = v$pathway_image, highlight.genes = v$selected_node_name, color.genes = colors, color_bar_psf_mode = v$color_bar_psf_mode, col_legend_title = v$col_legend_title, draw_color_bar = v$draw_color_bar, color_bar_lims = v$color_bar_lims) %>% 
        image_write(tempfile(fileext='png'), format = 'png')
    } else {
    v$image_file <- kegg_node_mapper(group_graphics = group_graphics[[v$pathway_name]], kegg_pathway_graphics = v$graphical_data, pathway_name = v$pathway_name, pathway_image = v$pathway_image, highlight.genes = v$selected_node_name, show_changes = input$show_changes) %>% 
      image_write(tempfile(fileext='png'), format = 'png')
    }
    
    
  })
  
  #### image brush events ####
  observeEvent(input$image_brush, {
    if(!input$event_node_mode) {
      if(!is.null(input$multiple_choice)) {
        if(input$multiple_choice == "yes") {
          
          v$selected_node_name <- unique(rbind(v$selected_node_name, v$graphical_data$node_coords[which(input$image_brush[1:4]$xmax >= v$graphical_data$node_coords$x_end & input$image_brush[1:4]$xmin <= v$graphical_data$node_coords$x_start &
                                                                                                          input$image_brush[1:4]$ymax >= v$graphical_data$node_coords$y_end & input$image_brush[1:4]$ymin <= v$graphical_data$node_coords$y_start), c("node_name", "gr_name", "node_id", "entrez_id")]))
          
        } else {
          v$selected_node_name <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("node_name", "gr_name", "node_id", "entrez_id"))
          v$selected_node_name <- v$graphical_data$node_coords[which(input$image_brush[1:4]$xmax >= v$graphical_data$node_coords$x_end & input$image_brush[1:4]$xmin <= v$graphical_data$node_coords$x_start &
                                                                       input$image_brush[1:4]$ymax >= v$graphical_data$node_coords$y_end & input$image_brush[1:4]$ymin <= v$graphical_data$node_coords$y_start), c("node_name", "gr_name", "node_id", "entrez_id")]
        }
      } else {
        v$selected_node_name <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("node_name", "gr_name", "node_id", "entrez_id"))
        v$selected_node_name <- v$graphical_data$node_coords[which(input$image_brush[1:4]$xmax >= v$graphical_data$node_coords$x_end & input$image_brush[1:4]$xmin <= v$graphical_data$node_coords$x_start &
                                                                     input$image_brush[1:4]$ymax >= v$graphical_data$node_coords$y_end & input$image_brush[1:4]$ymin <= v$graphical_data$node_coords$y_start), c("node_name", "gr_name", "node_id", "entrez_id")]
      }
      
      # if(nrow(v$selected_node_name) < 3) {
      #   output$edge_search <- renderUI(actionButton("search_edge", label = "Search edge"))
      # } else {
      #   removeUI(selector = '#edge_search')
      # }
      
      hide('exp_change_panel')
      
      if(any(v$pathway_exp_colored, v$pathway_psf_colored)) {
        
        if(v$pathway_exp_colored) {
          colors <- v$psf_and_colors$exp_colors
        } else {
          colors <- v$psf_and_colors$psf_colors
        }
        
        v$image_file <- kegg_node_mapper(group_graphics = group_graphics[[v$pathway_name]], kegg_pathway_graphics = v$graphical_data, pathway_name = v$pathway_name, pathway_image = v$pathway_image, highlight.genes = v$selected_node_name, color.genes = colors, color_bar_psf_mode = v$color_bar_psf_mode, col_legend_title = v$col_legend_title, draw_color_bar = v$draw_color_bar, color_bar_lims = v$color_bar_lims) %>% 
          image_write(tempfile(fileext='png'), format = 'png')
      } else {
      v$image_file <- kegg_node_mapper(group_graphics = group_graphics[[v$pathway_name]], kegg_pathway_graphics = v$graphical_data, pathway_name = v$pathway_name, pathway_image = v$pathway_image, highlight.genes = v$selected_node_name, show_changes = input$show_changes) %>% 
        image_write(tempfile(fileext='png'), format = 'png')
      }
      
    }
  })
  
  #### image hover events ####
  output$hover_text <- renderUI({
    h5(textOutput("hover_data"))
  })
  
  output$hover_data <- renderText({
    if(length(which(input$image_hover[[1]] <= v$graphical_data$node_coords$x_end & input$image_hover[[1]] >= v$graphical_data$node_coords$x_start &
                    input$image_hover[[2]] <= v$graphical_data$node_coords$y_end & input$image_hover[[2]] >= v$graphical_data$node_coords$y_start)) > 0 ) {
      
      gr_name <- v$graphical_data$node_coords[which(input$image_hover[[1]] <= v$graphical_data$node_coords$x_end & input$image_hover[[1]] >= v$graphical_data$node_coords$x_start &
                                           input$image_hover[[2]] <= v$graphical_data$node_coords$y_end & input$image_hover[[2]] >= v$graphical_data$node_coords$y_start),10]
      
      if(any(v$pathway_exp_colored, v$pathway_psf_colored)) {
        
        node_id <- v$graphical_data$node_coords[which(input$image_hover[[1]] <= v$graphical_data$node_coords$x_end & input$image_hover[[1]] >= v$graphical_data$node_coords$x_start &
                                                        input$image_hover[[2]] <= v$graphical_data$node_coords$y_end & input$image_hover[[2]] >= v$graphical_data$node_coords$y_start),8]
        
        if(v$pathway_exp_colored) {
          paste0(gr_name , " | ", "FC = ", round(v$psf_and_colors$exp_values[node_id], digits = 3))
        } else {
          if(v$pathway_psf_colored) {
            paste0(gr_name , " | ", "FC = ", round(v$psf_and_colors$exp_values[node_id], digits = 3), " Signal = ", round(v$psf_and_colors$signal_values[node_id], digits = 3))
          }
        }
      } else {
        gr_name
      }
      
    }
  })
  
  
  #### edge drawing ####
  observeEvent(input$draw_edge, {
    
    if(v$allow_edge_draw) {
      
      
      if(any(v$pathway_exp_colored, v$pathway_psf_colored)) {
        
        if(v$pathway_exp_colored) {
          colors <- v$psf_and_colors$exp_colors
        } else {
          colors <- v$psf_and_colors$psf_colors
        }
        
        v$image_file <- kegg_node_mapper(group_graphics = group_graphics[[v$pathway_name]], kegg_pathway_graphics = v$graphical_data, pathway_name = v$pathway_name, pathway_image = v$pathway_image, highlight.genes = v$selected_node_name, color.genes = colors, color_bar_psf_mode = v$color_bar_psf_mode, col_legend_title = v$col_legend_title, draw_color_bar = v$draw_color_bar, color_bar_lims = v$color_bar_lims, edge_mapping = TRUE, show_changes = input$show_changes, edge_in_mode = input$ingoing_edge) %>% 
          image_write(tempfile(fileext='png'), format = 'png')
      } else {
        v$image_file <- kegg_node_mapper(group_graphics = group_graphics[[v$pathway_name]], kegg_pathway_graphics = v$graphical_data, pathway_name = v$pathway_name, pathway_image = v$pathway_image, highlight.genes = v$selected_node_name, edge_mapping = TRUE, show_changes = input$show_changes, edge_in_mode = input$ingoing_edge) %>% 
          image_write(tempfile(fileext='png'), format = 'png')
      }
      
      # v$image_file <- kegg_node_mapper(group_graphics = group_graphics[[v$pathway_name]], kegg_pathway_graphics = v$graphical_data, pathway_name = v$pathway_name, pathway_image = v$pathway_image, highlight.genes = v$selected_node_name, edge_mapping = TRUE, show_changes = input$show_changes, edge_in_mode = input$ingoing_edge) %>% 
      #   image_write(tempfile(fileext='png'), format = 'png')
    }
    
  })
  
  #### disconnected node(s) checker ####
  observeEvent(input$check_disconnected_nodes,{
    disconnected_nodes <- disconnected_node_detector(pathway_graph_processing(v$pathway)$graph)
    
    v$selected_node_name <- v$graphical_data$node_coords[which(v$graphical_data$node_coords$node_id %in% disconnected_nodes), c("node_name", "gr_name", "node_id", "entrez_id")]
    
    v$image_file <- kegg_node_mapper(group_graphics = group_graphics[[v$pathway_name]], kegg_pathway_graphics = v$graphical_data, pathway_name = v$pathway_name, pathway_image = v$pathway_image, highlight.genes = v$selected_node_name, show_changes = input$show_changes, highlight_color = "#a3297a", opacity = 0.5) %>% 
      image_write(tempfile(fileext='png'), format = 'png')
    
  })
  
  #### clear image marks ####
  observeEvent(c(input$app_mode, input$clear), {
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
    v$selected_node_name <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("node_name", "gr_name", "node_id", "entrez_id"))
    v$image_file <- kegg_node_mapper(group_graphics = group_graphics[[v$pathway_name]], kegg_pathway_graphics = v$graphical_data, pathway_name = v$pathway_name, pathway_image = v$pathway_image, highlight.genes = v$selected_node_name, show_changes = input$show_changes) %>% 
      image_write(tempfile(fileext='png'), format = 'png')
    
    v$visnet_list <- visnet_creator(v$graphical_data)
    
    v$pathway_exp_colored = FALSE
    
    v$pathway_psf_colored = FALSE
    
  })
  
  #### psf visualization ####
  observeEvent(input$psf_higlight, {
    
    if(input$psf_mode == TRUE) {
      
      v$selected_node_name <- v$graphical_data$node_coords[which(v$graphical_data$node_coords$node_id %in% v$ordered_nodes[v$ordered_node_id]), c("node_name", "gr_name", "node_id", "entrez_id")]
      
      v$image_file <- kegg_node_mapper(group_graphics = group_graphics[[v$pathway_name]], kegg_pathway_graphics = v$graphical_data, pathway_name = v$pathway_name, pathway_image = v$pathway_image, highlight.genes = v$selected_node_name, show_changes = input$show_changes) %>% 
        image_write(tempfile(fileext='png'), format = 'png')
      
      v$ordered_node_id <- v$ordered_node_id + 1
    }
    
  })
  
  observeEvent(input$psf_mode, {
    if(input$psf_mode == TRUE) {
      v$ordered_nodes <- setdiff(names(v$pathway$order$node.order), disconnected_node_detector(v$pathway$graph)) 
    }
    v$ordered_node_id <- 1
    v$selected_node_name <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("node_name", "gr_name", "node_id", "entrez_id"))
    v$image_file <- kegg_node_mapper(group_graphics = group_graphics[[v$pathway_name]], kegg_pathway_graphics = v$graphical_data, pathway_name = v$pathway_name, pathway_image = v$pathway_image, highlight.genes = v$selected_node_name, show_changes = input$show_changes) %>% 
      image_write(tempfile(fileext='png'), format = 'png')
  })
  
  #### rigth click dialog ####
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
      } else {
        if(nrow(v$selected_node_name) == 1) {
          show('edge_search_dialog')
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
            if(nrow(as.data.table(v$graphical_data$edge_coords)[from == v$selected_node_name[1,"node_id"] & to == v$selected_node_name[2,"node_id"] | from == v$selected_node_name[2,"node_id"] & to == v$selected_node_name[1,"node_id"]]) == 0) {
              show('add_edge_button')
              hide('dir_switch')
              hide('edge_delete_button')
              hide('edit_edge')
            } else {
              if(nrow(as.data.table(v$graphical_data$edge_coords)[from == v$selected_node_name[1,"node_id"] & to == v$selected_node_name[2,"node_id"] | from == v$selected_node_name[2,"node_id"] & to == v$selected_node_name[1,"node_id"]]) == 1) {
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
            show('transition_maker_button')
            show('event_association')
            show('delete_nodes')
          }
        }
      }
      
    } else {
      
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
          
          if(paste(v$selected_node_name[c(1,2),"node_id"], collapse = ";") %in% paste(v$graphical_data$edge_coords$from, v$graphical_data$edge_coords$to, sep = ";") |
             paste(v$selected_node_name[c(2,1),"node_id"], collapse = ";") %in% paste(v$graphical_data$edge_coords$from, v$graphical_data$edge_coords$to, sep = ";")
          ) {
            show('change_edge_weight_button')
          }
        }
      }
      
    }
    
    show('dialog_content')
    
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
      
      v$database_out <- node_pairs_generator(v$database_search_node, pathway_realted = TRUE, pathway_data = v$graphical_data$node_coords, direction = input$direction_state, direction_state = input$direction_side)
      
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
  
  #### database data updater based on diretion switcher ####
  observe({
    if(v$searching_state) {
      if(nrow(v$database_search_node) == 1) {
        v$database_out <- node_pairs_generator(v$database_search_node, pathway_realted = TRUE, pathway_data = v$graphical_data$node_coords, direction = input$direction_state, direction_state = input$direction_side)
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
          highlight_genes <- rbind(v$graphical_data$node_coords[which(v$graphical_data$node_coords$node_id %in% v$database_out[input$edge_table_rows_selected, 9]), c("node_name", "gr_name", "node_id", "entrez_id")],
                                   v$graphical_data$node_coords[which(v$graphical_data$node_coords$node_id %in% v$database_out[input$edge_table_rows_selected, 10]), c("node_name", "gr_name", "node_id", "entrez_id")])
          
          v$image_file <- kegg_node_mapper(group_graphics = group_graphics[[v$pathway_name]], kegg_pathway_graphics = v$graphical_data, pathway_name = v$pathway_name, pathway_image = v$pathway_image, highlight.genes = highlight_genes, edge_mapping = TRUE, advanced_edge_ampping = TRUE, show_changes = input$show_changes) %>%
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
    updateSelectizeInput(session, "type", selected = "")
    updateSelectizeInput(session, "subtype", selected = "")
    
    v$selected_interaction <- list(source = paste(v$database_out[input$edge_table_rows_selected, "INTERACTION_DATA_SOURCE"], "database"), selected_interaction = rbind(v$graphical_data$node_coords[which(v$graphical_data$node_coords$node_id %in% v$database_out[input$edge_table_rows_selected, 9]), c("gr_name", "node_id")],
                                                                                                                                                                       v$graphical_data$node_coords[which(v$graphical_data$node_coords$node_id %in% v$database_out[input$edge_table_rows_selected, 10]), c("gr_name", "node_id")]))
  })
  
  #### adding new edge from image ####
  observeEvent(input$add_edge, {
    v$edge_edit_stage <- FALSE
    updateActionButton(session, "create_edge", label = "Create edge")
    updateSelectizeInput(session, "interaction_source", selected = "image")
    updateSelectizeInput(session, "type", selected = "")
    updateSelectizeInput(session, "subtype", selected = "")
    show('edge_creating_panel')
    show('image_edge_direction_switch')
    
    v$selected_interaction <- list(source = "database", selected_interaction = v$selected_node_name[,c("gr_name", "node_id")])
    
    v$image_file <- kegg_node_mapper(group_graphics = group_graphics[[v$pathway_name]], kegg_pathway_graphics = v$graphical_data, pathway_name = v$pathway_name, pathway_image = v$pathway_image, highlight.genes = v$selected_node_name, edge_mapping = TRUE, advanced_edge_ampping = TRUE, show_changes = input$show_changes) %>%
      image_write(tempfile(fileext='png'), format = 'png')
    
  })
  
  #### edge direction switcher for image based edge adding ####
  observeEvent(input$image_edge_direction, {
    if(input$image_edge_direction) {
      v$selected_interaction <- list(source = "database", selected_interaction = v$selected_node_name[,c("gr_name", "node_id")])
      v$image_file <- kegg_node_mapper(group_graphics = group_graphics[[v$pathway_name]], kegg_pathway_graphics = v$graphical_data, pathway_name = v$pathway_name, pathway_image = v$pathway_image, highlight.genes = v$selected_node_name, edge_mapping = TRUE, advanced_edge_ampping = TRUE, show_changes = input$show_changes) %>%
        image_write(tempfile(fileext='png'), format = 'png')
    } else {
      v$selected_interaction <- list(source = "database", selected_interaction = v$selected_node_name[2:1,c("gr_name", "node_id")])
      v$image_file <- kegg_node_mapper(group_graphics = group_graphics[[v$pathway_name]], kegg_pathway_graphics = v$graphical_data, pathway_name = v$pathway_name, pathway_image = v$pathway_image, highlight.genes = v$selected_node_name[2:1,], edge_mapping = TRUE, advanced_edge_ampping = TRUE, show_changes = input$show_changes) %>%
        image_write(tempfile(fileext='png'), format = 'png')
    }
  }, ignoreInit = T)
  
  #### editing edge attributes ####
  observeEvent(input$edit_edge, {
    v$database_edit_node <- v$selected_node_name
    v$edge_edit_stage <- TRUE
    selected_edge <- as.data.table(v$graphical_data$edge_coords)[from == v$selected_node_name[1,"node_id"] & to == v$selected_node_name[2,"node_id"] | from == v$selected_node_name[2,"node_id"] & to == v$selected_node_name[1,"node_id"]]
    
    updateSelectizeInput(session, "type", selected = unname(unlist(graph::edgeData(v$pathway$grap, selected_edge$from, selected_edge$to, attr = "subtype1"))))
    updateSelectizeInput(session, "subtype", selected = unname(unlist(graph::edgeData(v$pathway$grap, selected_edge$from, selected_edge$to, attr = "subtype2"))))
    updateSelectizeInput(session, "interaction_source", selected = "")
    updateActionButton(session, "create_edge", label = "Edit edge")
    show('edge_creating_panel')
    
    selected_interaction <- rbind(v$graphical_data$node_coords[which(v$graphical_data$node_coords$node_id == selected_edge$from), c("gr_name", "node_id")],
                                  v$graphical_data$node_coords[which(v$graphical_data$node_coords$node_id == selected_edge$to), c("gr_name", "node_id")])
    
    v$selected_interaction <- list(source = "", selected_interaction = selected_interaction)
    
  })
  
  output$edge_info <- renderText({
    paste("<b>", paste0(v$selected_interaction$selected_interaction[1,"gr_name"], " to ", v$selected_interaction$selected_interaction[2,"gr_name"]), "</b>")
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
    
    if(input$type == "" | is.null(input$interaction_source)) {
      v$edge_adding_error <- "Fill all necessary fields"
      delay(3000, v$edge_adding_error <- "")
    } else {
      if(v$edge_edit_stage) {
        graph::edgeData(v$pathway$graph, v$selected_interaction$selected_interaction[1,"node_id"], v$selected_interaction$selected_interaction[2,"node_id"], attr = "subtype1") <- input$type
        graph::edgeData(v$pathway$graph, v$selected_interaction$selected_interaction[1,"node_id"], v$selected_interaction$selected_interaction[2,"node_id"], attr = "subtype2") <- input$subtype
        graph::edgeData(v$pathway$graph, v$selected_interaction$selected_interaction[1,"node_id"], v$selected_interaction$selected_interaction[2,"node_id"], attr = "change_info") <- "attr_change"
        data_source <- paste(input$interaction_source, collapse = ", ")
        graph::edgeData(v$pathway$graph, v$selected_interaction$selected_interaction[1,"node_id"], v$selected_interaction$selected_interaction[2,"node_id"], attr = "data_source") <- data_source
        v$pathway_change_alert <- v$pathway_change_alert + 1
        hide('edge_creating_panel')
      } else {
        data_source <- input$interaction_source
        data_source[which(data_source == "database")] <- v$selected_interaction$source
        data_source <- paste(data_source, collapse = ", ")
        if(!graph::isAdjacent(v$pathway$graph, from = v$selected_interaction$selected_interaction[1,"node_id"], to = v$selected_interaction$selected_interaction[2,"node_id"])) {
          v$pathway$graph = graph::addEdge(v$selected_interaction$selected_interaction[1,"node_id"], v$selected_interaction$selected_interaction[2,"node_id"], v$pathway$graph)
        }
        graph::edgeData(v$pathway$graph, v$selected_interaction$selected_interaction[1,"node_id"], v$selected_interaction$selected_interaction[2,"node_id"], attr = "subtype1") <- input$type
        graph::edgeData(v$pathway$graph, v$selected_interaction$selected_interaction[1,"node_id"], v$selected_interaction$selected_interaction[2,"node_id"], attr = "subtype2") <- input$subtype
        graph::edgeData(v$pathway$graph, v$selected_interaction$selected_interaction[1,"node_id"], v$selected_interaction$selected_interaction[2,"node_id"], attr = "existence") <- "added"
        graph::edgeData(v$pathway$graph, v$selected_interaction$selected_interaction[1,"node_id"], v$selected_interaction$selected_interaction[2,"node_id"], attr = "data_source") <- data_source
        v$pathway_change_alert <- v$pathway_change_alert + 1
        hide('edge_creating_panel')
        hide('image_edge_direction_switch')
      }
    }
    
  })
  
  #### edge direction switch ####
  observeEvent(input$edge_direction_switch, {
    index <- unique(which(v$graphical_data$edge_coords$from %in% v$selected_node_name[,"node_id"] & v$graphical_data$edge_coords$to %in% v$selected_node_name[,"node_id"]))
    
    from <- v$graphical_data$edge_coords$from[index]
    to <- v$graphical_data$edge_coords$to[index]
    
    v$pathway$graph <- edge_redirector(from = from, to = to, g = v$pathway$graph)
    
    v$pathway_change_alert <- v$pathway_change_alert + 1
  })
  
  #### delete edge ####
  observeEvent(input$edge_delete, {
    index <- unique(which(v$graphical_data$edge_coords$from %in% v$selected_node_name[,"node_id"] & v$graphical_data$edge_coords$to %in% v$selected_node_name[,"node_id"]))
    
    from <- v$graphical_data$edge_coords$from[index]
    to <- v$graphical_data$edge_coords$to[index]
    
    # v$pathway$graph <- graph::removeEdge(from = from, to = to, v$pathway$graph)
    
    graph::edgeData(v$pathway$graph, from = from, to = to, attr = "existence") <- "removed"
    
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
  observeEvent(input$associate_with_event_node, {
    if("bio_event" %in% unlist(graph::nodeData(v$pathway$graph, v$selected_node_name[,"node_id"], attr = "type"))) {
      event_ids <- names(which(unlist(graph::nodeData(v$pathway$graph, v$selected_node_name[,"node_id"], attr = "type")) == "bio_event"))
      gene_ids <- names(which(unlist(graph::nodeData(v$pathway$graph, v$selected_node_name[,"node_id"], attr = "type")) != "bio_event"))
      
      if(length(event_ids) > 0 & length(event_ids) < 2 & length(gene_ids) > 0) {
        v$pathway$graph <- graph::addEdge(gene_ids, event_ids, v$pathway$graph)
        graph::edgeData(v$pathway$graph, gene_ids, event_ids, attr = "data_source") <- "image"
        graph::edgeData(v$pathway$graph, gene_ids, event_ids, attr = "subtype1") <- "activation"
        graph::edgeData(v$pathway$graph, gene_ids, event_ids, attr = "subtype2") <- "indirect effect"
        
        v$pathway_change_alert <- v$pathway_change_alert + 1
      }
    } else {
      v$allow_edge_draw <- FALSE
      show("event_node_adding_panel")
    }
    
  })
  
  #### adding new node ####
  observeEvent(input$add_node, {
    
    if(is.null(input$image_brush[1:4])) {
      v$node_adding_error <- "Please specify node location"
      delay(3000, v$node_adding_error <- NULL)
    } else {
      new_node_num <- as.character(max(as.integer(v$pathway$graph@nodes)) + 1)
      v$pathway$graph <- graph::addNode(new_node_num, v$pathway$graph)
      graph::nodeData(v$pathway$graph, new_node_num, attr = "label") <- input$event_node_name
      graph::nodeData(v$pathway$graph, new_node_num, attr = "kegg.name") <- input$event_node_name
      graph::nodeData(v$pathway$graph, new_node_num, attr = "data_source") <- input$node_source
      graph::nodeData(v$pathway$graph, new_node_num, attr = "type") <- "bio_event"
      graph::nodeData(v$pathway$graph, new_node_num, attr = "kegg.id") <- new_node_num
      graph::nodeData(v$pathway$graph, new_node_num, attr = "kegg.gr.x") <- input$image_brush[[1]] + (input$image_brush[[2]] - input$image_brush[[1]])*0.5
      graph::nodeData(v$pathway$graph, new_node_num, attr = "kegg.gr.y") <- input$image_brush[[3]] + (input$image_brush[[4]] - input$image_brush[[3]])*0.5
      graph::nodeData(v$pathway$graph, new_node_num, attr = "kegg.gr.width") <- input$image_brush[[2]] - input$image_brush[[1]]
      graph::nodeData(v$pathway$graph, new_node_num, attr = "kegg.gr.height") <- input$image_brush[[4]] - input$image_brush[[3]]
      
      print(v$selected_node_name[,"node_id"])
      v$pathway$graph <- graph::addEdge(v$selected_node_name[,"node_id"], new_node_num, v$pathway$graph)
      graph::edgeData(v$pathway$graph, v$selected_node_name[,"node_id"], new_node_num, attr = "data_source") <- input$node_source
      graph::edgeData(v$pathway$graph, v$selected_node_name[,"node_id"], new_node_num, attr = "subtype1") <- "activation"
      graph::edgeData(v$pathway$graph, v$selected_node_name[,"node_id"], new_node_num, attr = "subtype2") <- "indirect effect"
      
      v$pathway_change_alert <- v$pathway_change_alert + 1
      hide('event_node_adding_panel')
      v$allow_edge_draw <- TRUE
      updateCheckboxInput(session, "event_node_mode", value = FALSE)
    }
    
    if(!(input$event_node_name %in% v$event_node_type)) {
      v$event_node_type <- unique(c(scan(file = "event_names.txt", what = "character", sep = "\t"), input$event_node_name))
      write(v$event_node_type, file = "event_names.txt", sep = "\t")
    }
    
  })
  
  observeEvent(input$close_node_panel, {
    hide('event_node_adding_panel')
    v$allow_edge_draw <- TRUE
    updateCheckboxInput(session, "event_node_mode", value = FALSE)
  })
  
  #### make node as a transition node ####
  observeEvent(input$make_transition,{
    
    graph::nodeData(v$pathway$graph, v$selected_node_name[,"node_id"], attr = "type") <- "transition_gene"
    
    graph::nodeData(v$pathway$graph, v$selected_node_name[,"node_id"], attr = "change_info") <- "class_change"
    
    graph::nodeData(v$pathway$graph, v$selected_node_name[,"node_id"], attr = "data_source") <- "image"
    
    v$pathway_change_alert <- v$pathway_change_alert + 1
  })
  
  #### delete node ####
  observeEvent(input$node_delete, {
    
    index <- unique(which(v$graphical_data$edge_coords$from %in% v$selected_node_name[,"node_id"] | v$graphical_data$edge_coords$to %in% v$selected_node_name[,"node_id"]))
    
    from <- v$graphical_data$edge_coords$from[index]
    to <- v$graphical_data$edge_coords$to[index]
    
    graph::edgeData(v$pathway$graph, from = from, to = to, attr = "existence") <- "removed"
    
    graph::nodeData(v$pathway$graph, n = v$selected_node_name[,"node_id"], attr = "existence") <- "removed"
    
    v$pathway_change_alert <- v$pathway_change_alert + 1
  })
  
  
  #### observer for saving pathway graph after any changies ####
  observeEvent(v$pathway_change_alert, {
    if(v$pathway_change_alert != 0) {
      v$pathway$sink.nodes <- psf::determine.sink.nodes(pathway_graph_processing(v$pathway))
      
      v$pathway$order <-  psf::order.nodes(pathway_graph_processing(v$pathway)$graph)
      
      v$pathway$graph <- psf::set.edge.impacts(v$pathway$graph)
      
      v$graphical_data <- graphical_data_generator(v$pathway)
      
      v$selected_node_name <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("node_name", "gr_name", "node_id", "entrez_id"))
      
      v$image_file <- kegg_node_mapper(group_graphics = group_graphics[[v$pathway_name]], kegg_pathway_graphics = v$graphical_data, pathway_name = v$pathway_name, pathway_image = v$pathway_image, show_changes = input$show_changes) %>%
        image_write(tempfile(fileext='png'), format = 'png')
      
      # saveRDS(v$pathway, file = paste0("edited_pathway_graphs/", v$pathway_name, ".RDS"))
    }
  }, suspended = F)
  
  #### downlod handler for edited all pathways ####
  output$download_pathway <- downloadHandler(
    filename = function(){
      paste0(v$pathway_name, ".RData")
    },
    content = function(file) {
      
      saved_graph <- v$pathway
      saved_pathway_name <- v$pathway_name
      saved_image <- v$pathway_image
      
      save(saved_graph, saved_pathway_name, saved_image, file = file)
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
      v$allow_graph_param_update <- FALSE
      
    } else {
      output$fc_table_load_error <- NULL
      v$exp_uploaded <- as.matrix(suppressWarnings(fread(input$file$datapath, stringsAsFactors = F)))
      
      rownames(v$exp_uploaded) <- v$exp_uploaded[,1]
      
      v$entrez_fc <- v$exp_uploaded[,-1, drop = F]
      
      print(v$entrez_fc)
      
      # psf_signal_calculator_and_coloring(entrez_fc = v$entrez_fc, pathway = v$pathway, pathway_name = v$pathway_name, update_mod = FALSE)
      
      if(isTRUE(tryCatch( { psf_signal_calculator_and_coloring(entrez_fc = v$entrez_fc, pathway = v$pathway, pathway_name = v$pathway_name, update_mod = FALSE) }
                          , error = function(e1) {TRUE}))) {
        
        output$fc_table_load_error <- renderText({paste("<font color=\"#ff0000\"><b>", "Something is wrong with a file", "</b></font>")})
        hide('vis_buttons')
        v$allow_graph_param_update <- FALSE
        
      } else {
        output$fc_table_load_error <- NULL
        v$psf_and_colors <- psf_signal_calculator_and_coloring(entrez_fc = v$entrez_fc, pathway = v$pathway,
                                                               pathway_name = v$pathway_name, update_mod = FALSE)
        v$allow_graph_param_update <- TRUE
        
        updateSelectizeInput(session, "node_influence", 
                             choices = v$graphical_data$node_coords$gr_name[which(v$graphical_data$node_coords$node_class == "gene")], 
                             selected = v$graphical_data$node_coords$gr_name[which(v$graphical_data$node_coords$node_class == "gene")][1], server = TRUE)
        
        show('vis_buttons')
        
        
      }
      
    }
    
  })
  
  #### FC mapping ####
  observeEvent(input$map_fc_values, {
    
    v$psf_and_colors <- psf_signal_calculator_and_coloring(entrez_fc = v$entrez_fc, pathway = v$pathway,
                                                           pathway_name = v$pathway_name, update_mod = FALSE)
    
    v$draw_color_bar <- TRUE
    
    v$color_bar_lims <- range(v$psf_and_colors$exp_values)
    
    v$col_legend_title <- "log FC"
    
    v$color_bar_psf_mode <- TRUE
    
    v$image_file <- kegg_node_mapper(group_graphics = group_graphics[[v$pathway_name]], kegg_pathway_graphics = v$graphical_data, pathway_name = v$pathway_name, pathway_image = v$pathway_image, color.genes = v$psf_and_colors$exp_colors, color_bar_psf_mode = v$color_bar_psf_mode, col_legend_title = v$col_legend_title, color_bar_lims = v$color_bar_lims, draw_color_bar = v$draw_color_bar) %>% 
      image_write(tempfile(fileext='png'), format = 'png')
    
    v$visnet_list <- visnet_creator(v$graphical_data, node_colors = v$psf_and_colors$exp_colors)
    
    v$pathway_exp_colored = TRUE
    
    v$pathway_psf_colored = FALSE
    
    v$mapping_value_type <- "FC"
  })
  
  #### node expression change observer ####
  observeEvent(input$change_exp, {
    
    updateNumericInput(session, "exp_slider", label = paste0(v$selected_node_name[,"gr_name"], " FC value"),
                       value = round(unlist(unname(graph::nodeData(v$psf_and_colors$psf_graph[[1]][[v$pathway_name]]$graph, v$selected_node_name[,"node_id"],attr = "expression"))), digits = 3), 
                       min = 0, max = 10, step = 0.001)
    
    show('exp_change_panel')
    
  })
  
  observeEvent(input$close_slider, {
    v$searching_state <- FALSE
    hide('exp_change_panel')
  })
  
  #### PSF calculation and node coloring ####
  observeEvent(input$psf_run, {
    
    v$psf_and_colors <- psf_signal_calculator_and_coloring(entrez_fc = v$entrez_fc, pathway = v$pathway,
                                                           pathway_name = v$pathway_name, update_mod = FALSE)
    
    v$draw_color_bar <- TRUE
    
    v$color_bar_lims <- range(v$psf_and_colors$signal_values)
    
    v$col_legend_title <- "Signal log value"
    
    v$color_bar_psf_mode <- TRUE
    
    
    v$image_file <- kegg_node_mapper(group_graphics = group_graphics[[v$pathway_name]], kegg_pathway_graphics = v$graphical_data, pathway_name = v$pathway_name, pathway_image = v$pathway_image, color.genes = v$psf_and_colors$psf_colors, color_bar_psf_mode = v$color_bar_psf_mode, col_legend_title = v$col_legend_title, color_bar_lims = v$color_bar_lims, draw_color_bar = v$draw_color_bar)
    
    
    sink_plot <- ggplot(v$psf_and_colors$sink_signals, aes(x=sink_name, y=signal, fill=signal)) +
      geom_boxplot(color="black", lwd=0.2, outlier.shape=NA) +
      geom_point(color = v$psf_and_colors$sink_signals$dot_color) +
      geom_jitter(color = v$psf_and_colors$sink_signals$dot_color, width = 0.2) +
      guides(fill=FALSE) +
      coord_flip() +
      theme_bw()
    
    ggsave("sink_boxplot.png", plot = sink_plot, device = "png", path = NULL,
           scale = 1, width = 400, height = magick::image_info(v$image_file)$height, units = "px",
           dpi = 96, limitsize = TRUE)
    
    gg_plot_img <- image_read('sink_boxplot.png')
    
    file.remove('sink_boxplot.png')
      
    image_append(c(v$image_file, gg_plot_img)) %>%
      image_write(tempfile(fileext='png'), format = 'png')  
    
    v$visnet_list <- visnet_creator(v$graphical_data, node_colors = v$psf_and_colors$psf_colors)
    
    v$pathway_psf_colored = TRUE
    
    v$pathway_exp_colored = FALSE
    
    v$mapping_value_type <- "PSF"
    
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
    
    v$image_file <- kegg_node_mapper(group_graphics = group_graphics[[v$pathway_name]], kegg_pathway_graphics = v$graphical_data, pathway_name = v$pathway_name, pathway_image = v$pathway_image, highlight.genes = v$selected_node_name, color.genes = v$psf_and_colors$psf_colors, color_bar_psf_mode = v$color_bar_psf_mode, col_legend_title = v$col_legend_title, color_bar_lims = v$color_bar_lims, draw_color_bar = v$draw_color_bar) %>% 
      image_write(tempfile(fileext='png'), format = 'png')
    
    v$visnet_list <- visnet_creator(v$graphical_data, node_colors = v$psf_and_colors$psf_colors)
    
    sink_gg_plot <- ggplot(v$psf_and_colors$sink_signals, aes(x=sink_name, y=signal, fill=signal)) +
      geom_boxplot(color="black", lwd=0.2, outlier.shape=NA) +
      geom_point(color = v$psf_and_colors$sink_signals$dot_color) +
      geom_jitter(color = v$psf_and_colors$sink_signals$dot_color, width = 0.2) +
      guides(fill=FALSE) +
      coord_flip() +
      theme_bw()
    
    v$sink_values_plot <- ggplotly(sink_gg_plot, width = 400, height = magick::image_info(v$pathway_image)$height)
    
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
    
    v$image_file <- kegg_node_mapper(group_graphics = group_graphics[[v$pathway_name]], kegg_pathway_graphics = v$graphical_data, pathway_name = v$pathway_name, pathway_image = v$pathway_image, highlight.genes = v$selected_node_name, color.genes = v$psf_and_colors$psf_colors, color_bar_psf_mode = v$color_bar_psf_mode, col_legend_title = v$col_legend_title, color_bar_lims = v$color_bar_lims, draw_color_bar = v$draw_color_bar) %>% 
      image_write(tempfile(fileext='png'), format = 'png')
    
    v$visnet_list <- visnet_creator(v$graphical_data, node_colors = v$psf_and_colors$psf_colors)
    
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
    
    return(list(src = v$image_file, contentType = "image/png", alt = "pathway_image"))
    
  }, deleteFile = TRUE)
  
  observe({
    v$visnet_list <- visnet_creator(v$graphical_data)
  })
  
  #### sink data plot rendering ####
  # output$sink_plot <- renderPlotly({
  #   v$sink_values_plot
  # })
  
  
  #### visnet rendering ####
  output$visnet <- renderVisNetwork({
    visNetwork(nodes = v$visnet_list$nodes, edges = v$visnet_list$edges, width = "100%", height = "800px") %>% 
      visIgraphLayout(layout = "layout_nicely") %>% 
      visInteraction(navigationButtons = TRUE, multiselect = T) %>%
      visExport()
  })
  
  #### visnet search ####
  observeEvent(input$search_go, {
    id <- v$graphical_data$node_coords$node_id[grep(input$network_search, v$graphical_data$node_coords$gr_name, ignore.case = T)]
    if(length(id) > 0) {
      visNetworkProxy("visnet") %>%
        visFocus(id = id[1], scale = 2)
    }
  })
  
})
