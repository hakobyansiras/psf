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


kegg_designer <- function(group_graphics, node_graphics, pathway_image,
                          psf_output, color_bar_psf_mode = F, col_legend_title, plot_type = "boxplot", adj = c(0.48, 1)) {
  
  color.genes <- psf_output$psf_colors
  color_bar_lims <- range(psf_output$mean_signal_values)
  sink_signals <- psf_output$sink_signals
  
  img <- magick::image_draw(pathway_image)
  
  
  ### node exp coloring
  node_graphics$x_center <- node_graphics$x_start + (node_graphics$x_end - node_graphics$x_start)/2
  
  node_graphics$y_center <- node_graphics$y_start + (node_graphics$y_end - node_graphics$y_start)/2
  
  rownames(color.genes) <- color.genes$node_id
  
  rownames(node_graphics) <- node_graphics$node_id
  
  sink_node_graphics <- node_graphics[which(node_graphics$sink),]
  
  
  if(any(node_graphics$node_id %in% color.genes$node_id)) {
    coloring_set <- node_graphics[color.genes$node_id,]
    graphics::rect(coloring_set$x_start,
         coloring_set$y_start - 1,
         coloring_set$x_end,
         coloring_set$y_end,
         # border = coloring_set$border_color,
         border = NA,
         # lty = coloring_set$lty_type,
         lwd=2,
         col = grDevices::adjustcolor( color.genes$col, alpha.f = 1)
    )
    
    graphics::text(x = coloring_set$x_center,
         y = coloring_set$y_start,
         labels = coloring_set$gr_name,
         col = color.genes$text_col, adj = c(0,0.2) + adj)
    
  }
  
  graphics::text(x = sink_node_graphics$x_end + 10,
       y = sink_node_graphics$y_center - 30, cex = 3,
       labels = rep("*", nrow(sink_node_graphics)),
       col = rep("#9ACD32", nrow(sink_node_graphics)), adj = c(0,0.2) + adj)
  
  
  ### scale color bar
  if(color_bar_psf_mode) {
    color_legend_maker(x = magick::image_info(img)$width - 230, y = 50, leg = 200, cols = c(pal1(10), pal2(10)), title = col_legend_title, lims = color_bar_lims, digits=3, prompt=FALSE,
                       lwd=4, outline=TRUE, subtitle = "", fsize = 1.3)
  } else {
    # color_legend_maker(x = magick::image_info(img)$width - 230, y = 50, leg = 200, cols = exp_pal(20), title = col_legend_title, lims = exp_color_all$exp_lims, digits=3, prompt=FALSE,
    #                    lwd=4, outline=TRUE, subtitle = "", fsize = 1.3)
    # temporarely turned off exp coloring legend
    color_legend_maker(x = magick::image_info(img)$width - 230, y = 50, leg = 200, cols = c(pal1(10), pal2(10)), title = col_legend_title, lims = color_bar_lims, digits=3, prompt=FALSE,
                       lwd=4, outline=TRUE, subtitle = "", fsize = 1.3)
  }
  
  text(x = c(magick::image_info(img)$width - 88, magick::image_info(img)$width - 30),
       y = c(70, 65), cex = c(1.5, 3),
       labels = c("Sink node", "*"),
       col = c("#000000", "#9ACD32"), adj = c(0,0.2) + adj)
  
  
  ## color grop nodes
  if(length(group_graphics) > 0 ) {
    lapply(group_graphics, function(z) {
      graphics::rect( z$kegg.gr.x-z$kegg.gr.width*0.5, 
            z$kegg.gr.y+z$kegg.gr.height*0.5, 
            z$kegg.gr.x+z$kegg.gr.width*0.5, 
            z$kegg.gr.y-z$kegg.gr.height*0.5, 
            border = "yellow", lty = "dashed", lwd=2)
    })
  }
  
  
  dev.off()
  
  if(plot_type == "boxplot") {
    
    if(ncol(psf_output$exp_values_all) == 1) {
      
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
    
  } else {
    
    png(filename = "sink_plot.png", width = 400, height = magick::image_info(img)$height, res = 96)
    
    gplots::heatmap.2(psf_output$sink_values_all, col = c(pal1(10), pal2(10)), trace = "none",
                      Rowv = FALSE, Colv = FALSE, 
                      lmat = rbind(c(5,3,4), c(2,1,1)),
                      lhei = c(0.1, 4),
                      lwid = c(0.1, 4, 0.75),
                      dendrogram = 'none', key = FALSE, 
                      trace = "none", ylab = "Sink values",
                      margin = c(10, 8)
    )
    
    dev.off()
  }
  
  sink_plot_img <- image_read('sink_plot.png')
  
  file.remove('sink_plot.png')
  
  combined_img <- image_append(c(img, sink_plot_img))
  
  return(combined_img)
  
}


color_legend_maker <- function (leg, cols, title = NULL, lims = c(0, 1), digits = 1, 
                                prompt = TRUE, lwd = 4, outline = TRUE, ...) 
{
  if (prompt) {
    cat("Click where you want to draw the bar\n")
    flush.console()
    x <- unlist(locator(1))
    y <- x[2]
    x <- x[1]
  }
  else {
    if (hasArg(x)) 
      x <- list(...)$x
    else x <- 0
    if (hasArg(y)) 
      y <- list(...)$y
    else y <- 0
  }
  if (hasArg(fsize)) 
    fsize <- list(...)$fsize
  else fsize <- 1
  if (hasArg(subtitle)) 
    subtitle <- list(...)$subtitle
  else subtitle <- NULL
  if (hasArg(direction)) 
    direction <- list(...)$direction
  else direction <- "rightwards"
  if (direction %in% c("rightwards", "leftwards")) {
    X <- x + cbind(0:(length(cols) - 1)/length(cols), 1:length(cols)/length(cols)) * 
      (leg)
    if (direction == "leftwards") {
      X <- X[nrow(X):1, ]
      if (!is.null(lims)) 
        lims <- lims[2:1]
    }
    Y <- cbind(rep(y, length(cols)), rep(y, length(cols)))
  }
  else if (direction %in% c("upwards", "downwards")) {
    Y <- y + cbind(0:(length(cols) - 1)/length(cols), 1:length(cols)/length(cols)) * 
      (leg)
    if (direction == "downwards") {
      X <- X[nrow(X):1, ]
      if (!is.null(lims)) 
        lims <- lims[2:1]
    }
    X <- cbind(rep(x, length(cols)), rep(x, length(cols)))
  }
  if (outline) 
    lines(c(X[1, 1], X[nrow(X), 2]), c(Y[1, 1], Y[nrow(Y), 
                                                  2]), lwd = lwd + 2, lend = 2)
  for (i in 1:length(cols)) lines(X[i, ], Y[i, ], col = cols[i], 
                                  lwd = lwd, lend = 2)
  if (direction %in% c("rightwards", "leftwards")) {
    if (!is.null(lims)) 
      text(x = x, y = y, round(lims[1], digits), pos = 3, 
           cex = fsize)
    if (!is.null(lims)) 
      text(x = x + leg, y = y, round(lims[2], digits), 
           pos = 3, cex = fsize)
    if (is.null(title)) 
      title <- "P(state=1)"
    text(x = (2 * x + leg)/2, y = y, title, pos = 3, cex = fsize)
    if (is.null(subtitle)) 
      text(x = (2 * x + leg)/2, y = y, paste("length=", 
                                             round(leg, 3), sep = ""), pos = 1, cex = fsize)
    else text(x = (2 * x + leg)/2, y = y, subtitle, pos = 1, 
              cex = fsize)
  }
  else if (direction %in% c("upwards", "downwards")) {
    if (!is.null(lims)) 
      text(x = x, y = y - 0.02 * diff(par()$usr[3:4]), 
           round(lims[1], digits), pos = 1, cex = fsize)
    if (!is.null(lims)) 
      text(x = x, y = y + leg + 0.02 * diff(par()$usr[3:4]), 
           round(lims[2], digits), pos = 3, cex = fsize)
    if (is.null(title)) 
      title <- "P(state=1)"
    text(x = x - 0.04 * diff(par()$usr[1:2]), y = (2 * y + 
                                                     leg)/2, title, pos = 3, cex = fsize, srt = 90)
    if (is.null(subtitle)) 
      text(x = x + 0.04 * diff(par()$usr[1:2]), y = (2 * 
                                                       y + leg)/2, paste("length=", round(leg, 3), 
                                                                         sep = ""), pos = 1, srt = 90, cex = fsize)
    else text(x = x + 0.04 * diff(par()$usr[1:2]), y = (2 * 
                                                          y + leg)/2, subtitle, pos = 1, cex = fsize, srt = 90)
  }
}
