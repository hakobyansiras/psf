color_code <- function(values, pal1, pal2, log_scale = TRUE) {
  
  if(log_scale) {
    center_val <- 0
  } else {
    center_val <- 1
  }
  
  if(all(values > center_val)) {
    calc_colors <- pal2(10)[cut(values[which(values > center_val)],10)]
  } else {
    if(all(values <= center_val)) {
      calc_colors <- pal1(10)[cut(values[which(values <= center_val)],10)]
    } else {
      calc_colors <- c(pal1(10)[cut(values[which(values <= center_val)],10)], 
                       pal2(10)[cut(values[which(values > center_val)],10)])
    }
  }
  
  return(calc_colors)
}


psf_signal_calculator_and_coloring <- function(entrez_fc, pathway, pathway_name, calculate_significance = FALSE, no_color_mode = F) {
  
  psf_graph <- psf.from.env.entrez.fc(entrez.fc = entrez_fc,
                                      kegg.collection = setNames(object = list(pathway), nm = pathway_name), 
                                      calculate.significance = calculate_significance, sum = FALSE)
  
  exp_values_all <- sapply(psf_graph, function(x) {
    
    unlist(graph::nodeData(x[[pathway_name]]$graph, attr = "expression"))[which(unlist(graph::nodeData(x[[pathway_name]]$graph, attr = "type")) == "gene")]
    
  })
  
  colnames(exp_values_all) <- colnames(entrez_fc)
  
  if(ncol(entrez_fc) == 1) {
    
    mean_exp_values <- round(log(drop(exp_values_all)[order(drop(exp_values_all))] + 0.00001), digits = 5)
    
  } else {
    
    ### averaging exp FC values across sample for pathway coloring 
    mean_exp_values <- rowMeans(sapply(psf_graph, function(x) {
      
      unlist(graph::nodeData(x[[pathway_name]]$graph, attr = "expression"))[which(unlist(graph::nodeData(x[[pathway_name]]$graph, attr = "type")) == "gene")]
      
    }))
    
    mean_exp_values <- round(log(mean_exp_values[order(mean_exp_values)] + 0.00001), digits = 5)
    
  }
  
  
  
  if(no_color_mode) {
    exp_colors <- NULL
  } else {
    
    exp_colors <- color_code(values = mean_exp_values, pal1 = pal1, pal2 = pal2)
    
    # exp_colors <- c(pal1(length(unname(mean_exp_values[which(mean_exp_values <= 1)])))[rank(mean_exp_values[which(mean_exp_values <= 1)])],
    #                 pal2(length(unname(mean_exp_values[which(mean_exp_values > 1)])))[rank(mean_exp_values[which(mean_exp_values > 1)])]
    # )
    
    # exp_colors <- pal(length(unname(mean_exp_values)))[rank(mean_exp_values)]
    
    exp_colors = data.frame(node_id = names(mean_exp_values)[c(which(mean_exp_values <= 0), which(mean_exp_values > 0))], 
                            col = exp_colors,
                            text_col = unname(sapply(exp_colors, function(x) {c( "black", "white")[  1+(sum( col2rgb(x) *c(299, 587,114))/1000 < 123) ]})),
                            stringsAsFactors = F
    )
  }
  
  
  signal_values_all <- sapply(psf_graph, function(x) {
    
    # unlist(graph::nodeData(x[[pathway_name]]$graph, attr = "signal"))[which(unlist(graph::nodeData(x[[pathway_name]]$graph, attr = "type")) == "gene")]
    
    unlist(graph::nodeData(x[[pathway_name]]$graph, attr = "signal"))[which(unlist(graph::nodeData(x[[pathway_name]]$graph, attr = "type")) != "map")]
    
  })
  
  signal_values_all <- round(log(signal_values_all) + 0.00001, digits = 5)
  
  colnames(signal_values_all) <- colnames(entrez_fc)
  
  ### averaging PSF signal values across sample for pathway coloring
  mean_signal_values <- rowMeans(sapply(psf_graph, function(x) {
    
    unlist(graph::nodeData(x[[pathway_name]]$graph, attr = "signal"))[which(unlist(graph::nodeData(x[[pathway_name]]$graph, attr = "type")) == "gene")]
    
  }))
  
  
  mean_signal_values <- round(log(mean_signal_values[order(mean_signal_values)] + 0.00001), digits = 5)
  
  if(no_color_mode) {
    psf_colors <- NULL
  } else {
    
    psf_colors <- color_code(values = mean_signal_values, pal1 = pal1, pal2 = pal2)
    
  }
  
  
  # psf_colors <- c(pal1(length(unname(mean_signal_values[which(mean_signal_values <= 1)])))[rank(mean_signal_values[which(mean_signal_values <= 1)])],
  #                 pal2(length(unname(mean_signal_values[which(mean_signal_values > 1)])))[rank(mean_signal_values[which(mean_signal_values > 1)])]
  # )
  
  # psf_colors <- pal(length(unname(mean_signal_values)))[rank(mean_signal_values)]
  
  psf_colors = data.frame(node_id = names(mean_signal_values)[c(which(mean_signal_values <= 0), which(mean_signal_values > 0))], 
                          col = psf_colors,
                          text_col = unname(sapply(psf_colors, function(x) {c( "black", "white")[  1+(sum( col2rgb(x) *c(299, 587,114))/1000 < 123) ]})),
                          stringsAsFactors = F
  )
  
  
  sink_names <- unlist(graph::nodeData(pathway$graph, pathway$sink.nodes, attr = "label"))
  
  sink_order_by_coord <- order(as.integer(unname(unlist(graph::nodeData(pathway$graph, pathway$sink.nodes, attr = "kegg.gr.y")))), decreasing = T)
  
  
  repeating_sinks <- table(sink_names)[which(table(sink_names) > 1)]
  
  if(length(repeating_sinks) > 0) {
    for(i in 1:length(repeating_sinks)) {
      
      sink_names[which(sink_names == names(repeating_sinks[i]))] <- sapply(1:repeating_sinks[i], function(x) {paste0(paste(rep(" ", x - 1), collapse = ""), names(repeating_sinks[i]))})
      
    }
  }
  
  # sink_names <- setNames(nm = names(sink_names), object = paste0(sink_names, "_", 1:length(sink_names)))
  
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
  
  
  sink_values_all <- signal_values_all[pathway$sink.nodes,]
  
  if(ncol(entrez_fc) == 1) {
    
    sink_values_all <- as.matrix(sink_values_all)
    
    rownames(sink_values_all) <- sink_names
    
    sink_values_all <- sink_values_all[rev(sink_names[sink_order_by_coord]),, drop = F]
    
    rownames(sink_values_all) <- gsub('_([0-9])', "", rownames(sink_values_all))
    
  } else {
    
    rownames(sink_values_all) <- sink_names
    
    sink_values_all <- sink_values_all[rev(sink_names[sink_order_by_coord]),]
    
    rownames(sink_values_all) <- gsub('_([0-9])', "", rownames(sink_values_all))
    
  }
  
  
  
  return(list(psf_colors = psf_colors, exp_colors = exp_colors, 
              mean_exp_values = mean_exp_values, mean_signal_values = mean_signal_values, 
              exp_values_all = exp_values_all, signal_values_all = signal_values_all, 
              sink_signals = signal_data, sink_values_all = sink_values_all,
              psf_graph = psf_graph))
  
}