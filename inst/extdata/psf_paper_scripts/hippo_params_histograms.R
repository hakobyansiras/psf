library(ggplot2)
library(hrbrthemes)
library(gridExtra)
library(igraph)
library(psf)

### loading unprocessed kgml files
kegg_collection_unporecessed <- generate.kegg.collection.from.kgml(list.files(system.file("extdata", "kgmls", package="psf"), full.names = T), sink.nodes = T)

### loading curated networks
load(system.file("extdata", "edited_pathways_new.RData", package="psf"))

#### network parameters calculation function ####
network_param_calculator <- function(old_pathway, new_pathway) {
  
  g_graphnel_new <- graph::removeNode(names(unlist(graph::nodeData(new_pathway$graph, attr = "kegg.type")))[which(unlist(graph::nodeData(new_pathway$graph, attr = "kegg.type")) == "map")], new_pathway$graph) 
  
  g_new <- igraph.from.graphNEL(g_graphnel_new)
  
  params_new <- list(Shortest_paths = shortest_paths_from_input_to_sink_nodes(new_pathway, sink_nodes = new_pathway$sink.nodes)[which(!is.infinite(shortest_paths_from_input_to_sink_nodes(new_pathway, sink_nodes = new_pathway$sink.nodes)), arr.ind = T)], 
                     Degree_in = igraph::degree(g_new, mode = "in"),
                     Degree_out = igraph::degree(g_new, mode = "out"),
                     Closeness_all = closeness(g_new, mode = "all"),
                     Betweenness_all = igraph::betweenness(g_new, directed = T),
                     Hub = hub_score(g_new, weights=NA, scale = T)$vector,
                     Autority = authority_score(g_new, weights=NA, scale = T)$vector,
                     Density = ecount(g_new)/(vcount(g_new)*(vcount(g_new)-1)),
                     Diameter = diameter(g_new, directed = T),
                     # clust_coeff = triad_census(g)[1], ## check with arsen which value to take
                     # degree_centrailty = centr_degree(g, mode = "total")$centralization,
                     # closness_centrality = centr_clo(g, mode = "total")$centralization, # closeness centrality is not well-defined for disconnected graphs
                     # betweennnes_centrality = centr_betw(g, directed=T, normalized=F)$centralization,
                     Eccentricity_centrality = eccentricity(g_new),
                     Disconnected_num = length(disconnected_node_detector(g_graphnel_new)),
                     Sink_num = length(new_pathway$sink.nodes)
  )
  
  
  g_graphnel_old <- graph::removeNode(names(unlist(graph::nodeData(old_pathway$graph, attr = "kegg.type")))[which(unlist(graph::nodeData(old_pathway$graph, attr = "kegg.type")) == "map")], old_pathway$graph) 
  
  g_old <- igraph.from.graphNEL(g_graphnel_old)
  
  unprocessed_pathway_sinks <- names(unlist(graph::nodeData(g_graphnel_old, attr = "label")))[which(unlist(graph::nodeData(g_graphnel_old, attr = "label")) %in% unlist(graph::nodeData(new_pathway$graph, new_pathway$sink.nodes, attr = "label")))]
  ### take curated sink nodes edit next line
  
  params_old <- list(Shortest_paths = shortest_paths_from_input_to_sink_nodes(old_pathway, sink_nodes = unprocessed_pathway_sinks)[which(!is.infinite(shortest_paths_from_input_to_sink_nodes(old_pathway, sink_nodes = unprocessed_pathway_sinks)), arr.ind = T)], 
                     Degree_in = igraph::degree(g_old, mode = "in"),
                     Degree_out = igraph::degree(g_old, mode = "out"),
                     Closeness_all = closeness(g_old, mode = "all"),
                     Betweenness_all = igraph::betweenness(g_old, directed = T),
                     Hub = hub_score(g_old, weights=NA, scale = T)$vector,
                     Autority = authority_score(g_old, weights=NA, scale = T)$vector,
                     Density = ecount(g_old)/(vcount(g_old)*(vcount(g_old)-1)),
                     Diameter = diameter(g_old, directed = T),
                     # clust_coeff = triad_census(g)[1], ## check with arsen which value to take
                     # degree_centrailty = centr_degree(g, mode = "total")$centralization,
                     # closness_centrality = centr_clo(g, mode = "total")$centralization, # closeness centrality is not well-defined for disconnected graphs
                     # betweennnes_centrality = centr_betw(g, directed=T, normalized=F)$centralization,
                     Eccentricity_centrality = eccentricity(g_old),
                     Disconnected_num = length(disconnected_node_detector(g_graphnel_old)),
                     Sink_num = length(old_pathway$sink.nodes)
  )
  
  return(list(new = params_new, old = params_old))
  
}

### calculation of raw and curated network parameters
params_new <- network_param_calculator(old_pathway = kegg_collection_unporecessed$Hippo_signaling_pathway, new_pathway = edited_pathways_new$Hippo_signaling_pathway)[[1]]
params_old <- network_param_calculator(old_pathway = kegg_collection_unporecessed$Hippo_signaling_pathway, new_pathway = edited_pathways_new$Hippo_signaling_pathway)[[2]]    


shortest_paths_df <- data.frame(group = c(rep("Raw", length(raw_shortest_paths)), rep("Curated", length(curated_shortest_paths))),
                                shortest_path_value = c(raw_shortest_paths, curated_shortest_paths), stringsAsFactors = F
)

betweenness_df <- data.frame(group = c(rep("Raw", length(params_old$Betweenness_all)), rep("Curated", length(params_new$Betweenness_all))),
                             betweenness = c(params_old$Betweenness_all, params_new$Betweenness_all), stringsAsFactors = F
)

closeness_df <- data.frame(group = c(rep("Raw", length(params_old$Closeness_all)), rep("Curated", length(params_new$Closeness_all))),
                           closeness = c(params_old$Closeness_all, params_new$Closeness_all), stringsAsFactors = F
)

eccentricity_df <- data.frame(group = c(rep("Raw", length(params_old$Eccentricity_centrality)), rep("Curated", length(params_new$Eccentricity_centrality))),
                              eccentricity = c(params_old$Eccentricity_centrality, params_new$Eccentricity_centrality), stringsAsFactors = F
)

### creating joint histograms for raw and curated network parameters
plot1 <- shortest_paths_df %>%
  ggplot(aes(x = shortest_path_value, fill = group)) +
  geom_histogram(color = "#e9ecef", alpha = 0.6, position = 'identity', bins = 10) +
  scale_fill_manual(values = c("#69b3a2", "#404080")) +
  theme_ipsum(base_family = "Helvetica", base_size = 14) +
  labs(x = "Shortest Paths", y = "Frequency", fill = "") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.title.x = element_text(size = 18, hjust = 0.5),
        axis.title.y = element_text(size = 18, hjust = 0.5),
        plot.title = element_text(size = 20, hjust = 0.5)) +
  coord_cartesian(clip = "off")

# Represent it
plot2 <- betweenness_df %>%
  ggplot(aes(x = betweenness, fill = group)) +
  geom_histogram(color = "#e9ecef", alpha = 0.6, position = 'identity', bins = 10) +
  scale_fill_manual(values = c("#69b3a2", "#404080")) +
  theme_ipsum(base_family = "Helvetica", base_size = 14) +
  labs(x = "Betweenness Centrality", y = "", fill = "") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, hjust = 0.5),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.title.x = element_text(size = 18, hjust = 0.5),
        axis.title.y = element_text(size = 18, hjust = 0.5),
        plot.title = element_text(size = 20, hjust = 0.5)) +
  coord_cartesian(clip = "off")  # To prevent axis labels from being clipped

# Represent it
plot3 <- closeness_df %>%
  ggplot(aes(x = closeness, fill = group)) +
  geom_histogram(color = "#e9ecef", alpha = 0.6, position = 'identity', bins = 10) +
  scale_fill_manual(values = c("#69b3a2", "#404080")) +
  theme_ipsum(base_family = "Helvetica", base_size = 14) +
  labs(x = "Closeness Centrality", y = "Frequency", fill = "") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, hjust = 0.5),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.title.x = element_text(size = 18, hjust = 0.5),
        axis.title.y = element_text(size = 18, hjust = 0.5),
        plot.title = element_text(size = 20, hjust = 0.5)) +
  coord_cartesian(clip = "off")  # To prevent axis labels from being clipped

# Represent it
plot4 <- eccentricity_df %>%
  ggplot(aes(x = eccentricity, fill = group)) +
  geom_histogram(color = "#e9ecef", alpha = 0.6, position = 'identity', bins = 10) +
  scale_fill_manual(values = c("#69b3a2", "#404080")) +
  theme_ipsum(base_family = "Helvetica", base_size = 14) +
  labs(x = "Eccentricity", y = "", fill = "") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, hjust = 0.5),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.title.x = element_text(size = 18, hjust = 0.5),
        axis.title.y = element_text(size = 18, hjust = 0.5),
        plot.title = element_text(size = 20, hjust = 0.5)) +
  coord_cartesian(clip = "off")  # To prevent axis labels from being clipped

hist_grid <- grid.arrange(plot1, plot2, plot3, plot4, ncol = 2)

ggsave("hippo_paprams_hists.pdf", hist_grid, height = 8.27, width = 11.69)
