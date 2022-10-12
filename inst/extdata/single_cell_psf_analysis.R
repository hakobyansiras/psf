#### Pathway activity analysis on GTEx single cell expression data, comparison of gene and pathway level UMAP clustering and feature analysis 

library(zellkonverter)
library(SingleCellExperiment)
library(Seurat)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(plotly)
library(psf)
library(biomaRt)
### loading curated pathway collection
load(system.file("extdata", "edited_pathways_new.RData", package="psf"))

### Downloading gene symbol to entrez id conversion data
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
gene_symbol_to_entrez <- getBM(attributes=c('entrezgene_id', 'hgnc_symbol'), mart = ensembl)
gene_symbol_to_entrez <- gene_symbol_to_entrez[which(!is.na(gene_symbol_to_entrez$entrezgene_id)),]
gene_symbol_to_entrez <- gene_symbol_to_entrez[which(!duplicated(gene_symbol_to_entrez$entrezgene_id)),]
gene_symbol_to_entrez <- setNames(object = gene_symbol_to_entrez$entrezgene_id, nm = gene_symbol_to_entrez$hgnc_symbol)

### loading GTEx single cell expression data (note that this data set is big and will take a lot of RAM)
gtex_single_cell <- zellkonverter::readH5AD("GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad")

gtex_coldata <- colData(gtex_single_cell)
gtex_coldata$tissue <- as.character(gtex_coldata$tissue)

### umap plotting function
seurat_umap_builder <- function(input_data, metadata, plot_title, group_by = "group", dim_num = 30, approx_state = TRUE, normalize = FALSE) {
  
  cbmc <- CreateSeuratObject(counts = input_data, meta.data = metadata)
  
  # perform visualization and clustering steps
  if(normalize) {
    cbmc <- NormalizeData(cbmc)
  }
  cbmc <- FindVariableFeatures(cbmc)
  cbmc <- ScaleData(cbmc)
  cbmc <- RunPCA(cbmc, verbose = FALSE, approx=approx_state)
  cbmc <- FindNeighbors(cbmc, dims = 1:dim_num)
  cbmc <- FindClusters(cbmc, resolution = 0.8, verbose = FALSE)
  cbmc <- RunUMAP(cbmc, dims = 1:dim_num)
  DimPlot(cbmc, label = TRUE, group.by = group_by) + 
    ggtitle(plot_title)
  
}

### Pseudo-bulk transformation of single cell data to avoid from many 0 expressions 
gtex_single_cell_aggregated_data <- lapply(unique(gtex_coldata$tissue), function(x) {
  
  tissue_samples <- sample(rownames(gtex_coldata)[which(gtex_coldata$tissue == x)], 5000)
  print(x)
  tissue_matrix <- Reduce(cbind,
                          lapply(0:99, function(y) {
                            
                            cell_range <- 1:50 + 50*y
                            
                            rowSums(exp(as.matrix(assay(gtex_single_cell[,tissue_samples[cell_range]]))))
                            
                          })
  )
  
  colnames(tissue_matrix) <- tissue_samples[1:100]
  
  tissue_matrix
})

gtex_single_cell_aggregated_data_matrix <- Reduce(cbind, gtex_single_cell_aggregated_data)

### calculating fold change values against global mean
gtex_single_cell_aggregated_data_fc <- as.matrix(gtex_single_cell_aggregated_data_matrix)/rowMeans(gtex_single_cell_aggregated_data_matrix)
### converting gene symbols to Entrez id for PSF analysis
rownames(gtex_single_cell_aggregated_data_fc) <- gene_symbol_to_entrez[rownames(gtex_single_cell_aggregated_data_fc)]

gtex_single_cell_aggregated_data_fc <- gtex_single_cell_aggregated_data_fc[which(!is.na(rownames(gtex_single_cell_aggregated_data_fc))),]

### Building UMAP plot for pseudo-bulk transformed gene level expression data
seurat_umap_builder(gtex_single_cell_aggregated_data_fc, as.data.frame(gtex_coldata[colnames(gtex_single_cell_aggregated_data_fc),]), plot_title = "Pseudobulk gene level UMAP", group_by = "tissue")


### building UMAP with only selected genes from 35 signaling pathways
genes_of_signaling_pathways <- unique(unlist(lapply(edited_pathways_new, function(x) {unname(unlist(graph::nodeData(x$graph, attr = "genes")))})))

genes_of_signaling_pathways_filtered <- genes_of_signaling_pathways[which(genes_of_signaling_pathways %in% rownames(gtex_single_cell_aggregated_data_fc))]

seurat_umap_builder(gtex_single_cell_aggregated_data_fc[genes_of_signaling_pathways_filtered,], as.data.frame(gtex_coldata[colnames(gtex_single_cell_aggregated_data_fc),]), plot_title = "Gtex subset aggregated", group_by = "tissue")

### Activity calculation with psf algorithm for 35 KEGG signaling pathways on pseudo-bulk transformed single cell fold change data

## psf calculation script (will run around several hours on single core)
gtex_single_cell_aggregated_psf <- psf.from.env.entrez.fc(entrez.fc = gtex_single_cell_aggregated_data_fc,
                                                          kegg.collection = edited_pathways_new, 
                                                          calculate.significance = F, sum = FALSE, return_only_signals = T)

## building psf matrix from pathway sink node psf values
gtex_single_cell_aggregated_psf_matrix <- Reduce(rbind,
                                                 lapply(names(edited_pathways_new), function(x) {
                                                   pathway_mat <- sapply(gtex_single_cell_aggregated_psf, function(y) {
                                                     y[[x]][edited_pathways_new[[x]]$sink.nodes]
                                                   })
                                                   
                                                   rownames(pathway_mat) <- paste(rownames(pathway_mat), x, sep = "_")
                                                   pathway_mat
                                                 })
)

## Building UMAP plot for calculated activity values
seurat_umap_builder(gtex_single_cell_aggregated_psf_matrix, as.data.frame(gtex_coldata[colnames(gtex_single_cell_aggregated_data_fc),]), plot_title = "Pathway activity UMAP", group_by = "tissue")



#### Cluster feature analysis ####
### detecting top cluster features
cbmc <- CreateSeuratObject(counts = gtex_single_cell_aggregated_psf_matrix, meta.data = as.data.frame(gtex_coldata[colnames(gtex_single_cell_aggregated_data_fc),]))

# perform visualization and clustering steps
# if(normalize) {
#   cbmc <- NormalizeData(cbmc)
# }
cbmc <- FindVariableFeatures(cbmc)
cbmc <- ScaleData(cbmc)
cbmc <- RunPCA(cbmc, verbose = FALSE, approx=TRUE)
cbmc <- FindNeighbors(cbmc, dims = 1:30)
cbmc <- FindClusters(cbmc, resolution = 0.8, verbose = FALSE)
cbmc <- RunUMAP(cbmc, dims = 1:30)
cbmc <- SetIdent(cbmc, value = cbmc@meta.data$tissue)

DimPlot(cbmc, label = TRUE) + 
  ggtitle("PSF umap")

### Calculating cluster specific top 20 features from PSF UMAP clusters with Wilcoxon Rank Sum test
clust_top_features <- lapply(unique(gtex_coldata$tissue), function(x) {
  head(FindMarkers(cbmc, ident.1 = x, logfc.threshold = 0.25, test.use = "wilcox", only.pos = FALSE), n = 20)
})

names(clust_top_features) <- unique(gtex_coldata$tissue)

lapply(clust_top_features, function(x) {
  sort(table(sapply(rownames(x), function(y) {
    paste(unlist(strsplit(y, split = "-"))[2:length(unlist(strsplit(y, split = "-")))], collapse = "_")
  })), decreasing = T)
})

### Ploting top cluster specific pathway frequncies fro skeletalmuscle and heart tissues
heart_top_features <- clust_top_features$heart

skeletalmuscle_top_features <- clust_top_features$skeletalmuscle

heart_top_features$pathway <- sapply(rownames(clust_top_features$heart), function(x) {paste(unlist(strsplit(x, split = "-"))[2:length(unlist(strsplit(x, split = "-")))], collapse = "-")})

heart_top_features <- heart_top_features[which(heart_top_features$avg_log2FC > 0),]

skeletalmuscle_top_features$pathway <- sapply(rownames(clust_top_features$skeletalmuscle), function(x) {paste(unlist(strsplit(x, split = "-"))[2:length(unlist(strsplit(x, split = "-")))], collapse = "-")})

skeletalmuscle_top_features <- skeletalmuscle_top_features[which(skeletalmuscle_top_features$avg_log2FC > 0), ]

plot_ly(x = names(table(skeletalmuscle_top_features$pathway)), y = unname(table(skeletalmuscle_top_features$pathway)), type = "bar")

plot_ly(x = names(table(heart_top_features$pathway)), y = unname(table(heart_top_features$pathway)), type = "bar")


### Extracting full list of cluster specific features
clust_full_features <- lapply(unique(gtex_coldata$tissue), function(x) {
  FindMarkers(cbmc, ident.1 = x, logfc.threshold = 0.25, test.use = "wilcox", only.pos = FALSE)
})

names(clust_full_features) <- unique(gtex_coldata$tissue)



### Making heatmaps of pseudo bulk fold change values for FoxO and Calcium signaling pathway genes. 
## FoxO signaling was the most frequnet upragulted pathway in skeletalmuscle tissue cluster
## Calcium signaling was the most frequnet upragulted pathway in heart tissue cluster

cell_annot <- data.frame(row.names = colnames(gtex_single_cell_aggregated_data_fc), category = gtex_coldata[colnames(gtex_single_cell_aggregated_data_fc),"tissue"])

# skeletalmuscle
gtex_pseudobulk_foxO_gene_exp_heatmap <- pheatmap(gtex_single_cell_aggregated_data_fc[intersect(rownames(gtex_single_cell_aggregated_data_fc), unique(unname(unlist(graph::nodeData(edited_pathways_new$FoxO_signaling_pathway$graph, attr = "genes"))))),], 
                                                  cluster_cols = T, cluster_rows = T, annotation_col = cell_annot, show_colnames = F,
                                                  labels_row = as.character(entrez_to_symbol[intersect(rownames(gtex_single_cell_aggregated_data_fc), unique(unname(unlist(graph::nodeData(edited_pathways_new$FoxO_signaling_pathway$graph, attr = "genes"))))),]),
                                                  fontsize_row = 5# filename = "FoxO_signaling_skeletal_muscle.png"
)

### selecting tissue associated upragulated genes
foxO_highligh_nodes <- unlist(lapply(graph::nodes(edited_pathways_new$FoxO_signaling_pathway$graph), function(x) {
  node_genes <- unname(unlist(graph::nodeData(edited_pathways_new$FoxO_signaling_pathway$graph, x,attr = "genes")))
  if(length(node_genes) > 0) {
    if(any(node_genes %in% intersect(rownames(gtex_single_cell_aggregated_data_fc), unique(unname(unlist(graph::nodeData(edited_pathways_new$FoxO_signaling_pathway$graph, attr = "genes")))))[gtex_pseudobulk_foxO_gene_exp_heatmap$tree_row$order[1:11]])) {
      x
    }
  }
}))

## highligthing them on the pathway
psf::plot_kegg_image_pathway(edited_pathways_new$FoxO_signaling_pathway, plot_type = "visnet", highlight_nodes = foxO_highligh_nodes)

# heart
gtex_pseudobulk_calcium_signaling_gene_exp_heatmap <- pheatmap(gtex_single_cell_aggregated_data_fc[intersect(rownames(gtex_single_cell_aggregated_data_fc), unique(unname(unlist(graph::nodeData(edited_pathways_new$Calcium_signaling_pathway$graph, attr = "genes"))))),], 
                                                               cluster_cols = T, cluster_rows = T, annotation_col = cell_annot, show_colnames = F,
                                                               labels_row = as.character(entrez_to_symbol[intersect(rownames(gtex_single_cell_aggregated_data_fc), unique(unname(unlist(graph::nodeData(edited_pathways_new$Calcium_signaling_pathway$graph, attr = "genes"))))),]),
                                                               fontsize_row = 5
)

### selecting tissue associated upragulated genes
calcium_signaling_highligh_nodes <- unlist(lapply(graph::nodes(edited_pathways_new$Calcium_signaling_pathway$graph), function(x) {
  node_genes <- unname(unlist(graph::nodeData(edited_pathways_new$Calcium_signaling_pathway$graph, x,attr = "genes")))
  if(length(node_genes) > 0) {
    if(any(node_genes %in% intersect(rownames(gtex_single_cell_aggregated_data_fc), unique(unname(unlist(graph::nodeData(edited_pathways_new$Calcium_signaling_pathway$graph, attr = "genes")))))[gtex_pseudobulk_calcium_signaling_gene_exp_heatmap$tree_row$order[1:9]])) {
      x
    }
  }
}))

## highligthing them on the pathway
psf::plot_kegg_image_pathway(edited_pathways_new$Calcium_signaling_pathway, plot_type = "visnet", highlight_nodes = calcium_signaling_highligh_nodes)



### making pathway plots with marked upragulated and downregulated sink nodes(terminal nodes) with red and blue colors
dir.create("umap_psf_cluster_pathway_plots")
lapply(names(clust_top_features), function(x) {
  
  pathway_sinks <- gsub("-", "_", rownames(clust_top_features[[x]]))
  
  tissue_pathways <- unique(
    sapply(pathway_sinks, function(y) {
      paste(unlist(strsplit(y, split = "_"))[2:length(unlist(strsplit(y, split = "_")))], collapse = "_")
    })
  )
  
  pathway_plots <- lapply(tissue_pathways, function(z) {
    psf::plot_kegg_image_pathway(edited_pathways_new[[z]], plot_type = "kegg", use_old_images = T,
                                 highlight_nodes = gsub(paste0("_", z), "", pathway_sinks[grep(z, pathway_sinks)]),
                                 highlight_color = ifelse(clust_top_features[[x]]$avg_log2FC[grep(z, pathway_sinks)] > 0, "red", "blue") 
    )
  })
  
  magick::image_write(Reduce(c, pathway_plots), format = "pdf", path = paste0("umap_psf_cluster_pathway_plots/", x, ".pdf"))
  
  # dev.off()
})



#### compariosn between gene cluster features enrichment and PSF cluster features ####
### detecting top cluster features
gene_cbmc_pseudobulk <- CreateSeuratObject(counts = gtex_single_cell_aggregated_data_fc, meta.data = as.data.frame(gtex_coldata[colnames(gtex_single_cell_aggregated_data_fc),]))

# perform visualization and clustering steps
# if(normalize) {
#   gene_cbmc_pseudobulk <- NormalizeData(gene_cbmc_pseudobulk)
# }
gene_cbmc_pseudobulk <- FindVariableFeatures(gene_cbmc_pseudobulk)
gene_cbmc_pseudobulk <- ScaleData(gene_cbmc_pseudobulk)
gene_cbmc_pseudobulk <- RunPCA(gene_cbmc_pseudobulk, verbose = FALSE, approx=TRUE)
gene_cbmc_pseudobulk <- FindNeighbors(gene_cbmc_pseudobulk, dims = 1:30)
gene_cbmc_pseudobulk <- FindClusters(gene_cbmc_pseudobulk, resolution = 0.8, verbose = FALSE)
gene_cbmc_pseudobulk <- RunUMAP(gene_cbmc_pseudobulk, dims = 1:30)
gene_cbmc_pseudobulk <- SetIdent(gene_cbmc_pseudobulk, value = gene_cbmc_pseudobulk@meta.data$tissue)

DimPlot(gene_cbmc_pseudobulk, label = TRUE) + 
  ggtitle("test")


### extracting cluster specific features(genes)
gene_clust_full_features_pseudobul <- lapply(unique(gtex_coldata$tissue), function(x) {
  FindMarkers(gene_cbmc_pseudobulk, ident.1 = x, logfc.threshold = 0.25, test.use = "wilcox", only.pos = FALSE)
})

names(gene_clust_full_features_pseudobul) <- unique(gtex_coldata$tissue)


## Enrichment of cluster specific features with WebGestaltR on kegg pathways
tissue_cluster_genestes_pseudobul <- lapply(names(gene_clust_full_features_pseudobul), function(x) {
  
  sig_genes <- rownames(gene_clust_full_features_pseudobul[[x]])[which(gene_clust_full_features_pseudobul[[x]]$p_val_adj < 0.05)]
  
  webgestalt_sets <- WebGestaltR::WebGestaltR(enrichMethod = "ORA", organism = "hsapiens", 
                                              enrichDatabase = "pathway_KEGG", 
                                              interestGeneType = "entrezgene", isOutput = F, 
                                              referenceSet = "genome", 
                                              interestGene = sig_genes)
})

names(tissue_cluster_genestes_pseudobul) <- names(gene_clust_full_features_pseudobul)


### comparing gene level enrichment and pathway level results
pathway_clust_top_features_ids <- lapply(clust_full_features, function(x) {
  pathway <- rownames(x)[which(x$p_val_adj < 0.05)]
  
  unique(sapply(pathway, function(y) {
    pathway_name <- paste(unlist(strsplit(y, split = "-"))[2:length(unlist(strsplit(y, split = "-")))], collapse = "_")
    unlist(strsplit(edited_pathways_new[[pathway_name]]$attrs$name, split = ":"))[2]
  }))
})


## overlap of ORA enriched and PSF deregulated KEGG pathways
pathway_overlap <- sapply(names(tissue_cluster_genestes_pseudobul), function(x) {
  sum(pathway_clust_top_features_ids[[x]] %in% tissue_cluster_genestes_pseudobul[[x]]$geneSet)/length(pathway_clust_top_features_ids[[x]])
})

## ploting proportion of PSF detected pathways also detected with the ORA method. 35 KEGG pathways were compared
plot_ly(
  x = names(pathway_overlap),
  y = unname(pathway_overlap),
  name = "Proportions",
  type = "bar"
)  %>% layout(yaxis = list(title = 'Proportion', tickfont = list(size = 15)),
              xaxis = list(tickfont = list(size = 15))
)

