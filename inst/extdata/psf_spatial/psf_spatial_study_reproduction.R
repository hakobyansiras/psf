library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(biomaRt)
library(psf)
suppressMessages(library(vesalius))
library(gsdensity)
library(ggrepel)
library(reshape2)
library(CelliD)



# Load curated KEGG signaling pathways
load(system.file("extdata", "kegg_curated_40_signalings.RData", package="psf"))
load(system.file("extdata", "kegg_sink_to_process.RData", package="psf"))


#### Data import, Seurat processing, PSF activity calcualtion ####
# Load the human melanoma spatial dataset. You can download the files from 10x genomics database
spatial_melanoma <- Load10X_Spatial("/path/to/spatial_dataset", 
                                    filename = "CytAssist_FFPE_Human_Skin_Melanoma_filtered_feature_bc_matrix.h5")

# Load pre-downloaded gene symbol to Entrez ID conversion data
load("gene_symbol_to_entrez.RData")

# Run pathway analysis
spatial_melanoma_psf <- spatial_psf_analysis(spatial_obj = spatial_melanoma, 
                                    pathway_collection = kegg_curated_40_signalings, 
                                    gene_symbol_to_entrez = gene_symbol_to_entrez, nthreads = 30)

mel_spatial_psf_mat <- as.matrix(spatial_melanoma_psf$spatial_psf_obj[["Spatial"]]$counts)

spatial_mouse_brain <- Load10X_Spatial("mouse_brain_spatial_data/", 
                                       filename = "V1_Mouse_Brain_Sagittal_Anterior_filtered_feature_bc_matrix.h5")

# Load mouse gene symbol to human gene conversion data
load(system.file("extdata", "mouse_symbol_to_human_entrez.RData", package="psf"))

spatial_mouse_brain_psf <- spatial_psf_analysis(spatial_obj = spatial_mouse_brain, pathway_collection = kegg_curated_40_signalings, 
                                                gene_symbol_to_entrez = mouse_symbol_to_human_entrez, 
                                                nthreads = 30, return_only_shiny_vars = FALSE)

mouse_brain_spatial_psf_mat <- as.matrix(spatial_mouse_brain_psf$spatial_psf_obj[["Spatial"]]$counts)


#### Vesalius clustering ####
vesalius_analysis <- function(counts, coords, nfeatures = 2000) {
  vesalius <- build_vesalius_assay(
    coordinates = coords, # spatial coordinates
    counts  = counts, # count matrix
    assay = "spatial_omics", # name you wish to give your assay
    verbose = TRUE # Do you want progress messages?
  )
  
  # mel_vesalius <- generate_embeddings(mel_vesalius,
  #                                 dim_reduction = "PCA",
  #                                 normalization = "log_norm",
  #                                 nfeatures = 709, # Setting number of features low for low run time
  #                                 verbose = FALSE)
  
  
  vesalius <- generate_embeddings(vesalius,
                                  dim_reduction = "UMAP",
                                  nfeatures = nfeatures, # Setting number of features low for low run time
                                  verbose = FALSE)
  
  vesalius <- regularise_image(vesalius, lambda = 1)
  vesalius <- smooth_image(vesalius, sigma = 5, iter = 10)
  vesalius <- equalize_image(vesalius, sleft = 5, sright = 5)
  
  
  vesalius <- segment_image(vesalius,
                            method = "kmeans",
                            col_resolution = 10,
                            verbose = FALSE)
  
  
  vesalius <- isolate_territories(vesalius, capture_radius = 0.05)
  spatial_plot <- territory_plot(vesalius, cex_pt = 3.5, split = F)
  
  
  vesalius <- identify_markers(vesalius, seed = NULL, query = NULL)
  deg <- get_markers(vesalius)
  
  
  return(list(vesalius = vesalius, clust_plot = spatial_plot, deg = deg))
}


mel_coords <- GetTissueCoordinates(object = spatial_melanoma_psf$spatial_psf_obj, scale = "lowres")
mel_coords$barcodes <- rownames(mel_coords)
colnames(mel_coords)[1:2] <- c("x", "y")
mel_coords <- mel_coords[,c(3,1,2)]

mel_vesalius_psf <- vesalius_analysis(counts = as.matrix(spatial_melanoma_psf$spatial_psf_obj@assays$Spatial@counts), coords = mel_coords, nfeatures = 709)
mel_vesalius_exp <- vesalius_analysis(counts = as.matrix(spatial_melanoma_psf$spatial_obj@assays$Spatial@counts), coords = mel_coords)



mouse_brain_coords <- GetTissueCoordinates(object = spatial_mouse_brain_psf$spatial_obj, scale = "lowres")
mouse_brain_coords$barcodes <- rownames(mouse_brain_coords)
colnames(mouse_brain_coords)[1:2] <- c("x", "y")
mouse_brain_coords <- mouse_brain_coords[,c(3,1,2)]

mouse_brain_vesalius_psf <- vesalius_analysis(counts = as.matrix(spatial_mouse_brain_psf$spatial_psf_obj@assays$Spatial@counts), coords = mouse_brain_coords, nfeatures = 709)
mouse_brain_vesalius_exp <- vesalius_analysis(counts = as.matrix(spatial_mouse_brain_psf$spatial_obj@assays$Spatial@counts), coords = mouse_brain_coords)


#### spatialGE clustering, spatially relevant pathway detection ####
## Seurat data preprocessing function for sptailGE analyis
seurat_spatial_spatialGE_processing <- function(seurat) {
  counts_df = as.data.frame(seurat[["Spatial"]]$counts)
  counts_df[['gene']] = rownames(counts_df)
  counts_df = counts_df[, c( 'gene', grep('gene', colnames(counts_df), value=T, invert=T) ) ]
  rownames(counts_df) = NULL
  
  # Create data frame from coordinate data
  # Add column with spot names (i.e., barcodes)
  # Put barcodes as first column
  # Remove row names
  coords_df = GetTissueCoordinates(seurat)
  coords_df[['barcode']] = rownames(coords_df)
  coords_df = coords_df[, c('barcode', 'x', 'y')]
  colnames(coords_df) = c("barcode", "imagerow", "imagecol")
  rownames(coords_df) = NULL
  
  # Make named list with count and coordinate data frames
  # Names should match
  counts_list = list(sample1=counts_df)
  coords_list = list(sample1=coords_df)
  
  # Create STlist
  STlist(rnacounts=counts_list, spotcoords=coords_list, samples=NULL) 
}


## SpatialGE pathway detection and clustering
spatialGE_analysis <- function(seurat, seurat_psf, gene.set.list) {
  # Create STlist
  st_list <- seurat_spatial_spatialGE_processing(seurat = seurat) 
  st_list <- transform_data(st_list, method='sct')
  
  stenrich_df <- STenrich(st_list, gene_sets=gene.set.list, 
                          reps=1000, 
                          num_sds=1.5, 
                          min_genes=5,
                          min_units=10, 
                          seed=12345)
  
  st_list <- STclust(st_list, ws=0.02, ks='dtc')
  
  
  
  psf_st_list <- seurat_spatial_spatialGE_processing(seurat = seurat_psf) 
  
  psf_st_list <- transform_data(psf_st_list, method='log')
  
  psf_st_list <- STclust(psf_st_list, ws=0.02, ks='dtc')
  
  return(list(stenrich_df = stenrich_df, exp_metada = st_list, psf_metadata = psf_st_list))
}


### Creating human gene sets based on 40 kegg signaling pathways
gene.set.list <- lapply(kegg_curated_40_signalings, function(x) {
  as.character(psf:::entrez_to_symbol[as.character(unique(unlist(graph::nodeData(x$graph, attr = "genes")))),])
})

## running spatialGE on human melanoma dataset
mel_spatilaGE <- spatialGE_analysis(seurat = spatial_melanoma_psf$spatial_obj, seurat_psf = spatial_melanoma_psf$spatial_psf_obj,
                                    gene.set.list = gene.set.list)


### Creating mouse gene sets based on 40 kegg signaling pathways
human_entrez_to_mouse_symbol <- setNames(object = names(mouse_symbol_to_human_entrez), nm = unname(mouse_symbol_to_human_entrez))
gene.set.list_mouse <- lapply(kegg_curated_40_signalings, function(x) {
  as.character(human_entrez_to_mouse_symbol[as.character(unique(unlist(graph::nodeData(x$graph, attr = "genes"))))])
})

## running spatialGE on mouse brain dataset
mouse_brain_spatialGE <- spatialGE_analysis(seurat = spatial_mouse_brain_psf$spatial_obj, seurat_psf = spatial_mouse_brain_psf$spatial_psf_obj,
                                            gene.set.list = gene.set.list_mouse)



#### GSDensity spatially relevant pathway detection ####
GSDesnity_pipline <- function(seurat_obj, gene.set.list, psf_matrix) {
  
  seurat_obj <- SCTransform(seurat_obj, assay = "Spatial", verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, assay = "SCT", verbose = FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:30)
  seurat_obj <- FindClusters(seurat_obj, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:30)
  
  
  ce <- compute.mca(object = seurat_obj, assay = "Spatial", slot = "Spatial")
  
  
  res <- computekldFix(coembed = ce, genes.use = intersect(rownames(ce),
                                                           rownames(seurat_obj)), n.grids = 100, gene.set.list = gene.set.list,
                       gene.set.cutoff = 20, n.times = 100)
  
  
  gene.set.deviated <- res[res$p.adj < 0.05, ]$gene.set
  
  
  # compute a nearest neighbor graph (edge list) in the MCA
  # space
  cells <- colnames(seurat_obj)
  el <- compute.nn.edges(coembed = ce, nn.use = 300)
  
  # We then compute the relevance between each cell and the
  # deviated gene sets
  
  cv.df <- run.rwr.list(el = el, gene_set_list = gene.set.list[gene.set.deviated],
                        cells = cells)
  
  cl.df <- compute.cell.label.df(cv.df)
  
  
  # An optional filtering step: we want to only keep the
  # terms with certain numbers of positive cells; here we use
  # 100
  
  positive.count <- apply(cl.df, MARGIN = 2, FUN = function(x) {
    length(x[x == "positive"])
  })
  print(positive.count)
  gene.set.deviated.2 <- names(positive.count[positive.count > 80])
  
  
  
  coords.df <- seurat_obj@images$slice1@coordinates[, c("imagerow",
                                                        "imagecol")]
  
  # compute the spatial relevance of gene sets
  
  # the 'weight_df' should have a format as the output of
  # 'run.rwr.list'; here we use the terms with at least 100
  # positive cells the parameter 'n' defines how to split the
  # spatial map. n = 10 means that 10 splits are made in each
  # dimension (total 100 grids) for the kde process
  
  spatial.klds <- compute.spatial.kld.df(spatial.coords = coords.df,
                                         weight_df = cv.df[, gene.set.deviated.2], n = 10)
  
  # Then we want to nominate gene sets: here we want to find
  # highly spatially related gene sets
  top.spatial.terms <- rownames(spatial.klds[spatial.klds$spatial.kld >
                                               quantile(spatial.klds$spatial.kld, 0.8), ])
  
  
  
  
  #### Spatial KLD analysis on PSF activity values ####
  psf_spatial.klds <- compute.spatial.kld.df(spatial.coords = coords.df,
                                             weight_df = t(psf_matrix), n = 10)
  
  # Then we want to nominate gene sets: here we want to find
  # highly spatially related gene sets
  psf_top.spatial.branches <- rownames(psf_spatial.klds[psf_spatial.klds$spatial.kld >
                                                          quantile(psf_spatial.klds$spatial.kld, 0.8, na.rm = T), ])
  
  
  # patchwork::wrap_plots(lapply(gene.set.deviated.2, function(x) {
  #   seurat_obj@meta.data[,x] <- cv.df[rownames(seurat_obj@meta.data), x]
  # 
  #   SpatialFeaturePlot(seurat_obj, features = x) +
  #     theme(legend.position = "top")
  # }))
  
  return(list(GSDensity_spatial_pathways = top.spatial.terms, psf_spatial_pathways = psf_top.spatial.branches,
              cv.df = cv.df, gene.set.deviated.2 = gene.set.deviated.2,
              exp_spatial.klds = spatial.klds, psf_spatial.klds = psf_spatial.klds
  ))
}

mel_GSdensity_results <- GSDesnity_pipline(seurat_obj = spatial_melanoma, gene.set.list = gene.set.list, 
                                           psf_matrix = mel_spatial_psf_mat)

mouse_brain_GSdensity_results <- GSDesnity_pipline(seurat_obj = spatial_mouse_brain, gene.set.list = gene.set.list_mouse, 
                                                   psf_matrix = mouse_brain_spatial_psf_mat)




#### Cluster specific pathway detection and comparison ####
cluster_based_pathway_comp <- function(idents_list, seurat_list, gene.set.list, cv.df, dir_name, gene.set.deviated.2 = gene.set.deviated.2, output_figures = T, gene_marker_fc = 0.5, psf_marker_fc = 0.25, GSDensity_quantile = 0.6) {
  
  if(output_figures) {
    dir.create(dir_name)
  }
  
  per_cluster_altered_pathways <- lapply(names(idents_list), function(x) {
    
    exp_seurat <- seurat_list$spatial_obj
    psf_seurat <- seurat_list$spatial_psf_obj
    Idents(exp_seurat) <- idents_list[[x]]
    Idents(psf_seurat) <- idents_list[[x]]
    
    altered_pathways <- lapply(as.character(sort(unique(as.numeric(as.character(Idents(exp_seurat)))))), function(y) {
      gene_markers <- FindMarkers(exp_seurat, ident.1 = y, logfc.threshold = gene_marker_fc, test.use = "wilcox")
      ora_res <- fgsea::fora(pathways = gene.set.list, genes = rownames(gene_markers)[which(gene_markers$p_val_adj < 0.05)] , universe = rownames(exp_seurat))
      ORA_pathways <- ora_res[which(ora_res$padj < 0.05),]$pathway
      
      
      psf_markers <- FindMarkers(psf_seurat, ident.1 = y, logfc.threshold = 0, test.use = "wilcox")
      psf_markers <- psf_markers[which(psf_markers$p_val_adj < 0.05),]
      psf_markers$Pathway_name <- sapply(rownames(psf_markers), function(z) {gsub("-", "_", unlist(strsplit(z, split = "; "))[2])})
      
      per_pathway_total_FC <- sapply(unique(psf_markers$Pathway_name), function(p) {
        sum(abs(psf_markers$avg_log2FC[which(psf_markers$Pathway_name == p)]))
      })
      
      
      psf_markers <- psf_markers[which(abs(psf_markers$avg_log2FC) > psf_marker_fc),]
      
      
      if(nrow(psf_markers) > 0) {
        altered_pathway_branches <- table(psf_markers$Pathway_name)
        
        altered_pathway_branch_ratios <- data.frame(pathway = names(altered_pathway_branches), 
                                                    n_sinks = pathwya_sink_nums[names(altered_pathway_branches)],
                                                    n_altered_sinks = as.integer(unname(altered_pathway_branches)),
                                                    ratio = as.numeric(unname(altered_pathway_branches)/pathwya_sink_nums[names(altered_pathway_branches)]),
                                                    PSF_FC_sum = per_pathway_total_FC[names(altered_pathway_branches)]
        )
      }
      
      
      PSF_pathways <- unique(psf_markers$Pathway_name)
      # if (length(PSF_pathways) == 0) {
      #   PSF_pathways <- character()
      # }
      
      
      # psf_mean_markers <- FindMarkers(mel_sp_mean_psf_obj, ident.1 = y, logfc.threshold = 0.25, test.use = "wilcox")
      # PSF_mean_pathways <- gsub("-", "_", rownames(psf_mean_markers)[which(psf_mean_markers$p_val_adj < 0.05)])
      # 
      # psf_sum_markers <- FindMarkers(mel_sp_sum_psf_obj, ident.1 = y, logfc.threshold = 0.25, test.use = "wilcox")
      # PSF_sum_pathways <- gsub("-", "_", rownames(psf_sum_markers)[which(psf_sum_markers$p_val_adj < 0.05)])
      
      
      sp_metadata <- exp_seurat@meta.data
      sp_metadata$current_idnets <- idents_list[[x]][rownames(sp_metadata)]
      
      jsd.df <- compute.spec(cell_df = cv.df[, gene.set.deviated.2], 
                             metadata = sp_metadata, # each row is a cell; columns include partition information
                             cell_group = "current_idnets" # 'cell_group' should use a column name in the metadata as the input
      )
      
      
      sorted_vals <- sort(jsd.df[,y])
      GSdensity_pathways <- names(which(sorted_vals > quantile(sorted_vals, GSDensity_quantile)))
      
      if(nrow(psf_markers) > 0) {
        altered_pathway_branch_ratios <- cbind(altered_pathway_branch_ratios, 
                                               ORA = altered_pathway_branch_ratios$pathway %in% ORA_pathways, 
                                               # PSF_mean = altered_pathway_branch_ratios$pathway %in% PSF_mean_pathways, 
                                               # PSF_sum = altered_pathway_branch_ratios$pathway %in% PSF_sum_pathways, 
                                               GSdensity = altered_pathway_branch_ratios$pathway %in% GSdensity_pathways, 
                                               any_overlap_with_ORA_GSdensity = altered_pathway_branch_ratios$pathway %in% unique(c(ORA_pathways, GSdensity_pathways)))
      } else {
        altered_pathway_branch_ratios <- NULL
      }
      
      list(gene_markers = gene_markers, ora_res = ora_res, psf_markers = psf_markers, jsd.df = jsd.df, 
           ORA = ORA_pathways, PSF = PSF_pathways,
           # PSF_mean = PSF_mean_pathways, PSF_sum = PSF_sum_pathways, 
           GSdensity = GSdensity_pathways, altered_pathway_branch_ratios = altered_pathway_branch_ratios)
    })
    names(altered_pathways) <- paste0(x, "_", as.character(sort(unique(as.numeric(as.character(Idents(exp_seurat)))))))
    
    if(output_figures) {
      ggsave(SpatialDimPlot(psf_seurat, label = TRUE, label.size = 6, pt.size.factor = 1.5, stroke = 0) + ggtitle(x), filename = paste0(dir_name, "/", x, "_umap.pdf"))
      
      
      pdf(file = paste0(dir_name, "/", x, "_venns.pdf"), paper='A4r', height = 8.27, width = 11.69, onefile = T)
      
      layout(matrix(1:16, ncol = 4, nrow = 4, byrow = T))
      lapply(names(altered_pathways), function(z) {
        venn::venn(altered_pathways[[z]][c("PSF", "ORA", "GSdensity")],
                   zcolor = "style", box = F, ilabels = "counts", lty = 1, ilcs = 1.3, sncs = 1.5, lwd = 1.3)
        title(main = z, line = -0.5, cex.main = 1.5)
      })
      
      dev.off()
      
      qpdf::pdf_combine(
        input = c(paste0(dir_name, "/", x, "_umap.pdf"), paste0(dir_name, "/", x, "_venns.pdf")),
        output = paste0(dir_name, "/", x, ".pdf")
      )
      
      file.remove(c(paste0(dir_name, "/", x, "_umap.pdf"), paste0(dir_name, "/", x, "_venns.pdf")))
    }
    
    altered_pathways
    
  })
  names(per_cluster_altered_pathways) <- names(idents_list)
  return(per_cluster_altered_pathways)
}

mel_vesalius_exp_markers <- setNames(object = mel_vesalius_exp$vesalius@territories$Territory, nm = mel_vesalius_exp$vesalius@territories$barcodes)
mel_vesalius_exp_markers[which(mel_vesalius_exp_markers == "isolated")] <- max(as.integer(mel_vesalius_exp_markers), na.rm = T) + 1
mel_vesalius_exp_markers <- setNames(nm = names(mel_vesalius_exp_markers), object = factor(as.integer(mel_vesalius_exp_markers)))

mel_vesalius_psf_markers <- setNames(object = mel_vesalius_psf$vesalius@territories$Territory, nm = mel_vesalius_psf$vesalius@territories$barcodes)
mel_vesalius_psf_markers[which(mel_vesalius_psf_markers == "isolated")] <- max(as.integer(mel_vesalius_psf_markers), na.rm = T) + 1
mel_vesalius_psf_markers <- setNames(nm = names(mel_vesalius_psf_markers), object = factor(as.integer(mel_vesalius_psf_markers)))


melanoma_sp_idents <- list(
  seurat_exp = Idents(spatial_melanoma_psf$spatial_obj),
  seurat_psf = Idents(spatial_melanoma_psf$spatial_psf_obj),
  spatialGE_exp = setNames(object = mel_spatilaGE$exp_metada@spatial_meta$sample1$stclust_spw0.02_dsplFalse, 
                           nm = mel_spatilaGE$exp_metada@spatial_meta$sample1$libname),
  spatialGE_psf = setNames(object = mel_spatilaGE$psf_metadata@spatial_meta$sample1$stclust_spw0.02_dsplFalse, 
                           nm = mel_spatilaGE$psf_metadata@spatial_meta$sample1$libname),
  vesalius_exp = mel_vesalius_exp_markers,
  vesalius_psf = mel_vesalius_psf_markers
)

mel_sp_cluster_specific_pathways_low <- cluster_based_pathway_comp(idents_list = melanoma_sp_idents, seurat_list = spatial_melanoma_psf, gene.set.list = gene.set.list,
                                                               cv.df = mel_GSdensity_results$cv.df, dir_name = "mel_overlap_vis_plots/", gene.set.deviated.2 = mel_GSdensity_results$gene.set.deviated.2, 
                                                               output_figures = T, gene_marker_fc = 0.5, psf_marker_fc = 0.25, GSDensity_quantile = 0.6)

mel_sp_cluster_specific_pathways_high <- cluster_based_pathway_comp(idents_list = melanoma_sp_idents, seurat_list = spatial_melanoma_psf, gene.set.list = gene.set.list,
                                                                   cv.df = mel_GSdensity_results$cv.df, dir_name = "mel_overlap_vis_plots/", gene.set.deviated.2 = mel_GSdensity_results$gene.set.deviated.2, 
                                                                   output_figures = T, gene_marker_fc = 1, psf_marker_fc = 0.5, GSDensity_quantile = 0.8)


mouse_brain_vesalius_exp_markers <- setNames(object = mouse_brain_vesalius_exp$vesalius@territories$Territory, nm = mouse_brain_vesalius_exp$vesalius@territories$barcodes)
mouse_brain_vesalius_exp_markers[which(mouse_brain_vesalius_exp_markers == "isolated")] <- max(as.integer(mouse_brain_vesalius_exp_markers), na.rm = T) + 1
mouse_brain_vesalius_exp_markers <- setNames(nm = names(mouse_brain_vesalius_exp_markers), object = factor(as.integer(mouse_brain_vesalius_exp_markers)))

mouse_brain_vesalius_psf_markers <- setNames(object = mouse_brain_vesalius_psf$vesalius@territories$Territory, nm = mouse_brain_vesalius_psf$vesalius@territories$barcodes)
mouse_brain_vesalius_psf_markers[which(mouse_brain_vesalius_psf_markers == "isolated")] <- max(as.integer(mouse_brain_vesalius_psf_markers), na.rm = T) + 1
mouse_brain_vesalius_psf_markers <- setNames(nm = names(mouse_brain_vesalius_psf_markers), object = factor(as.integer(mouse_brain_vesalius_psf_markers)))


mouse_brain_sp_idents <- list(
  seurat_exp = Idents(spatial_mouse_brain_psf$spatial_obj),
  seurat_psf = Idents(spatial_mouse_brain_psf$spatial_psf_obj),
  spatialGE_exp = setNames(object = mouse_brain_spatialGE$exp_metada@spatial_meta$sample1$stclust_spw0.02_dsplFalse, 
                           nm = mouse_brain_spatialGE$exp_metada@spatial_meta$sample1$libname),
  spatialGE_psf = setNames(object = mouse_brain_spatialGE$psf_metadata@spatial_meta$sample1$stclust_spw0.02_dsplFalse, 
                           nm = mouse_brain_spatialGE$psf_metadata@spatial_meta$sample1$libname),
  vesalius_exp = mouse_brain_vesalius_exp_markers,
  vesalius_psf = mouse_brain_vesalius_psf_markers
)


mouse_brain_sp_cluster_specific_pathways_low <- cluster_based_pathway_comp(idents_list = mouse_brain_sp_idents, seurat_list = spatial_mouse_brain_psf, gene.set.list = gene.set.list_mouse,
                                                                       cv.df = mouse_brain_GSdensity_results$cv.df, dir_name = "mouse_brain_overlap_vis_plots", gene.set.deviated.2 = mouse_brain_GSdensity_results$gene.set.deviated.2, 
                                                                       output_figures = T, gene_marker_fc = 0.5, psf_marker_fc = 0.25, GSDensity_quantile = 0.6)

mouse_brain_sp_cluster_specific_pathways_high <- cluster_based_pathway_comp(idents_list = mouse_brain_sp_idents, seurat_list = spatial_mouse_brain_psf, gene.set.list = gene.set.list_mouse,
                                                                           cv.df = mouse_brain_GSdensity_results$cv.df, dir_name = "mouse_brain_overlap_vis_plots", gene.set.deviated.2 = mouse_brain_GSdensity_results$gene.set.deviated.2, 
                                                                           output_figures = T, gene_marker_fc = 1, psf_marker_fc = 0.5, GSDensity_quantile = 0.8)



pathway_num_boxplot <- function(clust_specific_pathway_list, gene.set.list, seurat_obj, plot_name = "") {
  
  n_pathways_by_methods <- Reduce(rbind,
                                  lapply(names(clust_specific_pathway_list), function(x) {
                                    Reduce(rbind,
                                           lapply(names(clust_specific_pathway_list[[x]]), function(z) {
                                             
                                             y <- clust_specific_pathway_list[[x]][[z]]
                                             
                                             sorted_vals <- sort(y$jsd.df[,unlist(strsplit(z, split = "_"))[3]])
                                             GSdensity_pathways <- names(which(sorted_vals > quantile(sorted_vals, 0.6)))
                                             
                                             low_stringency <- data.frame(
                                               N_Genesets = c(
                                                 length(y$PSF),
                                                 length(GSdensity_pathways),
                                                 length(y$ORA),
                                                 length(intersect(y$PSF, GSdensity_pathways)),
                                                 length(intersect(y$PSF, y$ORA)),
                                                 length(intersect(GSdensity_pathways, y$ORA))
                                               ),
                                               Functional_method = factor(c("PSF", "GSdensity", "ORA", "PSF_GSdensity", "PSF_ORA", "GSdensity_ORA"), 
                                                                          levels = c("PSF", "GSdensity", "ORA", "PSF_GSdensity", "PSF_ORA", "GSdensity_ORA")),
                                               Clustering_method = gsub("_", " ", x), threshold_stringency = "low"
                                             )
                                             
                                             
                                             sorted_vals <- sort(y$jsd.df[,unlist(strsplit(z, split = "_"))[3]])
                                             GSdensity_pathways <- names(which(sorted_vals > quantile(sorted_vals, 0.8)))
                                             
                                             gene_markers <- y$gene_markers[which(abs(y$gene_markers$avg_log2FC) > 1),]
                                             ora_res <- fgsea::fora(pathways = gene.set.list, genes = rownames(gene_markers)[which(gene_markers$p_val_adj < 0.05)] , universe = rownames(seurat_obj))
                                             ORA_pathways <- ora_res[which(ora_res$padj < 0.05),]$pathway
                                             
                                             PSF_pathways <- unique(y$psf_markers$Pathway_name[which(abs(y$psf_markers$avg_log2FC) > 0.5)])
                                             
                                             
                                             high_stringency <- data.frame(
                                               N_Genesets = c(
                                                 length(PSF_pathways),
                                                 length(GSdensity_pathways),
                                                 length(ORA_pathways),
                                                 length(intersect(PSF_pathways, GSdensity_pathways)),
                                                 length(intersect(PSF_pathways, ORA_pathways)),
                                                 length(intersect(GSdensity_pathways, ORA_pathways))
                                               ),
                                               Functional_method = factor(c("PSF", "GSdensity", "ORA", "PSF_GSdensity", "PSF_ORA", "GSdensity_ORA"), 
                                                                          levels = c("PSF", "GSdensity", "ORA", "PSF_GSdensity", "PSF_ORA", "GSdensity_ORA")),
                                               Clustering_method = gsub("_", " ", x), threshold_stringency = "high"
                                             )
                                             
                                             rbind(low_stringency,high_stringency)
                                             
                                           }))
                                  }))
  
  n_pathways_by_methods$threshold_stringency <- factor(n_pathways_by_methods$threshold_stringency, levels = c("low", "high"))
  n_pathways_by_methods$Clustering_method <- gsub("vesalius", "Vesalius", n_pathways_by_methods$Clustering_method)
  n_pathways_by_methods$Clustering_method <- gsub("seurat", "Seurat", n_pathways_by_methods$Clustering_method)
  n_pathways_by_methods$Clustering_method <- factor(n_pathways_by_methods$Clustering_method, levels = c("Seurat exp", "Seurat psf", "spatialGE exp", "spatialGE psf",
                                                                                                        "Vesalius exp", "Vesalius psf"
  ))
  library(reshape2)
  library(ggpubr)
  
  ggboxplot(n_pathways_by_methods, x = "Functional_method", y = "N_Genesets",
            color = "threshold_stringency") +
    labs(title = "", x = "Methods", y = "N sig pathways", fill = "threshold_stringency") +
    theme_minimal() +
    facet_wrap(~Clustering_method, nrow = 3) +
    theme(axis.text.x = element_text(size = 14, angle = 45, vjust=1, hjust=1), axis.title.x = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"),
          legend.text=element_text(size=16), legend.title = element_text(size=16), strip.text = element_text(size=20),
          plot.title = element_text(size = 22, face = "bold")) + ggtitle(plot_name)
}


pathway_num_boxplot(clust_specific_pathway_list = mel_sp_cluster_specific_pathways_low, gene.set.list = gene.set.list, 
                    seurat_obj = spatial_melanoma_psf$spatial_obj, plot_name = "Human melanoma")


pathway_num_boxplot(clust_specific_pathway_list = mouse_brain_sp_cluster_specific_pathways_low, gene.set.list = gene.set.list_mouse, 
                    seurat_obj = spatial_mouse_brain_psf$spatial_obj, plot_name =  "Mouse brain")




#### Plotting clusterings on spatial tissue maps ####
ggsave(patchwork::wrap_plots(
  lapply(names(melanoma_sp_idents)[grep("exp", names(melanoma_sp_idents))], function(x) {
    seurat_obj <- spatial_melanoma_psf$spatial_psf_obj
    Idents(seurat_obj) <- melanoma_sp_idents[[x]]
    
    plot_title <- gsub("_", " ", x)
    plot_title <- gsub("vesalius", "Vesalius", plot_title)
    plot_title <- gsub("seurat", "Seurat", plot_title)
    
    SpatialDimPlot(seurat_obj, label = T, label.size = 5) + ggtitle(plot_title)
  }), nrow = 3
), filename = "mel_overlap_vis_plots/exp_spatial_plots.pdf", width = 135, height = 382, units = "mm")

ggsave(patchwork::wrap_plots(
  lapply(names(melanoma_sp_idents)[grep("psf", names(melanoma_sp_idents))], function(x) {
    seurat_obj <- spatial_melanoma_psf$spatial_psf_obj
    Idents(seurat_obj) <- melanoma_sp_idents[[x]]
    
    plot_title <- gsub("_", " ", x)
    plot_title <- gsub("vesalius", "Vesalius", plot_title)
    plot_title <- gsub("seurat", "Seurat", plot_title)
    
    SpatialDimPlot(seurat_obj, label = T, label.size = 5) + ggtitle(plot_title)
  }), nrow = 3
), filename = "mel_overlap_vis_plots/psf_spatial_plots.pdf", width = 135, height = 382, units = "mm")



ggsave(patchwork::wrap_plots(
  lapply(names(mouse_brain_sp_idents)[grep("exp", names(mouse_brain_sp_idents))], function(x) {
    seurat_obj <- spatial_mouse_brain_psf$spatial_psf_obj
    Idents(seurat_obj) <- mouse_brain_sp_idents[[x]]
    
    plot_title <- gsub("_", " ", x)
    plot_title <- gsub("vesalius", "Vesalius", plot_title)
    plot_title <- gsub("seurat", "Seurat", plot_title)
    
    SpatialDimPlot(seurat_obj, label = T, label.size = 5) + ggtitle(plot_title)
  }), nrow = 3
), filename = "mouse_brain_overlap_vis_plots/exp_spatial_plots.pdf", width = 135, height = 382, units = "mm")

ggsave(patchwork::wrap_plots(
  lapply(names(mouse_brain_sp_idents)[grep("psf", names(mouse_brain_sp_idents))], function(x) {
    seurat_obj <- spatial_mouse_brain_psf$spatial_psf_obj
    Idents(seurat_obj) <- mouse_brain_sp_idents[[x]]
    
    plot_title <- gsub("_", " ", x)
    plot_title <- gsub("vesalius", "Vesalius", plot_title)
    plot_title <- gsub("seurat", "Seurat", plot_title)
    
    SpatialDimPlot(seurat_obj, label = T, label.size = 5) + ggtitle(plot_title)
  }), nrow = 3
), filename = "mouse_brain_overlap_vis_plots/psf_spatial_plots.pdf", width = 135, height = 382, units = "mm")



#### Expression and psf based clustering comparison Sankey diagrams ####
## ggplot colors function to color sankey diagram nodes
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

## sankey plot building function
sankey_builder <- function(links, node_colors = NULL) {
  
  links <- links[which(links$value > 0),]
  
  nodes <- data.frame(
    name=c(as.character(links$source), 
           as.character(links$target)) %>% unique()
  )
  
  if(is.null(node_colors)) {
    color_scale_js <- NULL
  } else {
    color_scale_js <- paste0(
      'd3.scaleOrdinal()',
      '.domain([', paste0('"', nodes$name, '"', collapse = ", "), '])',
      '.range([', paste0('"', node_colors[nodes$name], '"', collapse = ", "), '])'
    )
  }
  
  
  
  links$IDsource <- match(links$source, nodes$name)-1 
  links$IDtarget <- match(links$target, nodes$name)-1
  
  # Make the Network
  if(is.null(color_scale_js)) {
    sn <- sankeyNetwork(Links = links, Nodes = nodes,
                        Source = "IDsource", Target = "IDtarget",
                        Value = "value", NodeID = "name",
                        sinksRight=FALSE, fontSize = 20)
  } else {
    sn <- sankeyNetwork(Links = links, Nodes = nodes,
                        Source = "IDsource", Target = "IDtarget",
                        Value = "value", NodeID = "name", colourScale = color_scale_js,
                        sinksRight=FALSE, fontSize = 20)
  }
  
  
  sn$x$nodes <-
    sn$x$nodes %>% 
    mutate(is_source_node = name %in% links$source)
  
  htmlwidgets::onRender(
    sn,
    '
  function(el,x) {
  d3.select(el)
    .selectAll(".node text")
    .filter(function(d) { return d.is_source_node; })
    .attr("x", x.options.nodeWidth - 16)
    .attr("text-anchor", "end");
  
  d3.select(el)
    .selectAll(".node text")
    .filter(function(d) { return !d.is_source_node; })
    .attr("x", x.options.nodeWidth)
    .attr("text-anchor", "start");
  }
  '
  )
  
}


overalap_df_builder <- function(ident_list) {
  
  ident_df = data.frame(barcodes = names(ident_list$seurat_exp), row.names = names(ident_list$seurat_exp),
                        seurat_exp = as.character(ident_list$seurat_exp[names(ident_list$seurat_exp)]),
                        vesalius_exp = as.character(ident_list$vesalius_exp[names(ident_list$seurat_exp)]), 
                        spatialGE_exp = as.character(ident_list$spatialGE_exp[names(ident_list$seurat_exp)]),
                        seurat_psf = as.character(ident_list$seurat_psf[names(ident_list$seurat_exp)]), 
                        vesalius_psf = as.character(ident_list$vesalius_psf[names(ident_list$seurat_exp)]), 
                        spatialGE_psf = as.character(ident_list$spatialGE_psf[names(ident_list$seurat_exp)])
  )
  
  vesalius_to_seurat_exp <- Reduce(rbind, 
                                   lapply(unique(ident_df$vesalius_exp)[order(as.integer(unique(ident_df$vesalius_exp)))], function(x) {
                                     data.frame(source = paste0("Vesalius_", x),
                                                target = paste0("Seurat_", names(table(ident_df$seurat_exp[which(ident_df$vesalius_exp == x)]))),
                                                value = as.integer(unname(table(ident_df$seurat_exp[which(ident_df$vesalius_exp == x)])))
                                     )
                                   })
  )
  
  
  seurat_to_spatialGE_exp <- Reduce(rbind, 
                                    lapply(unique(ident_df$seurat_exp)[order(as.integer(unique(ident_df$seurat_exp)))], function(x) {
                                      data.frame(source = paste0("Seurat_", x),
                                                 target = paste0("spatialGE_", names(table(ident_df$spatialGE_exp[which(ident_df$seurat_exp == x)]))),
                                                 value = as.integer(unname(table(ident_df$spatialGE_exp[which(ident_df$seurat_exp == x)])))
                                      )
                                    })
  )
  
  
  vesalius_to_seurat_psf <- Reduce(rbind, 
                                   lapply(unique(ident_df$vesalius_psf)[order(as.integer(unique(ident_df$vesalius_psf)))], function(x) {
                                     data.frame(source = paste0("Vesalius_", x),
                                                target = paste0("Seurat_", names(table(ident_df$seurat_psf[which(ident_df$vesalius_psf == x)]))),
                                                value = as.integer(unname(table(ident_df$seurat_psf[which(ident_df$vesalius_psf == x)])))
                                     )
                                   })
  )
  
  
  seurat_to_spatialGE_psf <- Reduce(rbind, 
                                    lapply(unique(ident_df$seurat_psf)[order(as.integer(unique(ident_df$seurat_psf)))], function(x) {
                                      data.frame(source = paste0("Seurat_", x),
                                                 target = paste0("spatialGE_", names(table(ident_df$spatialGE_psf[which(ident_df$seurat_psf == x)]))),
                                                 value = as.integer(unname(table(ident_df$spatialGE_psf[which(ident_df$seurat_psf == x)])))
                                      )
                                    })
  )
  
  seurat_exp_to_psf <- Reduce(rbind,
                              lapply(unique(ident_df$seurat_exp)[order(as.integer(unique(ident_df$seurat_exp)))], function(x) {
                                data.frame(source = paste0("Seurat_exp_", x),
                                           target = paste0("Seurat_psf_", names(table(ident_df$seurat_psf[which(ident_df$seurat_exp == x)]))),
                                           value = as.integer(unname(table(ident_df$seurat_psf[which(ident_df$seurat_exp == x)])))
                                )
                              })
  )
  
  vesalius_exp_to_psf <- Reduce(rbind,
                                lapply(unique(ident_df$vesalius_exp)[order(as.integer(unique(ident_df$vesalius_exp)))], function(x) {
                                  data.frame(source = paste0("Vesalius_exp_", x),
                                             target = paste0("Vesalius_psf_", names(table(ident_df$vesalius_psf[which(ident_df$vesalius_exp == x)]))),
                                             value = as.integer(unname(table(ident_df$vesalius_psf[which(ident_df$vesalius_exp == x)])))
                                  )
                                })
  )
  
  
  spatialGE_exp_to_psf <- Reduce(rbind,
                                 lapply(unique(ident_df$spatialGE_exp)[order(as.integer(unique(ident_df$spatialGE_exp)))], function(x) {
                                   data.frame(source = paste0("spatialGE_exp_", x),
                                              target = paste0("spatialGE_psf_", names(table(ident_df$spatialGE_psf[which(ident_df$spatialGE_exp == x)]))),
                                              value = as.integer(unname(table(ident_df$spatialGE_psf[which(ident_df$spatialGE_exp == x)])))
                                   )
                                 })
  )
  
  
  clust_name_convert <- setNames(object = c( "Seurat_",  "Seurat_", "spatialGE_", "spatialGE_", "Vesalius_", "Vesalius_"), names(ident_list))
  
  exp_clust_colors <- Reduce(c, lapply(names(ident_list)[grep("exp", names(ident_list))], function(x) {
    setNames(object = gg_color_hue(length(levels(ident_list[[x]]))), 
             nm = paste0(clust_name_convert[x], levels(ident_list[[x]])))
  }))
  
  psf_clust_colors <- Reduce(c, lapply(names(ident_list)[grep("psf", names(ident_list))], function(x) {
    setNames(object = gg_color_hue(length(levels(ident_list[[x]]))), 
             nm = paste0(clust_name_convert[x], levels(ident_list[[x]])))
  }))
  
  
  clust_name_convert_exp_psf <- setNames(object = c( "Seurat_exp_",  "Seurat_psf_", "spatialGE_exp_", "spatialGE_psf_", "Vesalius_exp_", "Vesalius_psf_"), names(ident_list))
  
  psf_and_exp_clust_colors <- Reduce(c, lapply(names(ident_list), function(x) {
    setNames(object = gg_color_hue(length(levels(ident_list[[x]]))), 
             nm = paste0(clust_name_convert_exp_psf[x], levels(ident_list[[x]])))
  }))
  
  
  return(
    list(exp_clust_links = rbind(vesalius_to_seurat_exp, seurat_to_spatialGE_exp), exp_clust_colors = exp_clust_colors,
         psf_clust_links = rbind(vesalius_to_seurat_psf, seurat_to_spatialGE_psf), psf_clust_colors = psf_clust_colors,
         seurat_exp_to_psf = seurat_exp_to_psf, vesalius_exp_to_psf = vesalius_exp_to_psf, 
         spatialGE_exp_to_psf = spatialGE_exp_to_psf, psf_and_exp_clust_colors = psf_and_exp_clust_colors, ident_df = ident_df
    )
  )
}



library(networkD3)

mel_clust_overlaps <- overalap_df_builder(ident_list = melanoma_sp_idents)

sankey_builder(links = mel_clust_overlaps$exp_clust_links, node_colors = mel_clust_overlaps$exp_clust_colors)
sankey_builder(links = mel_clust_overlaps$psf_clust_links, node_colors = mel_clust_overlaps$psf_clust_colors)

# mel_exp_to_psf_shortened_labels <- lapply(mel_clust_overlaps[grep("exp_to_psf", names(mel_clust_overlaps))], function(x) {
#   node_colors <- mel_clust_overlaps$psf_and_exp_clust_colors[unique(c(x$source, x$target))]
#   names(node_colors) <- gsub("spatialGE_|Vesalius_|Seurat_", "", names(node_colors))
#   
#   x$source <- gsub("spatialGE_|Vesalius_|Seurat_", "", x$source)
#   x$target <- gsub("spatialGE_|Vesalius_|Seurat_", "", x$target)
#   c(paste0(x$source, " [", x$value, "] ", x$target),
#     paste0(":", names(node_colors), " ", unname(node_colors))
#   )
# })

sankey_builder(links = mel_clust_overlaps$seurat_exp_to_psf, node_colors = mel_clust_overlaps$psf_and_exp_clust_colors)
sankey_builder(links = mel_clust_overlaps$vesalius_exp_to_psf, node_colors = mel_clust_overlaps$psf_and_exp_clust_colors)
sankey_builder(links = mel_clust_overlaps$spatialGE_exp_to_psf, node_colors = mel_clust_overlaps$psf_and_exp_clust_colors)


mouse_brain_clust_overlaps <- overalap_df_builder(ident_list = mouse_brain_sp_idents)

sankey_builder(links = mouse_brain_clust_overlaps$exp_clust_links, node_colors = mouse_brain_clust_overlaps$exp_clust_colors)
sankey_builder(links = mouse_brain_clust_overlaps$psf_clust_links, node_colors = mouse_brain_clust_overlaps$psf_clust_colors)

# mosue_brain_exp_to_psf_shortened_labels <- lapply(mouse_brain_clust_overlaps[grep("exp_to_psf", names(mouse_brain_clust_overlaps))], function(x) {
#   node_colors <- mouse_brain_clust_overlaps$psf_and_exp_clust_colors[unique(c(x$source, x$target))]
#   names(node_colors) <- gsub("spatialGE_|Vesalius_|Seurat_", "", names(node_colors))
#   
#   x$source <- gsub("spatialGE_|Vesalius_|Seurat_", "", x$source)
#   x$target <- gsub("spatialGE_|Vesalius_|Seurat_", "", x$target)
#   c(paste0(x$source, " [", x$value, "] ", x$target),
#     paste0(":", names(node_colors), " ", unname(node_colors))
#   )
# })

sankey_builder(links = mouse_brain_clust_overlaps$seurat_exp_to_psf, node_colors = mouse_brain_clust_overlaps$psf_and_exp_clust_colors)
sankey_builder(links = mouse_brain_clust_overlaps$vesalius_exp_to_psf, node_colors = mouse_brain_clust_overlaps$psf_and_exp_clust_colors)
sankey_builder(links = mouse_brain_clust_overlaps$spatialGE_exp_to_psf, node_colors = mouse_brain_clust_overlaps$psf_and_exp_clust_colors)


#### ARI index calculation ####
library(mclust)

## Human melanoma
mclust::adjustedRandIndex(mel_clust_overlaps$ident_df$seurat_exp, mel_clust_overlaps$ident_df$vesalius_exp)
mclust::adjustedRandIndex(mel_clust_overlaps$ident_df$seurat_exp, mel_clust_overlaps$ident_df$spatialGE_exp)
mclust::adjustedRandIndex(mel_clust_overlaps$ident_df$vesalius_exp, mel_clust_overlaps$ident_df$spatialGE_exp)

mclust::adjustedRandIndex(mel_clust_overlaps$ident_df$seurat_psf, mel_clust_overlaps$ident_df$vesalius_psf)
mclust::adjustedRandIndex(mel_clust_overlaps$ident_df$seurat_psf, mel_clust_overlaps$ident_df$spatialGE_psf)
mclust::adjustedRandIndex(mel_clust_overlaps$ident_df$vesalius_psf, mel_clust_overlaps$ident_df$spatialGE_psf)


mclust::adjustedRandIndex(mel_clust_overlaps$ident_df$seurat_exp, mel_clust_overlaps$ident_df$seurat_psf)
mclust::adjustedRandIndex(mel_clust_overlaps$ident_df$spatialGE_exp, mel_clust_overlaps$ident_df$spatialGE_psf)
mclust::adjustedRandIndex(mel_clust_overlaps$ident_df$vesalius_exp, mel_clust_overlaps$ident_df$vesalius_psf)




## Mouse brain
mclust::adjustedRandIndex(mouse_brain_clust_overlaps$ident_df$seurat_exp, mouse_brain_clust_overlaps$ident_df$vesalius_exp)
mclust::adjustedRandIndex(mouse_brain_clust_overlaps$ident_df$seurat_exp, mouse_brain_clust_overlaps$ident_df$spatialGE_exp)
mclust::adjustedRandIndex(mouse_brain_clust_overlaps$ident_df$vesalius_exp, mouse_brain_clust_overlaps$ident_df$spatialGE_exp)

mclust::adjustedRandIndex(mouse_brain_clust_overlaps$ident_df$seurat_psf, mouse_brain_clust_overlaps$ident_df$vesalius_psf)
mclust::adjustedRandIndex(mouse_brain_clust_overlaps$ident_df$seurat_psf, mouse_brain_clust_overlaps$ident_df$spatialGE_psf)
mclust::adjustedRandIndex(mouse_brain_clust_overlaps$ident_df$vesalius_psf, mouse_brain_clust_overlaps$ident_df$spatialGE_psf)


mclust::adjustedRandIndex(mouse_brain_clust_overlaps$ident_df$seurat_exp, mouse_brain_clust_overlaps$ident_df$seurat_psf)
mclust::adjustedRandIndex(mouse_brain_clust_overlaps$ident_df$spatialGE_exp, mouse_brain_clust_overlaps$ident_df$spatialGE_psf)
mclust::adjustedRandIndex(mouse_brain_clust_overlaps$ident_df$vesalius_exp, mouse_brain_clust_overlaps$ident_df$vesalius_psf)


#### Total pathway disturbance for overlapping and PSF specific pathways ####
pdf(file = paste0("psf_branch_ratios_per_method.pdf"), height = 8.27, width = 13, onefile = T)

patchwork::wrap_plots(lapply(names(mel_sp_cluster_specific_pathways_high)[c(1,3,5,2,4,6)], function(x) {
  psf_total_fcs <- Reduce(rbind, lapply(names(mel_sp_cluster_specific_pathways_high[[x]]), function(y) {
    if(!is.null(mel_sp_cluster_specific_pathways_high[[x]][[y]]$altered_pathway_branch_ratios)) {
      cbind(mel_sp_cluster_specific_pathways_high[[x]][[y]]$altered_pathway_branch_ratios, Ident = gsub("spatialGE_|vesalius_|seurat_", "", y))
    }
  }))
  
  psf_total_fcs$class <- ifelse(psf_total_fcs$any_overlap_with_ORA_GSdensity, "Overlapping", "Non overlapping")
  
  psf_total_fcs$PSF_FC_sum <- log(psf_total_fcs$PSF_FC_sum)
  
  plot_title <- gsub("_", " ", x)
  plot_title <- gsub("vesalius", "Vesalius", plot_title)
  plot_title <- gsub("seurat", "Seurat", plot_title)
  
  ggboxplot(psf_total_fcs, x = "Ident", y = "PSF_FC_sum",
            color = "class") +
    labs(title = "", x = "Ident", y = "Total log pathway FC", fill = "class") +
    theme_minimal() +
    # facet_wrap(~method, nrow = 2, drop = T) +
    theme(axis.text.x = element_text(size = 14, angle = 45, vjust=1, hjust=1), axis.title.x = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"),
          legend.text=element_text(size=16), legend.title = element_text(size=16), strip.text = element_text(size=20),
          plot.title = element_text(size = 22, face = "bold")) + ggtitle(plot_title)
  
  # boxplot(log(psf_total_fcs$PSF_FC_sum) ~ psf_total_fcs$class, 
  #         main = x, ylab = "Total pathway FC", xlab = "Group")
}), guides = 'collect') + plot_layout(axes = 'collect')



patchwork::wrap_plots(lapply(names(mouse_brain_sp_cluster_specific_pathways_high)[c(1,3,5,2,4,6)], function(x) {
  psf_total_fcs <- Reduce(rbind, lapply(names(mouse_brain_sp_cluster_specific_pathways_high[[x]]), function(y) {
    if(!is.null(mouse_brain_sp_cluster_specific_pathways_high[[x]][[y]]$altered_pathway_branch_ratios)) {
      cbind(mouse_brain_sp_cluster_specific_pathways_high[[x]][[y]]$altered_pathway_branch_ratios, Ident = gsub("spatialGE_|vesalius_|seurat_", "", y))
    }
  }))
  
  psf_total_fcs$class <- ifelse(psf_total_fcs$any_overlap_with_ORA_GSdensity, "Overlapping", "Non overlapping")
  
  psf_total_fcs$PSF_FC_sum <- log(psf_total_fcs$PSF_FC_sum)
  
  plot_title <- gsub("_", " ", x)
  plot_title <- gsub("vesalius", "Vesalius", plot_title)
  plot_title <- gsub("seurat", "Seurat", plot_title)
  
  ggboxplot(psf_total_fcs, x = "Ident", y = "PSF_FC_sum",
            color = "class") +
    labs(title = "", x = "Ident", y = "Total log pathway FC", fill = "class") +
    theme_minimal() +
    # facet_wrap(~method, nrow = 2, drop = T) +
    theme(axis.text.x = element_text(size = 14, angle = 45, vjust=1, hjust=1), axis.title.x = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"),
          legend.text=element_text(size=16), legend.title = element_text(size=16), strip.text = element_text(size=20),
          plot.title = element_text(size = 22, face = "bold")) + ggtitle(plot_title)
  
  # boxplot(log(psf_total_fcs$PSF_FC_sum) ~ psf_total_fcs$class, 
  #         main = x, ylab = "Total pathway FC", xlab = "Group")
}), guides = 'collect') + plot_layout(axes = 'collect')

dev.off()



#### Comparison of spatially aware pathway analysis methods ####

## plotting overlaps between three methods in venn diagramm
pdf(file = "spatial_pathway_overlaps.pdf", paper='A4r', height = 8.27, width = 11.69, onefile = T)
layout(matrix(1:2, nrow = 1))
mel_spatially_relevant_branches <- rownames(mel_GSdensity_results$psf_spatial.klds[mel_GSdensity_results$psf_spatial.klds$spatial.kld > quantile(mel_GSdensity_results$psf_spatial.klds$spatial.kld, 0.8, na.rm = T), ])
mel_spatially_relevant_branches <- mel_spatially_relevant_branches[which(mel_spatially_relevant_branches != "NA")]
mel_spatially_relevant_pathways <- unique(sapply(mel_spatially_relevant_branches, function(z) {gsub("-", "_", unlist(strsplit(z, split = "; "))[2])}))
venn::venn(list(PSF = mel_spatially_relevant_pathways,
                GSDensity = rownames(mel_GSdensity_results$exp_spatial.klds[mel_GSdensity_results$exp_spatial.klds$spatial.kld > quantile(mel_GSdensity_results$exp_spatial.klds$spatial.kld, 0.8, na.rm = T), ]),
                spatialGE = mel_spatilaGE$stenrich_df$sample1$gene_set[which(mel_spatilaGE$stenrich_df$sample1$adj_p_value < 0.05)]
),zcolor = "style", box = F, ilabels = "counts", lty = 1, ilcs = 1.3, sncs = 1.5, lwd = 1.3)
title(main = "Human melanoma", line = -2, cex.main = 1.5)

mb_spatially_relevant_branches <- rownames(mouse_brain_GSdensity_results$psf_spatial.klds[mouse_brain_GSdensity_results$psf_spatial.klds$spatial.kld > quantile(mouse_brain_GSdensity_results$psf_spatial.klds$spatial.kld, 0.8, na.rm = T), ])
mb_spatially_relevant_branches <- mb_spatially_relevant_branches[which(mb_spatially_relevant_branches != "NA")]
mb_spatially_relevant_pathways <- unique(sapply(mb_spatially_relevant_branches, function(z) {gsub("-", "_", unlist(strsplit(z, split = "; "))[2])}))
venn::venn(list(PSF = mb_spatially_relevant_pathways,
                GSDensity = rownames(mouse_brain_GSdensity_results$exp_spatial.klds[mouse_brain_GSdensity_results$exp_spatial.klds$spatial.kld > quantile(mouse_brain_GSdensity_results$exp_spatial.klds$spatial.kld, 0.8, na.rm = T), ]),
                spatialGE = mouse_brain_spatialGE$stenrich_df$sample1$gene_set[which(mouse_brain_spatialGE$stenrich_df$sample1$adj_p_value < 0.05)]
), zcolor = "style", box = F, ilabels = "counts", lty = 1, ilcs = 1.3, sncs = 1.5, lwd = 1.3)
title(main = "Mouse brain", line = -2, cex.main = 1.5)
dev.off()


## checking affected sink rations of overlapping and PSF specific pathways in melanoma dataset
mean(
  sapply(intersect(mel_spatially_relevant_pathways, c(rownames(mel_GSdensity_results$exp_spatial.klds[mel_GSdensity_results$exp_spatial.klds$spatial.kld > quantile(mel_GSdensity_results$exp_spatial.klds$spatial.kld, 0.8, na.rm = T), ]),
                                                      spatialGE_env$mel_stenrich_df$sample1$gene_set[which(spatialGE_env$mel_stenrich_df$sample1$adj_p_value < 0.05)])
  ), function(x) {
    sum(grepl(x, mel_spatially_relevant_branches))/pathwya_sink_nums[x]
  })
)

mean(
  sapply(setdiff(mel_spatially_relevant_pathways, c(rownames(mel_GSdensity_results$exp_spatial.klds[mel_GSdensity_results$exp_spatial.klds$spatial.kld > quantile(mel_GSdensity_results$exp_spatial.klds$spatial.kld, 0.8, na.rm = T), ]),
                                                    spatialGE_env$mel_stenrich_df$sample1$gene_set[which(spatialGE_env$mel_stenrich_df$sample1$adj_p_value < 0.05)])
  ), function(x) {
    sum(grepl(x, mel_spatially_relevant_branches))/pathwya_sink_nums[x]
  })
)



## Plotting some example of spatially relevant pathways of human melanoma for GSDensity and PSF
mel_overlapping_spatial_pathways <- Reduce(intersect, 
                                           list(PSF = mel_spatially_relevant_pathways,
                                                GSDensity = rownames(mel_GSdensity_results$exp_spatial.klds[mel_GSdensity_results$exp_spatial.klds$spatial.kld > quantile(mel_GSdensity_results$exp_spatial.klds$spatial.kld, 0.8, na.rm = T), ]),
                                                spatialGE = spatialGE_env$mel_stenrich_df$sample1$gene_set[which(spatialGE_env$mel_stenrich_df$sample1$adj_p_value < 0.05)])
)


mel_sp_psf_mat <- as.matrix(psf_spatial$spatial_psf_obj@assays$Spatial@counts)


mel_selected_spatial_sinks <- setNames(object = c("14; Cytosolic-DNA-sensing-pathway", "53; NF-kappa-B-signaling-pathway", 
                                                  "42; PPAR-signaling-pathway", "36; TNF-signaling-pathway", "19; Toll-like-receptor-signaling-pathway"), 
                                       nm = mel_overlapping_spatial_pathways)


pdf(file = "mel_psf_gsdensity_featureplots.pdf", paper='A4r', height = 8.27, width = 9, onefile = T)
patchwork::wrap_plots(lapply(mel_overlapping_spatial_pathways[c(1,5)], function(x) {
  seurat_obj <- psf_spatial$spatial_obj
  seurat_obj@meta.data[,"GSDensity"] <- mel_GSdensity_results$cv.df[rownames(seurat_obj@meta.data), x]
  
  fc <- mel_spatial_psf_mat[mel_selected_spatial_sinks[x],rownames(seurat_obj@meta.data)]
  log_fc <- log2(fc)
  log_fc[which(log_fc < 0)] <- 0
  seurat_obj@meta.data[,"PSF"] <- log_fc
  
  GSDensity_plot <- SpatialFeaturePlot(seurat_obj, features = "GSDensity") +
    theme(legend.position = "right", legend.text=element_text(size=10), legend.title = element_text(size=11), plot.title = element_text(size = 22, face = "bold")) + ggtitle(gsub("_", " ", x))
  psf_plot <- SpatialFeaturePlot(seurat_obj, features = "PSF") +
    theme(legend.position = "right", legend.text=element_text(size=10), legend.title = element_text(size=11))
  
  GSDensity_plot + psf_plot
}), nrow = 2)
dev.off()


## Plotting PSF specific spatially activated pathway branch examples
pdf(file = "mel_psf_specific_featureplots.pdf", paper='A4r', height = 8.27, width = 9, onefile = T)
patchwork::wrap_plots(lapply(c("80; Hippo-signaling-pathway", "56; ErbB-signaling-pathway"), function(x) {
  seurat_obj <- psf_spatial$spatial_obj
  
  fc <- mel_spatial_psf_mat[x,rownames(seurat_obj@meta.data)]
  log_fc <- log2(fc)
  # log_fc[which(log_fc < 0)] <- 0
  seurat_obj@meta.data[, "Log2 PSF"] <- log_fc
  
  SpatialFeaturePlot(seurat_obj, features = "Log2 PSF") +
    theme(legend.position = "right", legend.text=element_text(size=10), legend.title = element_text(size=11), plot.title = element_text(size = 18)) + ggtitle(unlist(strsplit(x, split = "; "))[2])
}), nrow = 1)
dev.off()




## Plotting some example of spatially relevant pathways of mouse brain for GSDensity and PSF
mouse_brain_overlapping_spatial_pathways <- intersect(mb_spatially_relevant_pathways, 
                                                      rownames(mouse_brain_GSdensity_results$exp_spatial.klds[mouse_brain_GSdensity_results$exp_spatial.klds$spatial.kld > quantile(mouse_brain_GSdensity_results$exp_spatial.klds$spatial.kld, 0.8, na.rm = T), ]))


mouse_brain_selected_spatial_sinks <- setNames(object = c("58; Fc_gamma_R_mediated_phagocytosis", "49; Wnt_signaling_pathway"), 
                                               nm = mouse_brain_overlapping_spatial_pathways)


pdf(file = "mouse_brain_psf_gsdensity_featureplots.pdf", paper='A4r', height = 8.27, width = 9, onefile = T)
patchwork::wrap_plots(lapply(mouse_brain_overlapping_spatial_pathways, function(x) {
  seurat_obj <- spatial_mouse_brain_psf$spatial_obj
  seurat_obj@meta.data[,"GSDensity"] <- mouse_brain_GSdensity_results$cv.df[rownames(seurat_obj@meta.data), x]
  
  fc <- mouse_brain_spatial_psf_mat[mouse_brain_selected_spatial_sinks[x],rownames(seurat_obj@meta.data)]
  log_fc <- log2(fc)
  log_fc[which(log_fc < 0)] <- 0
  seurat_obj@meta.data[,"PSF"] <- log_fc
  
  GSDensity_plot <- SpatialFeaturePlot(seurat_obj, features = "GSDensity") +
    theme(legend.position = "right", legend.text=element_text(size=10), legend.title = element_text(size=11), plot.title = element_text(size = 22, face = "bold")) + ggtitle(gsub("_", " ", x))
  psf_plot <- SpatialFeaturePlot(seurat_obj, features = "PSF") +
    theme(legend.position = "right", legend.text=element_text(size=10), legend.title = element_text(size=11))
  
  GSDensity_plot + psf_plot
}), nrow = 2)
dev.off()

#### Cluster borders ligand receptor interaction detection in human melanoma dataset ####
### Customized linked dimplot shiny app for manual spot selection
modified_linked_dimplot <-function (object, dims = 1:2, reduction = NULL, image = NULL, 
                                    group.by = NULL, alpha = c(0.1, 1), combine = TRUE) 
{
  ui <- miniPage(gadgetTitleBar(title = "LinkedDimPlot", left = miniTitleBarButton(inputId = "reset", label = "Reset")), 
                 miniContentPanel(
                   actionButton("record_selection", label = "Record"),
                   fillRow(
                     plotOutput(outputId = "spatialplot", height = "100%", click = clickOpts(id = "spclick", clip = TRUE),hover = hoverOpts(id = "sphover", delay = 10, nullOutside = TRUE)), 
                     plotOutput(outputId = "dimplot", height = "100%", 
                                brush = brushOpts(id = "brush", delay = 10, clip = TRUE, resetOnNew = FALSE), 
                                click = clickOpts(id = "dimclick",clip = TRUE), 
                                hover = hoverOpts(id = "dimhover", delay = 10, nullOutside = TRUE)), height = "97%"
                   ), 
                   verbatimTextOutput(outputId = "info")
                   
                 )
  )
  image <- image %||% Seurat:::DefaultImage(object = object)
  cells.use <- Cells(x = object[[image]])
  reduction <- reduction %||% DefaultDimReduc(object = object)
  dims <- dims[1:2]
  dims <- paste0(Key(object = object[[reduction]]), dims)
  group.by <- group.by %||% "ident"
  group.data <- FetchData(object = object, vars = group.by, 
                          cells = cells.use)
  coords <- GetTissueCoordinates(object = object[[image]])
  embeddings <- Embeddings(object = object[[reduction]])[cells.use, 
                                                         dims]
  plot.data <- cbind(coords, group.data, embeddings)
  plot.data$selected_ <- FALSE
  Idents(object = object) <- group.by
  server <- function(input, output, session) {
    click <- reactiveValues(pt = NULL, invert = FALSE)
    plot.env <- reactiveValues(data = plot.data, alpha.by = NULL, selected_cells = c(), clust_num = 0)
    observeEvent(eventExpr = input$done, handlerExpr = {
      plots <- list(plot.env$spatialplot, plot.env$dimplot)
      if (combine) {
        plots <- wrap_plots(plots, ncol = 2)
      }
      stopApp(returnValue = plots)
    })
    observeEvent(eventExpr = input$reset, handlerExpr = {
      click$pt <- NULL
      click$invert <- FALSE
      session$resetBrush(brushId = "brush")
    })
    observeEvent(eventExpr = input$brush, handlerExpr = click$pt <- NULL)
    observeEvent(eventExpr = input$spclick, handlerExpr = {
      click$pt <- input$spclick
      click$invert <- TRUE
      
      clicked <- nearPoints(df = plot.data, coordinfo = if (click$invert) {
        Seurat:::InvertCoordinate(x = click$pt)
      }
      else {
        click$pt
      }, threshold = 10, maxpoints = 1)
      
      
      plot.env$selected_cells <- c(plot.env$selected_cells, rownames(clicked))
      
    })
    
    observeEvent(input$record_selection, {
      
      write(plot.env$selected_cells, file = paste0("cluster", plot.env$clust_num))
      
      plot.env$clust_num = plot.env$clust_num + 1
      plot.env$selected_cells <- c()
      
    })
    
    observeEvent(eventExpr = input$dimclick, handlerExpr = {
      click$pt <- input$dimclick
      click$invert <- FALSE
    })
    observeEvent(eventExpr = c(input$brush, # input$spclick, 
                               input$dimclick), handlerExpr = {
                                 plot.env$data <- if (is.null(x = input$brush)) {
                                   clicked <- nearPoints(df = plot.data, coordinfo = if (click$invert) {
                                     Seurat:::InvertCoordinate(x = click$pt)
                                   }
                                   else {
                                     click$pt
                                   }, threshold = 10, maxpoints = 1)
                                   if (nrow(x = clicked) == 1) {
                                     cell.clicked <- rownames(x = clicked)
                                     group.clicked <- plot.data[cell.clicked, group.by, 
                                                                drop = TRUE]
                                     idx.group <- which(x = plot.data[[group.by]] == 
                                                          group.clicked)
                                     plot.data[idx.group, "selected_"] <- TRUE
                                     plot.data
                                   }
                                   else {
                                     plot.data
                                   }
                                 }
                                 else if (input$brush$outputId == "dimplot") {
                                   brushedPoints(df = plot.data, brush = input$brush, 
                                                 allRows = TRUE)
                                 }
                                 else if (input$brush$outputId == "spatialplot") {
                                   brushedPoints(df = plot.data, brush = Seurat:::InvertCoordinate(x = input$brush), 
                                                 allRows = TRUE)
                                 }
                                 plot.env$alpha.by <- if (any(plot.env$data$selected_)) {
                                   "selected_"
                                 }
                                 else {
                                   NULL
                                 }
                               })
    output$spatialplot <- renderPlot(expr = {
      plot.env$spatialplot <- SingleSpatialPlot(data = plot.env$data, 
                                                image = object[[image]], col.by = group.by, pt.size.factor = 1.6, 
                                                crop = TRUE, alpha.by = plot.env$alpha.by) + 
        scale_alpha_ordinal(range = alpha) + NoLegend()
      plot.env$spatialplot
    })
    output$dimplot <- renderPlot(expr = {
      plot.env$dimplot <- SingleDimPlot(data = plot.env$data, 
                                        dims = dims, col.by = group.by, alpha.by = plot.env$alpha.by) + 
        scale_alpha_ordinal(range = alpha) + guides(alpha = "none")
      plot.env$dimplot
    })
    output$info <- renderPrint(expr = {
      cell.hover <- rownames(x = nearPoints(df = plot.data, 
                                            coordinfo = if (is.null(x = input[["sphover"]])) {
                                              input$dimhover
                                            }
                                            else {
                                              Seurat:::InvertCoordinate(x = input$sphover)
                                            }, threshold = 10, maxpoints = 1))
      if (length(x = cell.hover) == 1) {
        paste(cell.hover, paste("Group:", plot.data[cell.hover, 
                                                    group.by, drop = TRUE]), collapse = "<br />")
      }
      else {
        NULL
      }
    })
    
  }
  runGadget(app = ui, server = server)
}


### load previously manually selected spot groups
load("psf_manual_clusters.RData")

manual_cluster_labels <- setNames(object = rep("non_selected", length(rownames(spatial_melanoma_psf$spatial_psf_obj@meta.data))), nm = rownames(spatial_melanoma_psf$spatial_psf_obj@meta.data))
manual_cluster_labels[names(psf_manual_clusters)] <- unname(psf_manual_clusters)


spatial_melanoma_psf$spatial_psf_obj@meta.data$manual_selection <- manual_cluster_labels[rownames(spatial_melanoma_psf$spatial_psf_obj@meta.data)]


cluster_selection_list <- lapply(unique(manual_cluster_labels), function(x) {
  names(manual_cluster_labels)[which(manual_cluster_labels == x)]
})

names(cluster_selection_list) <- unique(manual_cluster_labels)
names(cluster_selection_list) <- gsub("cluster", "c", names(cluster_selection_list))


### selecteing spot groups of interest
manual_selection_subset <- cluster_selection_list[c("c4", "c2", "c18", "c12", "c10", "c13", "c16", "c20")]
manual_selection_subset_vec <- unlist(lapply(names(manual_selection_subset), function(x) {
  setNames(object = rep(x, length(manual_selection_subset[[x]])), nm = manual_selection_subset[[x]])
}))


### visualizing selected spot groups
spatial_melanoma_psf$spatial_psf_obj@meta.data$custom_selection <- as.character(spatial_melanoma_psf$spatial_psf_obj@meta.data$seurat_clusters)
spatial_melanoma_psf$spatial_psf_obj@meta.data[names(manual_selection_subset_vec),"custom_selection"] <- unname(manual_selection_subset_vec)

cluster_selection_list_subset <- sapply(unique(spatial_melanoma_psf$spatial_psf_obj@meta.data$custom_selection), function(x) {
  rownames(spatial_melanoma_psf$spatial_psf_obj@meta.data)[which(spatial_melanoma_psf$spatial_psf_obj@meta.data$custom_selection == x)]
})


### Spatial plot of 
SpatialDimPlot(spatial_melanoma_psf$spatial_psf_obj, label = F, label.size = 6,
               cols.highlight =  c(gg_color_hue(9), c("#fcc244", "#fcc244", "#e33a14", "#e33a14", "#006400", "#006400", "#fcc244", "#006400")),
               cells.highlight = cluster_selection_list_subset[order(names(cluster_selection_list_subset))])



### performingDetecting significantly differentially activated pathways branches between border and core spot groups for each cluster
psf_clust_comparisons <- list(c_12_10 = c("c12", "c10"),
                              c_18_4 = c("c18", "c4"),
                              c_20_10 = c("c20", "c10"),
                              c_13_16 = c("c13", "c16")
)

spatial_melanoma_psf$spatial_psf_obj@meta.data$manual_selection <- gsub("cluster", "c", spatial_melanoma_psf$spatial_psf_obj@meta.data$manual_selection)

psf_adjesent_clust_comparison <- lapply(psf_clust_comparisons, function(y) {
  
  clust_0 <- mel_spatial_psf_mat[,rownames(spatial_melanoma_psf$spatial_psf_obj@meta.data)[which(spatial_melanoma_psf$spatial_psf_obj@meta.data$manual_selection == y[1])]]
  clust_1 <- mel_spatial_psf_mat[,rownames(spatial_melanoma_psf$spatial_psf_obj@meta.data)[which(spatial_melanoma_psf$spatial_psf_obj@meta.data$manual_selection == y[2])]]
  
  clust_0_and_1_psf <- cbind(clust_0, clust_1)
  
  gr <- factor(c(rep("clus0", times=ncol(clust_0)), rep("clus1", times=ncol(clust_1))), levels=c("clus0", "clus1"))
  
  p.val <- apply(clust_0_and_1_psf, 1, function(x) {
    fit <- lm(x~gr)
    summary(fit)$coefficients[2,4]
  })
  
  stat_df <- data.frame(sink = rownames(clust_0_and_1_psf), log2_fc = log2(rowMeans(clust_0)/rowMeans(clust_1)),  p_val = p.val, fdr = p.adjust(p.val, method = "fdr"))
  
  stat_df <- stat_df[order(stat_df$fdr),]
  
  stat_df[which(stat_df$fdr < 0.05),]
  
})

### annotating significant branches
psf_adjesent_clust_comparison <- lapply(psf_adjesent_clust_comparison, function(x) {
  if(nrow(x) > 0) {
    cbind(x, sink_id = kegg_sink_to_process[x$sink, "sink_id"], 
          sink_name = kegg_sink_to_process[x$sink, "Sink"],
          pathway = kegg_sink_to_process[x$sink, "Pathway_name"],
          process = kegg_sink_to_process[x$sink, "Process"], 
          downstream_pathway = kegg_sink_to_process[x$sink, "Pathway"])
  } else {
    x
  }
})


### Detecting receptor ligand connections between KEGG signaling pathways
determine.input.nodes <- function(pathway) {
  dfs <- graphnel_to_df(pathway)
  return(setdiff(dfs$edge_table$from, dfs$edge_table$to))
}

kegg_pathway_connection_df <- Reduce(rbind, lapply(names(kegg_curated_40_signalings), function(x) {
  
  sink_genes <- graph::nodeData(kegg_curated_40_signalings[[x]]$graph, 
                                n = kegg_curated_40_signalings[[x]]$sink.nodes, attr = "genes")
  
  connections <- lapply(names(kegg_curated_40_signalings), function(y) {
    
    input_genes <- graph::nodeData(kegg_curated_40_signalings[[y]]$graph, 
                                   n = determine.input.nodes(kegg_curated_40_signalings[[y]]), attr = "genes")
    
    Reduce(rbind, lapply(names(sink_genes), function(z) {
      Reduce(rbind, lapply(names(input_genes), function(h) {
        if(length(intersect(sink_genes[[z]], input_genes[[h]])) > 0) {
          data.frame(sink = paste0(z, "; ", x), input_node = h, input_pathway = y)
          # paste0(z, " > ", h)
        }
      }))
      
    }))   
    
  })
  
  connections[!sapply(connections, is.null)]
  
  Reduce(rbind, connections)
  
}))


### Checking if there are interacting significant pathways in adjacent border spot groups of two clusters in melanoma
adjacent_border_pairs <- list(
  c("c_18_4", "c_12_10"),
  c("c_13_16", "c_20_10")
)


lapply(adjacent_border_pairs, function(x) {
  filtered_connection_df <- kegg_pathway_connection_df[which(kegg_pathway_connection_df$sink %in% rownames(psf_adjesent_clust_comparison[[x[1]]])),]
  filtered_connection_df[which(filtered_connection_df$input_pathway %in% psf_adjesent_clust_comparison[[x[2]]]$pathway),]
})


lapply(adjacent_border_pairs, function(x) {
  filtered_connection_df <- kegg_pathway_connection_df[which(kegg_pathway_connection_df$sink %in% rownames(psf_adjesent_clust_comparison[[x[2]]])),]
  filtered_connection_df[which(filtered_connection_df$input_pathway %in% psf_adjesent_clust_comparison[[x[1]]]$pathway),]
})




#### Calculation of PSF with different number of pseudo bulked cells on GETX single cell data ####
library(zellkonverter)
library(SingleCellExperiment)
library(psf)
library(biomaRt)

load(system.file("extdata", "kegg_curated_40_signalings.RData", package="psf"))

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

gtex_exp_mat <- as.matrix(assay(gtex_single_cell))

rownames(gtex_exp_mat) <- gene_symbol_to_entrez[rownames(gtex_exp_mat)]
# gtex_exp_mat <- gtex_exp_mat[which(!is.na(rownames(gtex_exp_mat))),]  


# chemokine_signaling_genes <- unique(unlist(graph::nodeData(kegg_curated_40_signalings$Chemokine_signaling_pathway$graph, attr = "genes"))) 

signaling_pathway_genes <- unique(unlist(lapply(kegg_curated_40_signalings, function(x) {
  unique(unlist(graph::nodeData(x$graph, attr = "genes")))
})))


gtex_exp_mat <- gtex_exp_mat[intersect(rownames(gtex_exp_mat), signaling_pathway_genes),]

# n_cell_vec <- c(1, 2, 5, 10, 20, 40, 50, 100, 150, 200, 250, 300, 350, 400)
n_cell_vec <- c(1, seq(2, 50, 2))


pseudobulk_mat <- Reduce(cbind, lapply(n_cell_vec, function(z) {
  exp_mat <- Reduce(cbind, lapply(unique(gtex_coldata$tissue), function(x) {
    
    set.seed(1)
    tissue_samples <- sample(rownames(gtex_coldata)[which(gtex_coldata$tissue == x)], 1000)
    print(z)
    
    tissue_matrix <- Reduce(cbind,
                            lapply(0:(1000%/%z - 1), function(y) {
                              
                              cell_range <- 1:z + z*y
                              
                              rowSums(gtex_exp_mat[,tissue_samples[cell_range], drop = F])
                              
                            })
    )
    colnames(tissue_matrix) <- paste0(z, "_", x, 1:ncol(tissue_matrix))
    tissue_matrix
  }))
  (exp_mat + 1)/rowMeans(exp_mat + 1)
}))


pseudobalk_psf <- run_psf(entrez.fc = pseudobulk_mat, kegg.collection = kegg_curated_40_signalings, calculate.significance = F, ncores = 30)



pdf(file = "mean_psf_by_n_cells.pdf", height = 8.27, width = 18, onefile = T)

lapply(pseudobalk_psf, function(x) {
  layout(matrix(1:8, nrow = 2))
  lapply(unique(gtex_coldata$tissue), function(y) {
    tissue_psf_mat <- x$psf_activities[x$sink.nodes,grep(y, colnames(x$psf_activities))]
    print(x$attrs$title)
    boxplot(log(colMeans(tissue_psf_mat)) ~ factor(sapply(colnames(tissue_psf_mat), function(z) {unlist(strsplit(z, "_"))[1]}), levels = n_cell_vec),
            xlab = "N cells", ylab = "Mean log PSF", main = paste0(x$attrs$title, "_", y)
    )
  })
})
dev.off()