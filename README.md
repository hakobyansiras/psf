# PSF toolkit

<!-- badges: start -->
[![R-CMD-check](https://github.com/hakobyansiras/psf/workflows/R-CMD-check/badge.svg)](https://github.com/hakobyansiras/psf/actions)
<!-- badges: end -->

The **PSF Toolkit** is an R package developed for **topology-aware (TA) pathway analysis** of various types of omics data. The package includes interactive modules for **pathway editing** and **visualization**, facilitating pathway curation and results interpretation.

## Installation

```r
if (!require("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github("hakobyansiras/psf")
```

> **Note:** Due to potential network issues with GitHub servers, the package download may take longer than expected. If you encounter timeout errors, try increasing the timeout limit using:
> 
> ```r
> options(timeout = 400)
> ```

## Usage

### Running the Shiny App for Pathway Curation and Visualization

The PSF Toolkit includes a **Shiny** app for interactive pathway curation and visualization.

```r
library(psf)
run_shiny_app()
```

A detailed tutorial on pathway import, curation, and analysis is available [here](#).

## **Running PSF on Bulk RNA-seq Data with 40 Curated KEGG Signaling Pathways**

#### Loading Melanoma Expression Data (GSE112509)

```r
load(system.file("extdata", "melanoma_exp.RData", package="psf"))
```

#### Loading Curated KEGG Pathways

```r
load(system.file("extdata", "kegg_curated_40_signalings.RData", package="psf"))
```

#### Calculating Fold Change (FC) Values

```r
melanoma_fc <- (melanoma_deseq_normalized_counts + 1) / rowMeans(melanoma_deseq_normalized_counts)
```

#### Running PSF Analysis on 40 Pathways

```r
melanoma_psf <- run_psf(entrez.fc = melanoma_fc, kegg.collection = edited_pathways_new, 
                        calculate.significance = FALSE, ncores = 4)
```

#### Plotting the PI3K-Akt Signaling Pathway

```r
plot_pathway(melanoma_psf$PI3K_Akt_signaling_pathway, plot_type = "kegg", 
             color_nodes = "psf_activities", sample_id = "mean", log_norm = TRUE, 
             plot_sink_values = TRUE, use_old_images = TRUE)
```

#### Generating PDF Reports

```r
generate_psf_report(psf_list = melanoma_psf, folder_name = "example_psf_report", 
                    plot_type = "kegg", log_norm = TRUE, use_old_images = TRUE)
```

## **Running PSF on Spatial Transcriptomics Data**

> **Important:** Due to installation issues with the `Seurat` package, it is included in `Suggests` rather than `Imports`. The PSF package will still install even if `Seurat` fails to install. However, you must manually install `Seurat` to perform spatial transcriptomic analyses.

#### Loading a Spatial Dataset

```r
spatial_melanoma <- Load10X_Spatial("/path/to/spatial_dataset", 
                                    filename = "CytAssist_FFPE_Human_Skin_Melanoma_filtered_feature_bc_matrix.h5")
```

#### Downloading Gene Symbol to Entrez ID Conversion Data

```r
load("gene_symbol_to_entrez.RData")
```

Alternatively, you can generate it using **biomaRt**:

```r
# ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
# gene_symbol_to_entrez <- getBM(attributes = c('entrezgene_id', 'hgnc_symbol'), mart = ensembl)
# gene_symbol_to_entrez <- gene_symbol_to_entrez[!is.na(gene_symbol_to_entrez$entrezgene_id), ]
# gene_symbol_to_entrez <- gene_symbol_to_entrez[!duplicated(gene_symbol_to_entrez$entrezgene_id), ]
# gene_symbol_to_entrez <- setNames(gene_symbol_to_entrez$entrezgene_id, gene_symbol_to_entrez$hgnc_symbol)
```

#### Loading KEGG Signaling Pathways

```r
load(system.file("extdata", "kegg_curated_40_signalings.RData", package="psf"))
load(system.file("extdata", "kegg_sink_to_process.RData", package="psf"))
```

#### Running Pathway Analysis

```r
psf_spatial <- spatial_psf_analysis(spatial_obj = spatial_melanoma, 
                                    pathway_collection = kegg_curated_40_signalings, 
                                    gene_symbol_to_entrez = gene_symbol_to_entrez, nthreads = 30)
```

#### Launching the Spatial PSF Browser App

```r
run_psf_spatial_browser(psf_spatial_results = psf_spatial)
```

#### App Usage

A short demo video demonstrating the PSF spatial browser is available: [PSF Spatial Browser Demo](https://www.youtube.com/watch?v=lHTgYBA374o)

You can also explore a **demo Shiny app**: [PSF Spatial Browser](https://apps.armlifebank.am/PSF_spatial_browser/)

The app includes a **dedicated Help tab** containing detailed instructions on how to use it.

## **Running PSF on Telomere Maintenance Mechanism (TMM) Pathways**

#### Loading Melanoma Expression Data

```r
load(system.file("extdata", "melanoma_exp.RData", package="psf"))
```

#### Loading TMM Pathway

```r
load(system.file("extdata", "tmm_pathway.RData", package="psf"))
```

#### Calculating Fold Change (FC) Values

```r
melanoma_fc <- (melanoma_deseq_normalized_counts + 1) / rowMeans(melanoma_deseq_normalized_counts)
```

#### Running PSF Analysis on TMM Pathway

```r
psf_tmm_result_mel <- psf.from.env.entrez.fc(entrez.fc = melanoma_fc,
                                             kegg.collection = tmm_pathway, 
                                             calculate.significance = FALSE, sum = FALSE, split = FALSE, 
                                             return_only_signals = FALSE, tmm_mode = TRUE, tmm_updated_mode = TRUE)
```

#### Plotting the TMM Pathway for a Sample

```r
plot_tmm_pathway(pathway = psf_tmm_result_mel$X031$tmm, no_color_mode = FALSE,
                 mapping_data_type = "signal", log_norm = FALSE, layout = "layout_nicely")
```

For more details, please refer to the package documentation and tutorials.