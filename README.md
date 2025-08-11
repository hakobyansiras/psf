# PSF toolkit

<!-- badges: start -->
[![R-CMD-check](https://github.com/hakobyansiras/psf/workflows/R-CMD-check/badge.svg)](https://github.com/hakobyansiras/psf/actions)
<!-- badges: end -->

The **PSF Toolkit** is an R package developed for **topology-aware (TA) pathway analysis** of various types of omics data. The package includes interactive modules for **pathway editing** and **visualization**, facilitating pathway curation and results interpretation.

- Website: https://hakobyansiras.github.io/PSFC
- Get started: see the Articles tab on the website.

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

A detailed tutorial on pathway import, curation, and analysis is available on the website.

## Build documentation site locally

```r
# Install development dependencies
install.packages(c("roxygen2", "pkgdown"), repos = "https://cloud.r-project.org")

# Generate Rd files (optional, already present)
roxygen2::roxygenise()

# Build site into docs/
pkgdown::build_site()
```

## **Running PSF on Bulk RNA-seq Data with 40 Curated KEGG Signaling Pathways**

```r
library(psf)

# Load melanoma DESeq2 normalized expression data (GSE112509)
load(system.file("extdata", "melanoma_exp.RData", package="psf"))

# Load curated KEGG pathway collection
load(system.file("extdata", "kegg_curated_40_signalings.RData", package="psf"))

# Calculate fold change (FC) values against the global mean
melanoma_fc <- (melanoma_deseq_normalized_counts + 1) / rowMeans(melanoma_deseq_normalized_counts + 1)

# Run PSF analysis on 40 pathways for all samples
melanoma_psf <- run_psf(entrez.fc = melanoma_fc, kegg.collection = kegg_curated_40_signalings, 
                        calculate.significance = FALSE, ncores = 4)

# Plot KEGG image-based PI3K-Akt signaling pathway with color-coded nodes
plot_pathway(melanoma_psf$PI3K_Akt_signaling_pathway, plot_type = "kegg", 
             color_nodes = "psf_activities", sample_id = "mean", log_norm = TRUE, 
             plot_sink_values = TRUE, use_old_images = TRUE)

# Generate PDF reports for PSF-analyzed pathways
generate_psf_report(psf_list = melanoma_psf, folder_name = "example_psf_report", 
                    plot_type = "kegg", log_norm = TRUE, use_old_images = TRUE)
```

## **Running PSF on Spatial Transcriptomics Data**

> **Important:** Due to installation issues with the `Seurat` package, it is included in `Suggests` rather than `Imports`. The PSF package will still install even if `Seurat` fails to install. However, you must manually install `Seurat` to perform spatial transcriptomic analyses.

```r
library(psf)
library(Seurat)

# Load the spatial dataset. The human melanoma dataset can be downloaded from the following link: https://www.10xgenomics.com/datasets/human-melanoma-if-stained-ffpe-2-standard. To import the data into Seurat, please download the Feature / barcode matrix HDF5 (filtered) and the Spatial imaging data files.
spatial_melanoma <- Load10X_Spatial("/path/to/spatial_dataset", 
                                    filename = "CytAssist_FFPE_Human_Skin_Melanoma_filtered_feature_bc_matrix.h5")

# Load pre-downloaded gene symbol to Entrez ID conversion data
load(system.file("extdata", "gene_symbol_to_entrez.RData", package="psf"))

# Load curated KEGG signaling pathways
load(system.file("extdata", "kegg_curated_40_signalings.RData", package="psf"))
load(system.file("extdata", "kegg_sink_to_process.RData", package="psf"))

# Run pathway analysis
psf_spatial <- spatial_psf_analysis(spatial_obj = spatial_melanoma, 
                                    pathway_collection = kegg_curated_40_signalings, 
                                    gene_symbol_to_entrez = gene_symbol_to_entrez, nthreads = 1)

# Launch the PSF spatial browser application
run_psf_spatial_browser(psf_spatial_results = psf_spatial)
```

### App Usage

A short demo video demonstrating the PSF spatial browser is available: [PSF Spatial Browser Demo](https://www.youtube.com/watch?v=lHTgYBA374o)

You can also explore a **demo Shiny app**: [PSF Spatial Browser](https://apps.armlifebank.am/PSF_spatial_browser/)

The app includes a **dedicated Help tab** containing detailed instructions on how to use it.

## **Running PSF on Telomere Maintenance Mechanism (TMM) Pathways**

```r
library(psf)

# Load melanoma DESeq2 normalized expression data
load(system.file("extdata", "melanoma_exp.RData", package="psf"))

# Load TMM pathway
load(system.file("extdata", "tmm_pathway.RData", package="psf"))

# Calculate fold change (FC) values against the global mean
melanoma_fc <- (melanoma_deseq_normalized_counts + 1) / rowMeans(melanoma_deseq_normalized_counts + 1)

# Run PSF analysis on the TMM pathway for all samples
psf_tmm_result_mel <- psf.from.env.entrez.fc(entrez.fc = melanoma_fc,
                                             kegg.collection = tmm_pathway, 
                                             calculate.significance = FALSE, sum = FALSE, split = FALSE, 
                                             return_only_signals = FALSE, tmm_mode = TRUE, tmm_updated_mode = TRUE)

# Plot the TMM pathway for a sample
plot_tmm_pathway(pathway = psf_tmm_result_mel$X031$tmm, no_color_mode = FALSE,
                 mapping_data_type = "signal", log_norm = FALSE, layout = "layout_nicely")
```

For more details, please refer to the package documentation and tutorials.