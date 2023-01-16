# PSF toolkit

<!-- badges: start -->
[![R-CMD-check](https://github.com/hakobyansiras/psf/workflows/R-CMD-check/badge.svg)](https://github.com/hakobyansiras/psf/actions)
<!-- badges: end -->

Pathway editor and signal flow calculator

# Installation

To install the package run this script in R
```
if (!require("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github("hakobyansiras/psf")
```

# Running PSF
This is a small example script to run PSF analysis on melanoma gene expression data.
```
library(psf)

### loading Melanoma deseq normalized expression data (GSE112509).
load(system.file("extdata", "melanoma_exp.RData", package="psf"))
### loading curated KEGG pathway collection.
load(system.file("extdata", "edited_pathways_new.RData", package="psf"))

### Calculating FC values against global mean.
melanoma_fc <- (melanoma_deseq_normalized_counts + 1)/rowMeans(melanoma_deseq_normalized_counts)

### runing PSF analysis on 36 pathways for 80 samples with 4 cores.
melanoma_psf <- run_psf(entrez.fc = melanoma_fc, kegg.collection = edited_pathways_new, 
                        calculate.significance = F, ncores = 4)

### plotting KEGG image based PI3K_Akt_signaling_pathway with colorcoded nodes. 
plot_pathway(melanoma_psf$PI3K_Akt_signaling_pathway, plot_type = "kegg", 
             color_nodes = "psf_activities", sample_id = "mean", log_norm = T, 
             plot_sink_values = TRUE, use_old_images = T)

### generating PDF reports for PSF analyzed pathways.
generate_psf_report(psf_list = melanoma_psf, folder_name = "example_psf_report", 
                    plot_type = "kegg", log_norm = T, use_old_images = T)
```

# Running interactive application
```run_shiny_app()```