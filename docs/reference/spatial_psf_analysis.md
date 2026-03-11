# Performs pathway activity analysis of Spatial transcriptomics data and subsequent clustering with Seurat clustering. Input data is a Seurat object.

Performs pathway activity analysis of Spatial transcriptomics data and
subsequent clustering with Seurat clustering. Input data is a Seurat
object.

## Usage

``` r
spatial_psf_analysis(
  spatial_obj,
  pathway_collection,
  gene_symbol_to_entrez,
  nthreads = 1,
  return_only_shiny_vars = TRUE
)
```

## Arguments

- spatial_obj:

  Seurat object of spatial transcriptomic dataset.

- pathway_collection:

  pathway collection list which will be used of PSF analysis

- gene_symbol_to_entrez:

  gene id conversion vector where names are the gene identifiers
  corresponding to count matrix row names and objects are entrez gene
  ids.

- nthreads:

  integer, the number of CPUs to use for calculation. Default value is
  1.

- return_only_shiny_vars:

  when set tot TRUE list of variables will be return which can be used
  to lunch PSF spatial browser shiny app. Default value is TRUE.
