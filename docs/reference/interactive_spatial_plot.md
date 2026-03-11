# Renders interactive plot of spatial tissue slice.

Renders interactive plot of spatial tissue slice.

## Usage

``` r
interactive_spatial_plot(
  seurat_object,
  feature = NULL,
  group = "seurat_clusters",
  log_norm = T
)
```

## Arguments

- seurat_object:

  Seurat object of spatial transcriptomics dataset.

- feature:

  If feature is provided spots of the tissue slice will be color coded
  by feature values.

- group:

  group variable from meta.data object based on which spots will be
  colored. Default values is "seurat_clusters".

- log_norm:

  log transform spot values before color mapping. Default value is TRUE.
