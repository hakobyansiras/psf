# Plots the TMM pathway with colored nodes and labels with interactive network.

Plots the TMM pathway with colored nodes and labels with interactive
network.

## Usage

``` r
plot_tmm_pathway(
  pathway,
  no_color_mode = F,
  mapping_data_type = "signal",
  log_norm = TRUE,
  layout = "layout_nicely"
)
```

## Arguments

- pathway:

  pathway object.

- no_color_mode:

  whetherto colorcode nodes based on node values or not. Default value
  is FALSE.

- mapping_data_type:

  type of node values to be visualized. When value type is specified
  pathway nodes will be color coded with expression FC values or PSF
  values and color legend will be added to the pathway plot. Possible
  values are c(NULL, "signal", "exp_fc"). Default value is signal.

- log_norm:

  log transform PSF and expression values before color mapping. Default
  value is TRUE.

- layout:

  layout type of the pathway. Default value is layout_nicely. Available
  values c("layout_nicely").
