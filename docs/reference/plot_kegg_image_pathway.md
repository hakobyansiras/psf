# Plots the pathway with colored nodes and labels

Plots the pathway with colored nodes and labels

## Usage

``` r
plot_kegg_image_pathway(
  pathway,
  no_color_mode = T,
  mapping_data_type = "signal",
  log_norm = TRUE,
  use_old_images = FALSE,
  plot_type = "kegg",
  highlight_nodes = NULL,
  highlight_color = "red",
  adj = c(0.48, 1),
  y_adj_text = 0,
  y_adj_sink = 0
)
```

## Arguments

- pathway:

  list object from kegg collection

- no_color_mode:

  when set to FALSE pathway nodes will be color coded with log
  expression FC values or log PSF values and color legend will be added
  to the pathway plot. Default value is TRUE

- mapping_data_type:

  type on node values to bi visualized. Possible values c("signal",
  "exp").

- log_norm:

  log transform PSF and expression values before color mapping. Default
  value is TRUE

- use_old_images:

  use_old_images use olde kegg images(for use with curated pathway
  collection).

- plot_type:

  network visualization type. Possible values c("kegg", "visnet"). When
  kegg option is used the bathway will be plotted over kegg png image.
  With visnet option function will plot interactive network with kegg
  layou.

- highlight_nodes:

  single value of node id or a vector of ids to be highlighted in
  plotted pathway. Default value is NULL

- highlight_color:

  Highlighed nodes color(s). Default values is "red"

- adj:

  two values in between 0 and 1 which specify the x and y adjustment of
  the node labels, with 0 for left/bottom, 1 for right/top, and 0.5 for
  centered. Default value is c(0.48, 1). On most devices values outside
  0 and 1 will also work. Only applicable for KEGG visualization type.

- y_adj_text:

  Numeric value to adjust Y position of node labels. Depending on R
  version this value should be adjusted to have accurate text alignment.
  Default value is 0.

- y_adj_sink:

  Numeric value to adjust Y position of sink node sign labels. Depending
  on R version this value should be adjusted to have accurate alignment.
  Default value is 0.
