# Generates pdf report with colored pathways and plots

Generates pdf report with colored pathways and plots

## Usage

``` r
generate_psf_report(
  psf_list,
  folder_name,
  plot_type = "kegg",
  log_norm = TRUE,
  y_adj_text = 0,
  y_adj_sink = 0,
  use_old_images = F
)
```

## Arguments

- psf_list:

  Output from run_psf function

- folder_name:

  name of the folder where pdf report(s) will be generated

- plot_type:

  Pathway visualization type. Possible values are c("kegg", "visnet").
  Default value is "kegg".

- log_norm:

  log transform PSF values before color mapping. Default value is TRUE.

- y_adj_text:

  Numeric value to adjust Y position of node labels. Depending on R
  version this value should be adjusted to have accurate text alignment.
  Default value is 0.

- y_adj_sink:

  Numeric value to adjust Y position of sink node sign labels. Depending
  on R version this value should be adjusted to have accurate alignment.
  Default value is 0.

- use_old_images:

  use olde kegg images(for use with curated pathway collection)
