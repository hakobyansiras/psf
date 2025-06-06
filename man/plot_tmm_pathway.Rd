% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fast_psf_functions.R
\name{plot_tmm_pathway}
\alias{plot_tmm_pathway}
\title{Plots the TMM pathway with colored nodes and labels with interactive network.}
\usage{
plot_tmm_pathway(
  pathway,
  no_color_mode = F,
  mapping_data_type = "signal",
  log_norm = TRUE,
  layout = "layout_nicely"
)
}
\arguments{
\item{pathway}{pathway object.}

\item{no_color_mode}{whetherto colorcode nodes based on node values or not. Default value is FALSE.}

\item{mapping_data_type}{type of node values to be visualized. When value type is specified pathway nodes will be color coded with expression FC values or PSF values and color legend will be added to the pathway plot. Possible values are c(NULL, "signal", "exp_fc"). Default value is signal.}

\item{log_norm}{log transform PSF and expression values before color mapping. Default value is TRUE.}

\item{layout}{layout type of the pathway. Default value is layout_nicely. Available values c("layout_nicely").}
}
\description{
Plots the TMM pathway with colored nodes and labels with interactive network.
}
