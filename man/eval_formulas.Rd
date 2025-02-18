% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fast_psf_functions.R
\name{eval_formulas}
\alias{eval_formulas}
\title{Calculates PSF formulas for each node in graphNEL object.}
\usage{
eval_formulas(
  g,
  node.ordering,
  sink.nodes,
  split = TRUE,
  sum = FALSE,
  tmm_mode = FALSE,
  tmm_update_mode = FALSE
)
}
\arguments{
\item{g}{graphNEL pathway object.}

\item{node.ordering}{order of nodes calculated with order.nodes function.}

\item{sink.nodes}{list of terminal (sink) nodes calculated with determine.sink.nodes function.}

\item{split}{logical, if true then the incoming signal will be proportionally splitted among the edges.}

\item{sum}{logical, default value is FALSE. When set to true pathway activity formulas will be calculated via addition, when set to false then activity formulas will be calculated via multiplication.}

\item{tmm_mode}{when set to true specific PSF configuration will be used for calculation of the pathway activity formulas described in https://www.frontiersin.org/articles/10.3389/fgene.2021.662464/full}

\item{tmm_update_mode}{when set to true specific PSF configuration will be used for calculation the pathway activity formulas described in https://www.frontiersin.org/articles/10.3389/fgene.2021.662464/full}
}
\description{
Calculates PSF formulas for each node in graphNEL object.
}
