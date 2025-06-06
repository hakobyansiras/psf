% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fast_psf_functions.R
\name{run_psf}
\alias{run_psf}
\title{Calculates pathway activity with PSF algorithm for provided kegg collection based on expression fold change data.}
\usage{
run_psf(
  entrez.fc,
  kegg.collection,
  calculate.significance = F,
  bst.steps = 200,
  shuffling_mode = "global",
  sum = FALSE,
  split = TRUE,
  tmm_mode = FALSE,
  tmm_update_mode = FALSE,
  ncores = 1
)
}
\arguments{
\item{entrez.fc}{expression fold change matrix with gene etrez id rownames.}

\item{kegg.collection}{the list of kegg pathways generated by generate.kegg.collection.from.kgml or generate.kegg.collection function.}

\item{calculate.significance}{logical, if true then function will also calculate significance for the PSF values by shuffling all the network nodes and checking if the resulted PSF values were calculated by chance.}

\item{bst.steps}{integer, the number of the iterations for shuffling and recalculating PSF values for the significance analysis.}

\item{shuffling_mode}{character, default value is "global". With this argument set of genes is selected for bootstrapping. When set to "global" all provided FC values will be used to generate random FC values for the pathway. When set to "local" FC values of pathway genes will be shuffled each time and used.}

\item{sum}{logical, default value is FALSE. When set to true pathway activity will be calculated via addition, when set to false then activity will be calculated via multiplication.}

\item{split}{logical, if true then the incoming signal will be proportionally spitted among the edges.}

\item{tmm_mode}{logical, default value is FALSE. When set to true specific PSF configuration will be used for the pathway activity calculation described in https://www.frontiersin.org/articles/10.3389/fgene.2021.662464/full}

\item{tmm_update_mode}{logical, default value is FALSE. When set to true specific PSF configuration will be used for the pathway activity calculation described in https://www.frontiersin.org/articles/10.3389/fgene.2021.662464/full}

\item{ncores}{integer, the number of CPUs to use for PSF calculation. Default value is 1.}
}
\value{
Kegg collection list with evaluated PSF formulas, calculated PSF values and significance values. Each object of the list corresponds to one pathway where "psf_activities", "p_val_mat", "exp_fc" matrices are added (each column of the matrix is one sample). This list object can be further used for pathway visualization and PSF report generation.
}
\description{
Calculates pathway activity with PSF algorithm for provided kegg collection based on expression fold change data.
}
