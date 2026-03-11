# Returns ordered list of the nodes by their influence on the signal of specified nodes.

Returns ordered list of the nodes by their influence on the signal of
specified nodes.

## Usage

``` r
calc_node_partial_influences(
  pathway,
  exp_fc_matrix = NULL,
  influenced_node,
  influence_direction = "any",
  node_combinations = 1,
  nproc = 1,
  get_influence_matrix = FALSE
)
```

## Arguments

- pathway:

  Pathway onbject.

- exp_fc_matrix:

  expression fold change matrix with gene etrez id rownames.

- influenced_node:

  sigle id of node or vector of node ids based on which partial
  influencse will be calculated.

- influence_direction:

  in which direction node affects the signal of target node. Possible
  values c("+", "-", "any"). Default value is "any".

- node_combinations:

  number of node combinations for partial influence analysis. Note:
  number of PSF calculations increases by the exponent of combinations(N
  nodes^combinations). Default value is 1.

- nproc:

  number of cpus to use for calculation. Default value is 1.

- get_influence_matrix:

  When set to true influence matrix will be returned instead of
  influencial nodes. Each column in influence matrix represent the log
  ratio of default and neutralized node(s) psf profile.
