# Calculates PSF formulas for each node in graphNEL object.

Calculates PSF formulas for each node in graphNEL object.

## Usage

``` r
eval_formulas(
  g,
  node.ordering,
  sink.nodes,
  split = TRUE,
  sum = FALSE,
  tmm_mode = FALSE,
  tmm_update_mode = FALSE
)
```

## Arguments

- g:

  graphNEL pathway object.

- node.ordering:

  order of nodes calculated with order.nodes function.

- sink.nodes:

  list of terminal (sink) nodes calculated with determine.sink.nodes
  function.

- split:

  logical, if true then the incoming signal will be proportionally
  splitted among the edges.

- sum:

  logical, default value is FALSE. When set to true pathway activity
  formulas will be calculated via addition, when set to false then
  activity formulas will be calculated via multiplication.

- tmm_mode:

  when set to true specific PSF configuration will be used for
  calculation of the pathway activity formulas described in
  https://www.frontiersin.org/articles/10.3389/fgene.2021.662464/full

- tmm_update_mode:

  when set to true specific PSF configuration will be used for
  calculation the pathway activity formulas described in
  https://www.frontiersin.org/articles/10.3389/fgene.2021.662464/full
