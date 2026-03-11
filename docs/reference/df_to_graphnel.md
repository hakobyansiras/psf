# Converts 2 data frames(node_table, edge_table) to graphNEL object

Converts 2 data frames(node_table, edge_table) to graphNEL object

## Usage

``` r
df_to_graphnel(node_table, edge_table)
```

## Arguments

- node_table:

  node data frame with specified columns. Mandatory columns: node_id,
  label, component_id_s, type. Non mandatory columns: psf_function, x,
  y.

- edge_table:

  edge data frame with specified columns. Mandatory columns: from, to,
  type. Non mandatory columns subtype, state, weight.
