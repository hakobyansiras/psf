# Add edge to GraphNEL graph

Add edge to GraphNEL graph

## Usage

``` r
add.kegg.edge(entry1, entry2, type, subtype, edge.attrs, g)
```

## Arguments

- entry1:

  interactor 1

- entry2:

  interactor 2

- type:

  KEGG edge relation type (eg. PPrel - protein protein relation, GErel-
  gene expression relation)

- subtype:

  interaction subypte

- edge.attrs:

  edge attributes

- g:

  GraphNEL object

## Value

GraphNEL object containing modified pathway
