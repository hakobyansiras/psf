# Map gene data onto a pathway graph

Map gene data onto a pathway graph

## Usage

``` r
map.gene.data(g, entrez.fc)
```

## Arguments

- g:

  Tha pathway graph of graphNEL class

- entrez.fc:

  gene expression fold change matrix with entrez gene rownames (A matrix
  data associated with the genes; rownames represent genes; a single
  gene-row may contain one or many data values;)

## Value

graphNEL object with the column-wise averaged gene data kept in the
nodedata attribute of the graph
