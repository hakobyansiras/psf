# Redirect the edge to new source and target nodes

Redirect the edge to new source and target nodes

## Usage

``` r
redirectEdge(prevSource, prevTarget, newSource, newTarget, g, edge.attrs)
```

## Arguments

- prevSource:

  The original source node

- prevTarget:

  The original target node

- newSource:

  The new source node

- newTarget:

  The new target node

- g:

  The pathway graph of graphNEL class

- edge.attrs:

  The list of edge attributes in the graph

## Details

A new edge between the new source and target nodes will be added to the
graph, the previous edge will be removed, and its attributes will be
transfered to the new edge.
