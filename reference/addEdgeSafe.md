# Add an edge between two nodes if no such edge exists, and if no reverse edge exists

Add an edge between two nodes if no such edge exists, and if no reverse
edge exists

## Usage

``` r
addEdgeSafe(sourceNode, targetNode, g)
```

## Arguments

- sourceNode:

  The source node

- targetNode:

  The target node

- g:

  The graph of class graphNEL

## Value

g The modified graph with added (or not) edge
