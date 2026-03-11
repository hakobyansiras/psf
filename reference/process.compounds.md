# Process protein compound interactions

Process protein compound interactions

## Usage

``` r
process.compounds(g, edge.attrs)
```

## Arguments

- g:

  The pathway graph of graphNEL class

- edge.attrs:

  The list of edge attributes in the pathway graph

## Details

This function will turn gene -\> gene interactions that occur via a
compound node into gene -\> compound -\> gene interactions.
