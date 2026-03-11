# Download kegg pathways of provided pathway id list from keggrest and generate kegg collection

Download kegg pathways of provided pathway id list from keggrest and
generate kegg collection

## Usage

``` r
generate.kegg.collection(pathway.id.list, out.dir, sink.nodes = T)
```

## Arguments

- pathway.id.list:

  charachter vector of pathway ids (hsa04151)

- out.dir:

  path to the directory where kgmls will be downloded

- sink.nodes:

  determin sink nodes, logical TRUE or FALSE
