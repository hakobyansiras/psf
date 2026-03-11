# Package index

## Pathway curation (Shiny)

- [`run_shiny_app()`](run_shiny_app.md) : Run KEGG interactive editor
  and visualization App
- [`pathway_shiny_vis()`](pathway_shiny_vis.md) : Run KEGG interactive
  visualization App for provided pathway
- [`generate.kegg.collection()`](generate.kegg.collection.md) : Download
  kegg pathways of provided pathway id list from keggrest and generate
  kegg collection
- [`generate.kegg.collection.from.kgml()`](generate.kegg.collection.from.kgml.md)
  : Generate kegg collection from kgml files
- [`download.KGML()`](download.KGML.md) : Download a KEGG pathway from
  the web
- [`parse.KGML()`](parse.KGML.md) : Parse KGML to a graph object
- [`get.pathway.attrs()`](get.pathway.attrs.md) : Get pathway general
  attributes

## Analysis (Shiny and R)

- [`run_psf()`](run_psf.md) : Calculates pathway activity with PSF
  algorithm for provided kegg collection based on expression fold change
  data.
- [`psf.from.env.entrez.fc()`](psf.from.env.entrez.fc.md) : Calculates
  pathway activity with PSF algorithm for provided kegg collection based
  on expression fold change data.
- [`set.edge.impacts()`](set.edge.impacts.md) : Sets edge impacts for
  PSF analysis
- [`determine.sink.nodes()`](determine.sink.nodes.md) : Detects terminal
  nodes of the pathway
- [`determine.input.nodes()`](determine.input.nodes.md) : Returns vector
  of node ids which do not have incoming edges but only outgoing.
- [`order.nodes()`](order.nodes.md) : Order graph node
- [`map.gene.data()`](map.gene.data.md) : Map gene data onto a pathway
  graph
- [`plot_pathway()`](plot_pathway.md) : Plots the pathway with colored
  nodes and labels
- [`plot_kegg_image_pathway()`](plot_kegg_image_pathway.md) : Plots the
  pathway with colored nodes and labels
- [`plot_tmm_pathway()`](plot_tmm_pathway.md) : Plots the TMM pathway
  with colored nodes and labels with interactive network.

## Reporting

- [`generate_psf_report()`](generate_psf_report.md) : Generates pdf
  report with colored pathways and plots
- [`calc_psf_and_generate_report_from_collection()`](calc_psf_and_generate_report_from_collection.md)
  : Calculates psf for given kegg pathway based on expression matrix and
  generates pdf report with colored pathways and plots

## Utilities (partial influence and graph)

- [`run_pi()`](run_pi.md) : Performs partial influence analysis which
  evaluates effect of each pathway node on specific node(s) of the
  pathway.
- [`calc_node_partial_influences()`](calc_node_partial_influences.md) :
  Returns ordered list of the nodes by their influence on the signal of
  specified nodes.
- [`graphnel_to_df()`](graphnel_to_df.md) : Converts graphNEL object to
  2 data frames(node_table, edge_table)
- [`df_to_graphnel()`](df_to_graphnel.md) : Converts 2 data
  frames(node_table, edge_table) to graphNEL object
- [`edge_data_frame_from_graph()`](edge_data_frame_from_graph.md) :
  Export data frame from graphNEL graph for edge data and its attributes
- [`update_edge_weights()`](update_edge_weights.md) : Import edge
  weights extracted(further edited) via edge_data_frame_from_graph
  function
- [`correctEdgeDirections()`](correctEdgeDirections.md) : Guess wrong
  directed binding interactions and reverse them
- [`isReverseDirection()`](isReverseDirection.md) : This function
  predicts if the edge diractions are wrong based on graphical position
  of the KEGG nodes
- [`out.edges()`](out.edges.md) : provide outgoing edges of the
  specified node
- [`edge.exists()`](edge.exists.md) : Check if the edge exists in the
  graph
- [`get.edge.type()`](get.edge.type.md) : Returns the general edge type
  (either activation or inhibition)
- [`process.compounds()`](process.compounds.md) : Process protein
  compound interactions
- [`process.groupNode()`](process.groupNode.md) : Extend the group node
  to its component gene nodes
- [`redirectEdge()`](redirectEdge.md) : Redirect the edge to new source
  and target nodes
- [`reverseEdge()`](reverseEdge.md) : Reverse edge direction
- [`remove.disconnected.nodes()`](remove.disconnected.nodes.md) : Remove
  nodes which do not have any interactions with other nodes
- [`add.kegg.edge()`](add.kegg.edge.md) : Add edge to GraphNEL graph
- [`add.kegg.edge.mut()`](add.kegg.edge.mut.md) : Add edge to GraphNEL
  graph
- [`addEdgeSafe()`](addEdgeSafe.md) : Add an edge between two nodes if
  no such edge exists, and if no reverse edge exists

## Use cases (Spatial and TMM)

- [`spatial_psf_analysis()`](spatial_psf_analysis.md) : Performs pathway
  activity analysis of Spatial transcriptomics data and subsequent
  clustering with Seurat clustering. Input data is a Seurat object.
- [`run_psf_spatial_browser()`](run_psf_spatial_browser.md) : Run
  spatial PSF browser App
- [`interactive_spatial_plot()`](interactive_spatial_plot.md) : Renders
  interactive plot of spatial tissue slice.
- [`plot_tmm_pathway()`](plot_tmm_pathway.md) : Plots the TMM pathway
  with colored nodes and labels with interactive network.
