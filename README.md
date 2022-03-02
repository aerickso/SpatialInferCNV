<img src="https://www.nds.ox.ac.uk/images/logos/secondary-logo" height="75" /> <img src="https://www.nds.ox.ac.uk/images/logos/primary-logo" height="75"/> 

# siCNV: Spatial InferCNV from Spatial Transcriptomics Data

Spatially resolved transcriptomics has emerged as a genome-wide analysis
of gene expression to explore tissues in an unsupervised manner. In this
study we infer genome-wide copy-number variations (CNV) from spatially
resolved mRNA profiles in situ. Gene expression has [previously been
used to infer CNVs](https://github.com/broadinstitute/infercnv) in
single cells, successfully identifying regions of chromosomal gain and
loss. Here we expand into a spatial modality, generating CNV calls in
each spatial region represented by barcoded spots.

We provide a R package via this github page, as well as scripts to
reproduce the main figures in the manuscript.

This code was tested using [R](https://www.r-project.org/) version 4.0.1
(2020-06-06), a Windows 10 Computer, 16GB RAM, and 4 CPUs (2.5 GHz). For
timely data-analyses of datasets comprising 2 or more Visium sections,
consider use of a high performance cluster. In our project, such
analyses were ran on the
[BMRC](https://www.medsci.ox.ac.uk/divisional-services/support-services-1/bmrc/cluster-usage),
with 10-20 CPUs, each 1.6 GHz and 16GB ram.

# Installation

``` r
#AUTH = '87a674f0a1c03d8ccc57bf23c5303695ec30b7ee'

#install_github('aerickso/SpatialInferCNV',
#                         auth_token = AUTH)
library(devtools)
install_github("aerickso/SpatialInferCNV")
library(SpatialInferCNV)
```

# Study Data

We provide data used in this study at the following [Mendeley
Repository](https://data.mendeley.com/v1/datasets/svw96g68dv/draft).
