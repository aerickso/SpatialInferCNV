# siCNV: Spatial InferCNV from Spatial Transcriptomics Data

Spatially resolved transcriptomics has emerged as a genome-wide analysis
of gene expression to explore tissues in an unsupervised manner. In this
study we infer genome-wide copy-number variations (CNV) from spatially
resolved mRNA profiles in situ. Gene expression has [previously been
used to infer CNVs](https://github.com/broadinstitute/infercnv) in
single cells, successfully identifying regions of chromosomal gain and
loss. Here we expand into a spatial modality, generating CNV calls in
each spatial region represented by barcoded spots.

We provide a R package via this github page, as well as [scripts to
reproduce the main
figures](https://github.com/aerickso/SpatialInferCNV/tree/main/FigureScripts)
in the manuscript.

This code was tested using [R](https://www.r-project.org/) version 4.0.1
(2020-06-06), a Windows 10 Computer, 16GB RAM, and 4 CPUs (2.5 GHz). For
timely data-analyses of datasets comprising 2 or more Visium sections,
consider use of a high performance cluster. In our project, such
analyses were ran on the
[BMRC](https://www.medsci.ox.ac.uk/divisional-services/support-services-1/bmrc/cluster-usage),
with 10-20 CPUs, each 1.6 GHz and 16GB ram.

# Installation of SpatialInferCNV Dependencies

``` r
install.packages("tidyverse")
install.packages("infercnv")
install.packages("Seurat")
install.packages("hdf5r")
install.packages("phylogram")
install.packages("ape")

library("tidyverse")
library("infercnv")
library("Seurat")
library("hdf5r")
library("phylogram")
library("ape")
```

# Installation - Dev

``` r
#Note: upon manuscript acceptance, the package will be made public, and thus the need for auth tokens will be removed and thus this code chunk will be deleted. For collaborators, you may need to generate your own new auth token.

#AUTH = '87a674f0a1c03d8ccc57bf23c5303695ec30b7ee'

install.packages("devtools")
library(devtools)
install_github('aerickso/SpatialInferCNV',
                         auth_token = AUTH)
library(SpatialInferCNV)
```

# Installation - Public

``` r
#This will be the final code chunk, with the above chunk deleted
install.packages("devtools")
library(devtools)
install_github("aerickso/SpatialInferCNV")
library(SpatialInferCNV)
```

# Package Functions

The package provides a number of functions, please read the function
documentation in each function located
[here](https://github.com/aerickso/SpatialInferCNV/tree/main/R).

# Study Data

We provide data used in this study at the following [Mendeley
Repository](https://data.mendeley.com/datasets/svw96g68dv/draft?a=3f263217-2bd3-4a3c-8125-5).
