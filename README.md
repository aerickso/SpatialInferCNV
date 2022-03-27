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

This code was tested using [R](https://www.r-project.org/) version
4.1.3, a Windows 11 Computer, 32GB RAM, and 12 CPUs (1.6 GHz).

For timely data-analyses of datasets comprising 2 or more Visium
sections, consider use of a high performance cluster. In our project,
the infercnv::run analyses steps were ran on the
[BMRC](https://www.medsci.ox.ac.uk/divisional-services/support-services-1/bmrc/cluster-usage),
with 10-20 CPUs, each 1.6 GHz and 16GB ram.

# Installation of SpatialInferCNV Dependencies

``` r
install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE))
     install.packages("BiocManager")
BiocManager::install("infercnv")
install.packages("tidyverse")
install.packages("Seurat")
install.packages("phylogram")
install.packages("ape")
install.packages("hdf5r")
```

# Installation - Dev

``` r
#Note: upon manuscript acceptance, the package will be made public, and thus the need for auth tokens will be removed and thus this code chunk will be deleted. For collaborators, you may need to generate your own new auth token.

AUTH = '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'
install.packages("devtools")
library(devtools)
install_github('aerickso/SpatialInferCNV',
                         auth_token = AUTH)
library(SpatialInferCNV)
```

# Installation - Public

``` r
#This will be the final code chunk, with the above chunk (Installation - Dev) deleted upon release
install.packages("devtools")
library(devtools)
install_github("aerickso/SpatialInferCNV")
library(SpatialInferCNV)
```

# Userguide

The package provides a number of functions, please read the function
documentation in each function located
[here](https://github.com/aerickso/SpatialInferCNV/blob/main/UserGuide/UserGuideDraft.md).

# Study Data

We provide data used in this study at the following [Mendeley
Repository](https://data.mendeley.com/v1/datasets/svw96g68dv/draft?a=3f263217-2bd3-4a3c-8125-8c517c3a9e29).
