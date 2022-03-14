# Setup

    library(tidyverse)

    ## -- Attaching packages --------------------------------------- tidyverse 1.3.1 --

    ## v ggplot2 3.3.5     v purrr   0.3.4
    ## v tibble  3.1.1     v dplyr   1.0.6
    ## v tidyr   1.1.3     v stringr 1.4.0
    ## v readr   2.0.1     v forcats 0.5.1

    ## -- Conflicts ------------------------------------------ tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

    library(infercnv)
    library(Seurat)

    ## Registered S3 method overwritten by 'spatstat.geom':
    ##   method     from
    ##   print.boxx cli

    ## Attaching SeuratObject

    library(hdf5r)

    ## 
    ## Attaching package: 'hdf5r'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     flatten_df

    library(SpatialInferCNV)

    ## Warning: replacing previous import 'phylogram::as.phylo' by 'ape::as.phylo' when
    ## loading 'SpatialInferCNV'

# Pre-processing clustered data

Importing previously downloaded Consensus\_PurestBenigns.csv (step 1),
and the Fig2\_forclustering.csv file created in step 2. We use this to
create an updated annotation file for infercnv::run.

    PurestBenigns_All <- read.csv("./Figure2_output/Patient 1/Consensus_PurestBenigns.csv")
    PurestBenigns_All$Histology <- "Purest Benigns"

    CorrectedBenigns_Consensus_AllCancer_ManualNodes_selected <- read.csv("./Mendeley/ProcessedFilesForFigures/Figure2/Step3/Inputs/Fig2_forclustering.csv")
    names(CorrectedBenigns_Consensus_AllCancer_ManualNodes_selected)[2] <- "Histology"

    Fig2a_ManualClusters <- rbind(CorrectedBenigns_Consensus_AllCancer_ManualNodes_selected, PurestBenigns_All)

    write.table(Fig2a_ManualClusters, "Fig2_ManualClusters_for_ClusteredPlot_and_HMM.tsv", 
                sep = "\t",
                quote = FALSE, 
                col.names = FALSE, 
                row.names = FALSE)

# Creating the inferCNV object (prior to run)

Now creating the object for the supervised clustered run.

    AllCancer_clustered <- infercnv::CreateInfercnvObject(raw_counts_matrix="./Organscale_Unsupervised_Consensus_AllCancer_Counts.tsv", 
                                                   gene_order_file="./siCNV_GeneOrderFile.tsv",
                                                   annotations_file="./Fig2_ManualClusters_for_ClusteredPlot_and_HMM.tsv",
                                                   delim="\t",
                                                   ref_group_names="Purest Benigns",
                                                                   chr_exclude = c("chrM"))

# Unsupervised Run - (Typically ran on cluster)

Now creating the object for the supervised clustered run. Note: this is
typically run

    AllCancer_clustered = infercnv::run(AllCancer_clustered,
                                                  cutoff=0.1,
                                                  out_dir="./Figure2_output/Figure2_step3/Outputs", 
                                                            num_threads = 20,
                                                  cluster_by_groups=TRUE, 
                                                  denoise=TRUE,
                                                  HMM=TRUE)

And here is the final output file infercnv.21\_denoised.png (order
rearranged in the manuscript figure 2).

![infercnv.21\_denoised.png
output](https://github.com/aerickso/SpatialInferCNV/blob/main/FigureScripts/Figure%202/Step3/infercnv.21_denoised.png)
