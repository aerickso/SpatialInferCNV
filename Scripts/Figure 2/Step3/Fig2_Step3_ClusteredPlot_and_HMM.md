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

Since we have already previously created the count matrices, we just
need to re-label the annotations to match the clone annotations in step
2, output an updated annotation file, and then use said file as input
(with the original count matrix from step 2) to generate the final runs

    PurestBenigns_All <- read.csv("./Mendeley/ProcessedFilesForFigures/Figure2/Step1/Inputs/Consensus_PurestBenigns_18112020.csv")
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

    AllCancer_clustered <- infercnv::CreateInfercnvObject(raw_counts_matrix="./Organscale_Unsupervised_Consensus_AllCancer_Counts.tsv", 
                                                   gene_order_file="./gene_position_27072020.tsv",
                                                   annotations_file="./Fig2_ManualClusters_for_ClusteredPlot_and_HMM.tsv",
                                                   delim="\t",
                                                   ref_group_names="Purest Benigns",
                                                                   chr_exclude = c("chrM"))

# Unsupervised Run - (Typically ran on cluster)

    AllCancer_clustered = infercnv::run(AllCancer_clustered,
                                                  cutoff=0.1,
                                                  out_dir="./Figure2_Step3/Outputs", 
                                                            num_threads = 20,
                                                  cluster_by_groups=TRUE, 
                                                  denoise=TRUE,
                                                  HMM=TRUE)
