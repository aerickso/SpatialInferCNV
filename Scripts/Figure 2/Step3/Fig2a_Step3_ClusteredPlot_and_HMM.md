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

\#Pre-processing clustered data

Since we have already previously created the count matrices, we just
need to re-label the annotations to match the clone annotations in step
2, output an updated annotation file, and then use said file as input
(with the original count matrix from step 2) to generate the final runs

    PurestBenigns_All <- read.csv("C:/Users/erick/Dropbox/Spatial Transcriptomics/SPACE_Data/April2019/InferCNV_R/InferCNV_Analyses/Consensus_AllBenigns_02112020/Outputs/Consensus_PurestBenigns_18112020.csv")
    PurestBenigns_All$Histology <- "Purest Benigns"

    CorrectedBenigns_Consensus_AllCancer_ManualNodes_selected <- read.csv("C:/Users/erick/Dropbox/Spatial Transcriptomics/SPACE_Data/April2019/InferCNV_R/GithubRepo_21022022/Figure2a_AllCancer_siCNV/Step 2  - ManualClustering and Fig 2e Visualization of Clones/Fig2a_forclustering_21022022.csv")
    names(CorrectedBenigns_Consensus_AllCancer_ManualNodes_selected)[2] <- "Histology"

    Fig2a_ManualClusters <- rbind(CorrectedBenigns_Consensus_AllCancer_ManualNodes_selected, PurestBenigns_All)

    write.table(Fig2a_ManualClusters, "Fig2a_ManualClustersOrdered_for_ClusteredPlot_and_HMM.tsv", 
                sep = "\t",
                quote = FALSE, 
                col.names = FALSE, 
                row.names = FALSE)

\#Unsupervised Test (prior to pushing to rescomp)

    AllCancer_clustered <- infercnv::CreateInfercnvObject(raw_counts_matrix="./Organscale_Unsupervised_Consensus_AllCancer_Counts_21022022.tsv", 
                                                   gene_order_file="./gene_position_27072020.tsv",
                                                   annotations_file="./Fig2a_ManualClustersOrdered_for_ClusteredPlot_and_HMM.tsv",
                                                   delim="\t",
                                                   ref_group_names="Purest Benigns",
                                                                   chr_exclude = c("chrM"))

    AllCancer_clustered = infercnv::run(AllCancer_clustered,
                                                  cutoff=0.1,
                                                  out_dir="./Outputs", 
                                                            num_threads = 20,
                                                  cluster_by_groups=TRUE, 
                                                  denoise=TRUE,
                                                  HMM=TRUE)
