\#Setup

    library(tidyverse)

    ## -- Attaching packages --------------------------------------- tidyverse 1.3.1 --

    ## v ggplot2 3.3.5     v purrr   0.3.4
    ## v tibble  3.1.1     v dplyr   1.0.6
    ## v tidyr   1.1.3     v stringr 1.4.0
    ## v readr   2.0.1     v forcats 0.5.1

    ## -- Conflicts ------------------------------------------ tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

    library(SpatialInferCNV)

    ## Registered S3 method overwritten by 'spatstat.geom':
    ##   method     from
    ##   print.boxx cli

    ## Warning: replacing previous import 'phylogram::as.phylo' by 'ape::as.phylo' when
    ## loading 'SpatialInferCNV'

We start by creating an empty working directory so that all downloaded
files are organized in one place. Something like:

    dir.create("Figure2_output")
    setwd("Figure2_output")

# Consensus Purest Benigns

    PurestBenigns_All <- read.csv("./Mendeley/ProcessedFilesForFigures/Figure2/Step1/Consensus_PurestBenigns.csv")

# Re-creating All Cancer Annotations

    H1_2_Annotations <- ImportHistologicalAnnotations("H1_2", "./Mendeley/Patient 1/Visium/H1_2/H1_2_Final_Consensus_Annotations.csv")
    H1_2_CancerSpots <- filter(H1_2_Annotations, Histology == "GG1")

    H1_4_Annotations <- ImportHistologicalAnnotations("H1_4", "./Mendeley/Patient 1/Visium/H1_4/H1_4_Final_Consensus_Annotations.csv")
    H1_4_CancerSpots <- filter(H1_4_Annotations, Histology == "GG2" | Histology ==  "GG4 Cribriform")

    H1_5_Annotations <- ImportHistologicalAnnotations("H1_5", "./Mendeley/Patient 1/Visium/H1_5/H1_5_Final_Consensus_Annotations.csv")
    H1_5_CancerSpots <- filter(H1_5_Annotations, Histology ==  "GG4 Cribriform")

    H2_1_Annotations <- ImportHistologicalAnnotations("H2_1", "./Mendeley/Patient 1/Visium/H2_1/H2_1_Final_Consensus_Annotations.csv")
    H2_1_CancerSpots <- filter(H2_1_Annotations, Histology == "GG2" | Histology ==  "GG4")

    H2_2_Annotations <- ImportHistologicalAnnotations("H2_2", "./Mendeley/Patient 1/Visium/H2_2/H2_2_Final_Consensus_Annotations.csv")
    H2_2_CancerSpots <- filter(H2_2_Annotations, Histology == "GG2")

    H2_5_Annotations <- ImportHistologicalAnnotations("H2_5", "./Mendeley/Patient 1/Visium/H2_5/H2_5_Final_Consensus_Annotations.csv")
    H2_5_CancerSpots <- filter(H2_5_Annotations, Histology == "GG4 Cribriform" | Histology == "Transition_State")

    rm(H1_2_Annotations,
       H1_4_Annotations,
       H1_5_Annotations,
       H2_1_Annotations,
       H2_2_Annotations,
       H2_5_Annotations)

    AllCancers <- rbind(H1_2_CancerSpots, H1_4_CancerSpots)
    AllCancers <- rbind(AllCancers, H1_5_CancerSpots)
    AllCancers <- rbind(AllCancers, H2_1_CancerSpots)
    AllCancers <- rbind(AllCancers, H2_2_CancerSpots)
    AllCancers <- rbind(AllCancers, H2_5_CancerSpots)

    names(AllCancers)[2] <- "Histology"

    rm(H1_2_CancerSpots,
       H1_4_CancerSpots,
       H1_5_CancerSpots,
       H2_1_CancerSpots,
       H2_2_CancerSpots,
       H2_5_CancerSpots)

    MergedAll <- rbind(PurestBenigns_All, AllCancers)

    rm(PurestBenigns_All)
    rm(AllCancers)

# Merging Cancer and Benign annotations with the ENSMBLIDs

    H2_1_ENSBMLID_Counts <- ImportCountData("H2_1", "./Mendeley/Patient 1/Visium/H2_1/filtered_feature_bc_matrix.h5")
    H2_1_Joined_Counts <- MergingCountAndAnnotationData("H2_1",MergedAll, H2_1_ENSBMLID_Counts)
    rm(H2_1_ENSBMLID_Counts)
    Counts_joined <- H2_1_Joined_Counts
    rm(H2_1_Joined_Counts)

    H1_5_ENSBMLID_Counts <- ImportCountData("H1_5", "./Mendeley/Patient 1/Visium/H1_5/filtered_feature_bc_matrix.h5")
    H1_5_Joined_Counts <- MergingCountAndAnnotationData("H1_5",MergedAll, H1_5_ENSBMLID_Counts)
    rm(H1_5_ENSBMLID_Counts)
    Counts_joined <- Counts_joined %>% full_join(H1_5_Joined_Counts, by = "Genes")
    rm(H1_5_Joined_Counts)

    H2_2_ENSBMLID_Counts <- ImportCountData("H2_2", "./Mendeley/Patient 1/Visium/H2_2/filtered_feature_bc_matrix.h5")
    H2_2_Joined_Counts <- MergingCountAndAnnotationData("H2_2",MergedAll, H2_2_ENSBMLID_Counts)
    rm(H2_2_ENSBMLID_Counts)
    Counts_joined <- Counts_joined %>% full_join(H2_2_Joined_Counts, by = "Genes")
    rm(H2_2_Joined_Counts)

    H1_2_ENSBMLID_Counts <- ImportCountData("H1_2", "./Mendeley/Patient 1/Visium/H1_2/filtered_feature_bc_matrix.h5")
    H1_2_Joined_Counts <- MergingCountAndAnnotationData("H1_2",MergedAll, H1_2_ENSBMLID_Counts)
    rm(H1_2_ENSBMLID_Counts)
    Counts_joined <- Counts_joined %>% full_join(H1_2_Joined_Counts, by = "Genes")
    rm(H1_2_Joined_Counts)

    H2_5_ENSBMLID_Counts <- ImportCountData("H2_5", "./Mendeley/Patient 1/Visium/H2_5/filtered_feature_bc_matrix.h5")
    H2_5_Joined_Counts <- MergingCountAndAnnotationData("H2_5",MergedAll, H2_5_ENSBMLID_Counts)
    rm(H2_5_ENSBMLID_Counts)
    Counts_joined <- Counts_joined %>% full_join(H2_5_Joined_Counts, by = "Genes")
    rm(H2_5_Joined_Counts)

    H1_4_ENSBMLID_Counts <- ImportCountData("H1_4", "./Mendeley/Patient 1/Visium/H1_4/filtered_feature_bc_matrix.h5")
    H1_4_Joined_Counts <- MergingCountAndAnnotationData("H1_4",MergedAll, H1_4_ENSBMLID_Counts)
    rm(H1_4_ENSBMLID_Counts)
    Counts_joined <- Counts_joined %>% full_join(H1_4_Joined_Counts, by = "Genes")
    rm(H1_4_Joined_Counts)

    V1_2_ENSBMLID_Counts <- ImportCountData("V1_2", "./Mendeley/Patient 1/Visium/V1_2/filtered_feature_bc_matrix.h5")
    V1_2_Joined_Counts <- MergingCountAndAnnotationData("V1_2",MergedAll, V1_2_ENSBMLID_Counts)
    rm(V1_2_ENSBMLID_Counts)
    Counts_joined <- Counts_joined %>% full_join(V1_2_Joined_Counts, by = "Genes")
    rm(V1_2_Joined_Counts)

# Joining all Counts

    Counts_joined <- Counts_joined %>% replace(., is.na(.), 0)
    Counts_joined <- Counts_joined %>% column_to_rownames(., var = "Genes")

    write.table(Counts_joined, "Organscale_Unsupervised_Consensus_AllCancer_Counts.tsv", sep = "\t")

    MergedAll_Final <- FinalAnnotations(MergedAll, Counts_joined)

    write.table(MergedAll_Final, "Organscale_Unsupervised_Consensus_AllCancer_Annotations.tsv", 
                sep = "\t",
                quote = FALSE, 
                col.names = FALSE, 
                row.names = FALSE)

# Creating the inferCNV object (prior to run)

    AllCancer_Unsupervised <- infercnv::CreateInfercnvObject(raw_counts_matrix="./Organscale_Unsupervised_Consensus_AllCancer_Counts.tsv", 
                                                   gene_order_file="./siCNV_GeneOrderFile.tsv",
                                                   annotations_file="./Organscale_Unsupervised_Consensus_AllCancer_Annotations.tsv",
                                                   delim="\t",
                                                   ref_group_names="Purest Benigns",
                                                                   chr_exclude = c("chrM"))

# Unsupervised Run - (Typically ran on cluster)

    AllCancer_Unsupervised = infercnv::run(AllCancer_Unsupervised,
                                                  cutoff=0.1,
                                                  out_dir="./Figure2_Step1/Outputs", 
                                                  cluster_by_groups=FALSE,
                                                  num_threads = 20, 
                                                  denoise=TRUE,
                                                  HMM=FALSE)
