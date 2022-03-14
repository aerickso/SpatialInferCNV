Now that we ran the previous unsupervised step, we next will identified
clones and run the final clustered inferCNVs to generate the clustered
figure panel image in 4c.

# Setup

Initializing libraries.

    library(SpatialInferCNV)

    ## Warning: replacing previous import 'phylogram::as.phylo' by 'ape::as.phylo' when
    ## loading 'SpatialInferCNV'

    library(phylogram)
    library(ape)

    ## 
    ## Attaching package: 'ape'

    ## The following object is masked from 'package:phylogram':
    ## 
    ##     as.phylo

    library(tidyverse)

    ## Registered S3 method overwritten by 'cli':
    ##   method     from         
    ##   print.boxx spatstat.geom

    ## -- Attaching packages --------------------------------------- tidyverse 1.3.1 --

    ## v ggplot2 3.3.5     v purrr   0.3.4
    ## v tibble  3.1.1     v dplyr   1.0.6
    ## v tidyr   1.1.3     v stringr 1.4.0
    ## v readr   2.0.1     v forcats 0.5.1

    ## -- Conflicts ------------------------------------------ tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

# Importing dendrogram

Next, we want to import this dendrogram file, this was created in the
previous step.

    SCC_for_clustering <- read.dendrogram(file = "./Figure4c_Step2/Outputs/infercnv.21_denoised.observations_dendrogram.txt")

    SCC_for_clustering_phylo <- as.phylo(SCC_for_clustering)

# Visualizing Tree

Next, we want to visualize the numbers associated with the nodes of
interest (clones). We output a large image file that allows us to
manually inspect which nodes (corresponding to clones) should be
selected.

    my.subtrees = subtrees(SCC_for_clustering_phylo)  # subtrees() to subset

    png("SCC_for_clustering_phylo.png",width=10000,height=2500, res = 300)
    plot(SCC_for_clustering_phylo,show.tip.label = FALSE)
    nodelabels(text=1:SCC_for_clustering_phylo$Nnode,node=1:SCC_for_clustering_phylo$Nnode+Ntip(SCC_for_clustering_phylo))
    dev.off()

We provide the following output image.

![infercnv.21\_denoised.png](https://github.com/aerickso/SpatialInferCNV/blob/main/FigureScripts/Figure%204/Figure4c_SCC/Step3/SCC_for_clustering_phylo.png)

# Clone selection

Next, view the output .png file, which provides a (albeit cluttered)
labeling of the dendrogram tree nodes. Manually select individual nodes
that correspond with a distinct subclonal grouping or signal, that will
be taken forward for re-clustering. This can be iteratively tweaked with
the next step + spatial visualization til optimal. We provide more
details
[here](https://github.com/aerickso/SpatialInferCNV/blob/main/FigureScripts/Figure%203/Figure3.md),
and provide the finalized selected clone nodes here.

We output a Figure4c\_SCC\_P6\_Clones.csv file, identifying the barcodes
and annotations for each clone for the next steps.

    #A - 1656 -  spots
    #B - 1322 -  spots
    #C - 1183 -  spots
    #D - 2  -  spots

    Node1656 <- SelectingSubTreeData(my.subtrees, 1656)
    Node1322 <- SelectingSubTreeData(my.subtrees, 1322)
    Node1183 <- SelectingSubTreeData(my.subtrees, 1183)
    Node2 <- SelectingSubTreeData(my.subtrees, 2)

    Merged <- rbind(Node1656, Node1322)
    Merged <- rbind(Merged, Node1183)
    Merged <- rbind(Merged, Node2)

    table(Merged$Node)

    Merged$Node <- ifelse(Merged$Node == "Node_1656" , "Clone_A", 
                         ifelse(Merged$Node == "Node_1322" , "Clone_B",
                         ifelse(Merged$Node == "Node_1183" , "Clone_C",
                         ifelse(Merged$Node == "Node_2" , "Clone_D",Merged$Node))))

    write.csv(Merged, "Figure4c_SCC_P6_Clones.csv", row.names = FALSE)

# Outputting the requisite files for infercnv::run

We import the files generated in step 2, with the updated clone
barcodes, and generate a new annotation file for input to infercnv::run.

    library(tidyverse)
    library(SpatialInferCNV)

    OriginalBarcodes <- read.table("./SCC_P6_BenignRef_and_Visium_Mapped_Annotations.tsv", sep = "\t")

    ClusteredBarcodes <- read.csv("./Figure4c_SCC_P6_Clones.csv")

    names(OriginalBarcodes)[1] <- "Barcode"
    names(OriginalBarcodes)[2] <- "Histology"

    UpdatedBarcodes <- left_join(OriginalBarcodes, ClusteredBarcodes)

    UpdatedBarcodes$Node <- ifelse(is.na(UpdatedBarcodes$Node), "PurestBenign_SCCPatient6", UpdatedBarcodes$Node)

    UpdatedBarcodes <- UpdatedBarcodes %>%
                          select(Barcode, Node) %>%
                          arrange(desc(Node))

    write.table(UpdatedBarcodes, 
                "Clustered_SCC_P6_BenignRef_and_Visium_Mapped_Annotations.tsv", 
                quote = FALSE, 
                col.names = FALSE, 
                row.names = FALSE, 
                sep = "\t")

# Creating the inferCNV object (prior to run)

We generate the infercnv object.

    SCC_P6_ForClusteringClones <- infercnv::CreateInfercnvObject(raw_counts_matrix="./SCC_P6_BenignRef_and_Visium_Mapped_Counts.tsv", 
                                                                    gene_order_file="./SCC_P6_BenignRef_and_Visium_GeneOrderFile.tsv",
                                                                    annotations_file="./Clustered_SCC_P6_BenignRef_and_Visium_Mapped_Annotations.tsv",
                                                                    delim="\t",
                                                                    ref_group_names="PurestBenign_SCCPatient6",
                                                                                    chr_exclude = c("chrM"))

# InferCNV Run - (Typically ran on cluster)

Running infercnv.

    SCC_P6_ForClusteringClones = infercnv::run(SCC_P6_ForClusteringClones,
                                                  cutoff=0.1,
                                                  out_dir="./Figure4c_Step3/Outputs", 
                                                  cluster_by_groups=TRUE,
                                                  num_threads = 20, 
                                                  denoise=TRUE,
                                                  HMM=TRUE)

InferCNV will output many files. We are primarily interested in the
final “infercnv.21\_denoised.png” file, corresponding to the one
provided in Figure 4c. These are reordered in the final figure.

![infercnv.21\_denoised.png](https://github.com/aerickso/SpatialInferCNV/blob/main/FigureScripts/Figure%204/Figure4c_SCC/Step3/infercnv.21_denoised.png)
