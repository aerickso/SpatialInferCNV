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

\#Importing dendrogram

    SCC_for_clustering <- read.dendrogram(file = "./infercnv.21_denoised.observations_dendrogram.txt")

    SCC_for_clustering_phylo <- as.phylo(SCC_for_clustering)

# Visualizing Tree - Option 3

    my.subtrees = subtrees(SCC_for_clustering_phylo)  # subtrees() to subset

    png("SCC_for_clustering_phylo.png",width=10000,height=2500, res = 300)
    plot(SCC_for_clustering_phylo,show.tip.label = FALSE)
    nodelabels(text=1:SCC_for_clustering_phylo$Nnode,node=1:SCC_for_clustering_phylo$Nnode+Ntip(SCC_for_clustering_phylo))
    dev.off()

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
