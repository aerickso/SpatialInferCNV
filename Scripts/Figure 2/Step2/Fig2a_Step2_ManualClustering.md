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

    Consensus_AllCancer_for_clustering <- read.dendrogram(file = "./infercnv.observations_dendrogram.txt")

    Consensus_AllCancer_for_clustering_phylo <- as.phylo(Consensus_AllCancer_for_clustering)

# Visualizing Tree

    my.subtrees = subtrees(Consensus_AllCancer_for_clustering_phylo)  # subtrees() to subset

    png("Consensus_AllCancer_forclustering_phylo.png",width=10000,height=2500, res = 300)
    plot(Consensus_AllCancer_for_clustering_phylo,show.tip.label = FALSE)
    nodelabels(text=1:Consensus_AllCancer_for_clustering_phylo$Nnode,node=1:Consensus_AllCancer_for_clustering_phylo$Nnode+Ntip(Consensus_AllCancer_for_clustering_phylo))
    dev.off()

    Merged$Node <- ifelse(Merged$Node == "Node_4617" , "Clone_J",
                         ifelse(Merged$Node == "Node_4446" , "Clone_I",
                         ifelse(Merged$Node == "Node_3617" , "Clone_I",
                         ifelse(Merged$Node == "Node_2934" , "Clone_F",
                         ifelse(Merged$Node == "Node_2893" , "Clone_E",
                         ifelse(Merged$Node == "Node_2832" , "Clone_E",
                         ifelse(Merged$Node == "Node_2769" , "Clone_E",
                         ifelse(Merged$Node == "Node_3114" , "Clone_H",
                         ifelse(Merged$Node == "Node_2991" , "Clone_B",
                         ifelse(Merged$Node == "Node_2525" , "Clone_G",
                         ifelse(Merged$Node == "Node_2284" , "Clone_C",
                         ifelse(Merged$Node == "Node_2078" , "Clone_D",
                         ifelse(Merged$Node == "Node_92" , "Clone_K",
                         ifelse(Merged$Node == "Node_3" , "Clone_A", Merged$Node))))))))))))))

    #Clone J - Node 4617
    #Clone I - Node 4446
    #Clone I - Node 3617
    #Clone F - Node 2934
    #Clone E - Node 2893
    #Clone E - Node 2832
    #Clone E - Node 2769
    #Clone H - Node 3114 
    #Clone B - Node 2991
    #Clone G - Node 2525
    #Clone C - Node 2284
    #Clone D - Node 2078
    #Clone K - Node 92 
    #Clone A - Node 3

    Node4617 <- SelectingSubTreeData(my.subtrees, 4617)
    Node4446 <- SelectingSubTreeData(my.subtrees, 4446)
    Node3617 <- SelectingSubTreeData(my.subtrees, 3617)
    Node2934 <- SelectingSubTreeData(my.subtrees, 2934)
    Node2893 <- SelectingSubTreeData(my.subtrees, 2893)
    Node2832 <- SelectingSubTreeData(my.subtrees, 2832)
    Node2769 <- SelectingSubTreeData(my.subtrees, 2769)
    Node3114 <- SelectingSubTreeData(my.subtrees, 3114)
    Node2991 <- SelectingSubTreeData(my.subtrees, 2991)
    Node2525 <- SelectingSubTreeData(my.subtrees, 2525)
    Node2284 <- SelectingSubTreeData(my.subtrees, 2284)
    Node2078 <- SelectingSubTreeData(my.subtrees, 2078)
    Node92 <- SelectingSubTreeData(my.subtrees, 92)
    Node3 <- SelectingSubTreeData(my.subtrees, 3)

    Merged <- rbind(Node4617, Node4446)
    Merged <- rbind(Merged, Node3617)
    Merged <- rbind(Merged, Node2934)
    Merged <- rbind(Merged, Node2893)
    Merged <- rbind(Merged, Node2832)
    Merged <- rbind(Merged, Node2769)
    Merged <- rbind(Merged, Node3114)
    Merged <- rbind(Merged, Node2991)
    Merged <- rbind(Merged, Node2525)
    Merged <- rbind(Merged, Node2284)
    Merged <- rbind(Merged, Node2078)
    Merged <- rbind(Merged, Node92)
    Merged <- rbind(Merged, Node3)

    table(Merged$Node)

    Merged$Node <- ifelse(Merged$Node == "Node_4617" , "Clone_J",
                         ifelse(Merged$Node == "Node_4446" , "Clone_I",
                         ifelse(Merged$Node == "Node_3617" , "Clone_I",
                         ifelse(Merged$Node == "Node_2934" , "Clone_F",
                         ifelse(Merged$Node == "Node_2893" , "Clone_E",
                         ifelse(Merged$Node == "Node_2832" , "Clone_E",
                         ifelse(Merged$Node == "Node_2769" , "Clone_E",
                         ifelse(Merged$Node == "Node_3114" , "Clone_H",
                         ifelse(Merged$Node == "Node_2991" , "Clone_B",
                         ifelse(Merged$Node == "Node_2525" , "Clone_G",
                         ifelse(Merged$Node == "Node_2284" , "Clone_C",
                         ifelse(Merged$Node == "Node_2078" , "Clone_D",
                         ifelse(Merged$Node == "Node_92" , "Clone_K",
                         ifelse(Merged$Node == "Node_3" , "Clone_A", Merged$Node))))))))))))))

    write.csv(Merged, "Fig2a_forclustering_21022022.csv", row.names = FALSE)

    H1_5_Merged <- Merged
    H1_5_Merged <- H1_5_Merged %>% mutate(section = substr(Barcode, 1, 4))
    H1_5_Merged$Barcode <- trimws(substr(H1_5_Merged$Barcode, 6, 100))
    H1_5_Merged$Barcode <- gsub("\\.", "\\-", H1_5_Merged$Barcode)
    H1_5_Clones_ForLoupeBrowser <- filter(H1_5_Merged, section == "H1_5") %>%
                                                select(Barcode, Node)
    write.csv(H1_5_Clones_ForLoupeBrowser, "Fig2e_H1_5_Clones_ForLoupeBrowser.csv", row.names = FALSE)

    H2_5_Merged <- Merged
    H2_5_Merged <- H2_5_Merged %>% mutate(section = substr(Barcode, 1, 4))
    H2_5_Merged$Barcode <- trimws(substr(H2_5_Merged$Barcode, 6, 100))
    H2_5_Merged$Barcode <- gsub("\\.", "\\-", H2_5_Merged$Barcode)
    H2_5_Clones_ForLoupeBrowser <- filter(H2_5_Merged, section == "H2_5") %>%
                                                select(Barcode, Node)
    write.csv(H2_5_Clones_ForLoupeBrowser, "Fig2e_H2_5_Clones_ForLoupeBrowser.csv", row.names = FALSE)

    H1_4_Merged <- Merged
    H1_4_Merged <- H1_4_Merged %>% mutate(section = substr(Barcode, 1, 4))
    H1_4_Merged$Barcode <- trimws(substr(H1_4_Merged$Barcode, 6, 100))
    H1_4_Merged$Barcode <- gsub("\\.", "\\-", H1_4_Merged$Barcode)
    H1_4_Clones_ForLoupeBrowser <- filter(H1_4_Merged, section == "H1_4") %>%
                                                select(Barcode, Node)
    write.csv(H1_4_Clones_ForLoupeBrowser, "Fig2e_H1_4_Clones_ForLoupeBrowser.csv", row.names = FALSE)

    H1_2_Merged <- Merged
    H1_2_Merged <- H1_2_Merged %>% mutate(section = substr(Barcode, 1, 4))
    H1_2_Merged$Barcode <- trimws(substr(H1_2_Merged$Barcode, 6, 100))
    H1_2_Merged$Barcode <- gsub("\\.", "\\-", H1_2_Merged$Barcode)
    H1_2_Clones_ForLoupeBrowser <- filter(H1_2_Merged, section == "H1_2") %>%
                                                select(Barcode, Node)
    write.csv(H1_2_Clones_ForLoupeBrowser, "Fig2e_H1_2_Clones_ForLoupeBrowser.csv", row.names = FALSE)

    H2_1_Merged <- Merged
    H2_1_Merged <- H2_1_Merged %>% mutate(section = substr(Barcode, 1, 4))
    H2_1_Merged$Barcode <- trimws(substr(H2_1_Merged$Barcode, 6, 100))
    H2_1_Merged$Barcode <- gsub("\\.", "\\-", H2_1_Merged$Barcode)
    H2_1_Clones_ForLoupeBrowser <- filter(H2_1_Merged, section == "H2_1") %>%
                                                select(Barcode, Node)
    write.csv(H2_1_Clones_ForLoupeBrowser, "Fig2e_H2_1_Clones_ForLoupeBrowser.csv", row.names = FALSE)

    H2_2_Merged <- Merged
    H2_2_Merged <- H2_2_Merged %>% mutate(section = substr(Barcode, 1, 4))
    H2_2_Merged$Barcode <- trimws(substr(H2_2_Merged$Barcode, 6, 100))
    H2_2_Merged$Barcode <- gsub("\\.", "\\-", H2_2_Merged$Barcode)
    H2_2_Clones_ForLoupeBrowser <- filter(H2_2_Merged, section == "H2_2") %>%
                                                select(Barcode, Node)
    write.csv(H2_2_Clones_ForLoupeBrowser, "Fig2e_H2_2_Clones_ForLoupeBrowser.csv", row.names = FALSE)
