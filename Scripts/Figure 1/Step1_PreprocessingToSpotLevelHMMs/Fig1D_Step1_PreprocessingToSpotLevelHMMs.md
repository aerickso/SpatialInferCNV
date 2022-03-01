# Setup

    library(tidyverse)
    library(infercnv)
    library(Seurat)
    library(hdf5r)
    library(readxl)
    library(SpatialInferCNV)

    ## Warning: replacing previous import 'phylogram::as.phylo' by 'ape::as.phylo' when
    ## loading 'SpatialInferCNV'

# Importing Patient 1, 1k array counts

    Counts_H2.1 <- ImportOriginalSTCountData("H2_1", "./Patient 1/1k_arrays/H2_1/180903_L11_CN63_D1_H2.1_EB_stdata.tsv")
    Counts_H2.2 <- ImportOriginalSTCountData("H2_2", "./Patient 1/1k_arrays/H2_2/180614_L11_53_C1_H2.2_EB_stdata.tsv")
    Counts_H2.3 <- ImportOriginalSTCountData("H2_3", "./Patient 1/1k_arrays/H2_3/180614_L11_53_D1_H2.3_EB_stdata.tsv")
    Counts_H2.4 <- ImportOriginalSTCountData("H2_4", "./Patient 1/1k_arrays/H2_4/180614_L11_53_E1_H2.4_EB_stdata.tsv")
    Counts_H2.5 <- ImportOriginalSTCountData("H2_5", "./Patient 1/1k_arrays/H2_5/180614_L11_55_C1_H2.5_EB_stdata.tsv")

    Counts_V2.1 <- ImportOriginalSTCountData("V2_1", "./Patient 1/1k_arrays/V2_1/180614_L11_55_D1_V2.1_EB_stdata.tsv")
    Counts_V2.2 <- ImportOriginalSTCountData("V2_2", "./Patient 1/1k_arrays/V2_2/180614_L11_55_E1_V2.2_EB_stdata.tsv")
    Counts_V2.3 <- ImportOriginalSTCountData("V2_3", "./Patient 1/1k_arrays/V2_3/180813_L11_CN60_C2_V2.3_EB_stdata.tsv")
    Counts_V2.4 <- ImportOriginalSTCountData("V2_4", "./Patient 1/1k_arrays/V2_4/180813_L11_CN60_D2_V2.4_EB_stdata.tsv")
    Counts_V2.5 <- ImportOriginalSTCountData("V2_5", "./Patient 1/1k_arrays/V2_5/180822_L11_CN57_C1_V2.5_EB_stdata.tsv")
    Counts_V2.6 <- ImportOriginalSTCountData("V2_6", "./Patient 1/1k_arrays/V2_6/180813_L11_CN15_C2_V2.6_EB_stdata.tsv")

    Counts_V1.1 <- ImportOriginalSTCountData("V1_1", "./Patient 1/1k_arrays/V1_1/180813_L11_CN15_D1_V1.1_EB_stdata.tsv")
    Counts_V1.2 <- ImportOriginalSTCountData("V1_2", "./Patient 1/1k_arrays/V1_2/180813_L11_CN15_E1_V1.2_EB_stdata.tsv")
    Counts_V1.3 <- ImportOriginalSTCountData("V1_3", "./Patient 1/1k_arrays/V1_3/180822_L11_CN56_C1_V1.3_EB_stdata.tsv")
    Counts_V1.4 <- ImportOriginalSTCountData("V1_4", "./Patient 1/1k_arrays/V1_4/180822_L11_CN56_D2_V1.4_EB_stdata.tsv")
    Counts_V1.5 <- ImportOriginalSTCountData("V1_5", "./Patient 1/1k_arrays/V1_5/180822_L11_CN56_E2_V1.5_EB_stdata.tsv")

    Counts_H1.1 <- ImportOriginalSTCountData("H1_1", "./Patient 1/1k_arrays/H1_1/180822_L11_CN57_E1_H1.1_EB_stdata.tsv")
    Counts_H1.2 <- ImportOriginalSTCountData("H1_2", "./Patient 1/1k_arrays/H1_2/180903_L11_CN51_C2_H1.2_EB_stdata.tsv")
    Counts_H1.3 <- ImportOriginalSTCountData("H1_3", "./Patient 1/1k_arrays/H1_3/180903_L11_CN61_D2_H1.3_EB_stdata.tsv")
    Counts_H1.4 <- ImportOriginalSTCountData("H1_4", "./Patient 1/1k_arrays/H1_4/180903_L11_CN51_E1_H1.4_EB_stdata.tsv")
    Counts_H1.5 <- ImportOriginalSTCountData("H1_5", "./Patient 1/1k_arrays/H1_5/180903_L11_CN63_C1_H1.5_EB_stdata.tsv")

# Importing Patient 1, 1k array spot barcode files

    H2_1_SelectedBarcodes <- ImportHistologicalOriginalSTSelections("H2_1", "./Patient 1/1k_arrays/H2_1/spot_data-selection-180903_L11_CN63_D1_P_H2.1_CY3_EB_aligned.tsv")
    H2_2_SelectedBarcodes <- ImportHistologicalOriginalSTSelections("H2_2", "./Patient 1/1k_arrays/H2_2/spot_data-selection-180614_L11_53_C1_H2.2_Cy3_EB_aligned.tsv")
    H2_3_SelectedBarcodes <- ImportHistologicalOriginalSTSelections("H2_3", "./Patient 1/1k_arrays/H2_3/spot_data-selection-180614_L11_53_D2_H2.3_Cy3_aligned.tsv")
    H2_4_SelectedBarcodes <- ImportHistologicalOriginalSTSelections("H2_4", "./Patient 1/1k_arrays/H2_4/spot_data-selection-180614_L11_53_E1_H2.4_Cy3_EB_aligned.tsv")
    H2_5_SelectedBarcodes <- ImportHistologicalOriginalSTSelections("H2_5", "./Patient 1/1k_arrays/H2_5/spot_data-selection-180614_L11_55_C1_H2.5_Cy3_aligned.tsv")

    V2_1_SelectedBarcodes <- ImportHistologicalOriginalSTSelections("V2_1", "./Patient 1/1k_arrays/V2_1/spot_data-selection-180614_L11_55_D1_V2.1_Cy3_EB_aligned.tsv")
    V2_2_SelectedBarcodes <- ImportHistologicalOriginalSTSelections("V2_2", "./Patient 1/1k_arrays/V2_2/spot_data-selection-180614_L11_55_E1_V2.2.tsv")
    V2_3_SelectedBarcodes <- ImportHistologicalOriginalSTSelections("V2_3", "./Patient 1/1k_arrays/V2_3/spot_data-selection-180813_L11_CN60_C2_P_V2.3_CY3_EB_aligned.tsv")
    V2_4_SelectedBarcodes <- ImportHistologicalOriginalSTSelections("V2_4", "./Patient 1/1k_arrays/V2_4/spot_data-selection-180813_L11_CN60_D2_P_V2.4_CY3_EB_aligned.tsv")
    V2_5_SelectedBarcodes <- ImportHistologicalOriginalSTSelections("V2_5", "./Patient 1/1k_arrays/V2_5/spot_data-selection-180822_L11_CN57_C1_P_V2.5_CY3_EB_aligned.tsv")
    V2_6_SelectedBarcodes <- ImportHistologicalOriginalSTSelections("V2_6", "./Patient 1/1k_arrays/V2_6/spot_data-selection-180813_L11_CN15_C2_P_V2.6_CY3_EB_aligned.tsv")

    V1_1_SelectedBarcodes <- ImportHistologicalOriginalSTSelections("V1_1", "./Patient 1/1k_arrays/V1_1/spot_data-selection-180813_L11_CN15_D1_P_V1.1_CY3_EB_aligned.tsv")
    V1_2_SelectedBarcodes <- ImportHistologicalOriginalSTSelections("V1_2", "./Patient 1/1k_arrays/V1_2/spot_data-selection-180813_L11_CN15_E1_P_V1.2_CY3_EB_aligned.tsv")
    V1_3_SelectedBarcodes <- ImportHistologicalOriginalSTSelections("V1_3", "./Patient 1/1k_arrays/V1_3/spot_data-selection-180822_L11_CN56_C1_P_V1.3_CY3_EB_aligned.tsv")
    V1_4_SelectedBarcodes <- ImportHistologicalOriginalSTSelections("V1_4", "./Patient 1/1k_arrays/V1_4/spot_data-selection-180822_L11_CN56_D1_P_V1.4_CY3_EB_aligned.tsv")
    V1_5_SelectedBarcodes <- ImportHistologicalOriginalSTSelections("V1_5", "./Patient 1/1k_arrays/V1_5/spot_data-selection-180822_L11_CN56_E2_P_V1.5_CY3_EB_aligned.tsv")

    H1_1_SelectedBarcodes <- ImportHistologicalOriginalSTSelections("H1_1", "./Patient 1/1k_arrays/H1_1/spot_data-selection-180822_L11_CN57_E1_P_H1.1_CY3_EB_aliged.tsv")
    H1_2_SelectedBarcodes <- ImportHistologicalOriginalSTSelections("H1_2", "./Patient 1/1k_arrays/H1_2/spot_data-selection-180903_L11_CN51_C2_P_H1.2_CY3_EB_aligned.tsv")
    H1_3_SelectedBarcodes <- ImportHistologicalOriginalSTSelections("H1_3", "./Patient 1/1k_arrays/H1_3/spot_data-selection-180903_L11_CN61_D2_P_H1.3_CY3_EB_aligned.tsv")
    H1_4_SelectedBarcodes <- ImportHistologicalOriginalSTSelections("H1_4", "./Patient 1/1k_arrays/H1_4/spot_data-selection-180903_L11_CN51_E1_P_H1.4_CY3_EB_aligned.tsv")
    H1_5_SelectedBarcodes <- ImportHistologicalOriginalSTSelections("H1_5", "./Patient 1/1k_arrays/H1_5/spot_data-selection-180903_L11_CN63_C1_P_H1.5_CY3_EB_aligned.tsv")

    save.image(file = "Counts_my_work_space.RData")

# Creating a barcode dataframe of all barcodes

    Barcodes_H2_1 <- H2_1_SelectedBarcodes %>% 
                      select(Barcode) %>%
                      mutate(Histology = "H2_1")
    Barcodes_H2_2 <- H2_2_SelectedBarcodes %>% 
                      select(Barcode) %>%
                      mutate(Histology = "H2_2")
    Barcodes_H2_3 <- H2_3_SelectedBarcodes %>% 
                      select(Barcode) %>%
                      mutate(Histology = "H2_3")
    Barcodes_H2_4 <- H2_4_SelectedBarcodes %>% 
                      select(Barcode) %>%
                      mutate(Histology = "H2_4")
    Barcodes_H2_5 <- H2_5_SelectedBarcodes %>% 
                      select(Barcode) %>%
                      mutate(Histology = "H2_5")

    Barcodes_V2_1 <- V2_1_SelectedBarcodes %>% 
                      select(Barcode) %>%
                      mutate(Histology = "V2_1")
    Barcodes_V2_2 <- V2_2_SelectedBarcodes %>% 
                      select(Barcode) %>%
                      mutate(Histology = "V2_2")
    Barcodes_V2_3 <- V2_3_SelectedBarcodes %>% 
                      select(Barcode) %>%
                      mutate(Histology = "V2_3")
    Barcodes_V2_4 <- V2_4_SelectedBarcodes %>% 
                      select(Barcode) %>%
                      mutate(Histology = "V2_4")
    Barcodes_V2_5 <- V2_5_SelectedBarcodes %>% 
                      select(Barcode) %>%
                      mutate(Histology = "V2_5")
    Barcodes_V2_6 <- V2_6_SelectedBarcodes %>% 
                      select(Barcode) %>%
                      mutate(Histology = "V2_6")

    Barcodes_V1_1 <- V1_1_SelectedBarcodes %>% 
                      select(Barcode) %>%
                      mutate(Histology = "V1_1")
    Barcodes_V1_2 <- V1_2_SelectedBarcodes %>% 
                      select(Barcode) %>%
                      mutate(Histology = "V1_2")
    Barcodes_V1_3 <- V1_3_SelectedBarcodes %>% 
                      select(Barcode) %>%
                      mutate(Histology = "V1_3")
    Barcodes_V1_4 <- V1_4_SelectedBarcodes %>% 
                      select(Barcode) %>%
                      mutate(Histology = "V1_4")
    Barcodes_V1_5 <- V1_5_SelectedBarcodes %>% 
                      select(Barcode) %>%
                      mutate(Histology = "V1_5")

    Barcodes_H1_1 <- H1_1_SelectedBarcodes %>% 
                      select(Barcode) %>%
                      mutate(Histology = "H1_1")
    Barcodes_H1_2 <- H1_2_SelectedBarcodes %>% 
                      select(Barcode) %>%
                      mutate(Histology = "H1_2")
    Barcodes_H1_3 <- H1_3_SelectedBarcodes %>% 
                      select(Barcode) %>%
                      mutate(Histology = "H1_3")
    Barcodes_H1_4 <- H1_4_SelectedBarcodes %>% 
                      select(Barcode) %>%
                      mutate(Histology = "H1_4")
    Barcodes_H1_5 <- H1_5_SelectedBarcodes %>% 
                      select(Barcode) %>%
                      mutate(Histology = "H1_5")

    JoinedBarcodes <- rbind(H2_1_SelectedBarcodes,
                            H2_2_SelectedBarcodes,
                            H2_3_SelectedBarcodes,
                            H2_4_SelectedBarcodes,
                            H2_5_SelectedBarcodes,
                            V2_1_SelectedBarcodes,
                            V2_2_SelectedBarcodes,
                            V2_3_SelectedBarcodes,
                            V2_4_SelectedBarcodes,
                            V2_5_SelectedBarcodes,
                            V2_6_SelectedBarcodes,
                            V1_1_SelectedBarcodes,
                            V1_2_SelectedBarcodes,
                            V1_3_SelectedBarcodes,
                            V1_4_SelectedBarcodes,
                            V1_5_SelectedBarcodes,
                            H1_1_SelectedBarcodes,
                            H1_2_SelectedBarcodes,
                            H1_3_SelectedBarcodes,
                            H1_4_SelectedBarcodes,
                            H1_5_SelectedBarcodes)

    rm(H2_1_SelectedBarcodes,
                            H2_2_SelectedBarcodes,
                            H2_3_SelectedBarcodes,
                            H2_4_SelectedBarcodes,
                            H2_5_SelectedBarcodes,
                            V2_1_SelectedBarcodes,
                            V2_2_SelectedBarcodes,
                            V2_3_SelectedBarcodes,
                            V2_4_SelectedBarcodes,
                            V2_5_SelectedBarcodes,
                            V2_6_SelectedBarcodes,
                            V1_1_SelectedBarcodes,
                            V1_2_SelectedBarcodes,
                            V1_3_SelectedBarcodes,
                            V1_4_SelectedBarcodes,
                            V1_5_SelectedBarcodes,
                            H1_1_SelectedBarcodes,
                            H1_2_SelectedBarcodes,
                            H1_3_SelectedBarcodes,
                            H1_4_SelectedBarcodes,
                            H1_5_SelectedBarcodes)


    saveRDS(JoinedBarcodes, file = "Fig1d_Barcodes.rds")

# Selecting and QCing count files

    JoinedCounts_H2.1 <- OriginalST_MergingCountAndAnnotationData(Barcodes_H2_1, Counts_H2.1)
    JoinedCounts_H2.2 <- OriginalST_MergingCountAndAnnotationData(Barcodes_H2_2, Counts_H2.2)
    JoinedCounts_H2.3 <- OriginalST_MergingCountAndAnnotationData(Barcodes_H2_3, Counts_H2.3)
    JoinedCounts_H2.4 <- OriginalST_MergingCountAndAnnotationData(Barcodes_H2_4, Counts_H2.4)
    JoinedCounts_H2.5 <- OriginalST_MergingCountAndAnnotationData(Barcodes_H2_5, Counts_H2.5)

    JoinedCounts_V2.1 <- OriginalST_MergingCountAndAnnotationData(Barcodes_V2_1, Counts_V2.1)
    JoinedCounts_V2.2 <- OriginalST_MergingCountAndAnnotationData(Barcodes_V2_2, Counts_V2.2)
    JoinedCounts_V2.3 <- OriginalST_MergingCountAndAnnotationData(Barcodes_V2_3, Counts_V2.3)
    JoinedCounts_V2.4 <- OriginalST_MergingCountAndAnnotationData(Barcodes_V2_4, Counts_V2.4)
    JoinedCounts_V2.5 <- OriginalST_MergingCountAndAnnotationData(Barcodes_V2_5, Counts_V2.5)
    JoinedCounts_V2.6 <- OriginalST_MergingCountAndAnnotationData(Barcodes_V2_6, Counts_V2.6)

    JoinedCounts_V1.1 <- OriginalST_MergingCountAndAnnotationData(Barcodes_V1_1, Counts_V1.1)
    JoinedCounts_V1.2 <- OriginalST_MergingCountAndAnnotationData(Barcodes_V1_2, Counts_V1.2)
    JoinedCounts_V1.3 <- OriginalST_MergingCountAndAnnotationData(Barcodes_V1_3, Counts_V1.3)
    JoinedCounts_V1.4 <- OriginalST_MergingCountAndAnnotationData(Barcodes_V1_4, Counts_V1.4)
    JoinedCounts_V1.5 <- OriginalST_MergingCountAndAnnotationData(Barcodes_V1_5, Counts_V1.5)

    JoinedCounts_H1.1 <- OriginalST_MergingCountAndAnnotationData(Barcodes_H1_1, Counts_H1.1)
    JoinedCounts_H1.2 <- OriginalST_MergingCountAndAnnotationData(Barcodes_H1_2, Counts_H1.2)
    JoinedCounts_H1.3 <- OriginalST_MergingCountAndAnnotationData(Barcodes_H1_3, Counts_H1.3)
    JoinedCounts_H1.4 <- OriginalST_MergingCountAndAnnotationData(Barcodes_H1_4, Counts_H1.4)
    JoinedCounts_H1.5 <- OriginalST_MergingCountAndAnnotationData(Barcodes_H1_5, Counts_H1.5)


    rm(Counts_H2.1,
       Counts_H2.2, 
       Counts_H2.3, 
       Counts_H2.4, 
       Counts_H2.5, 
       Counts_V2.1, 
       Counts_V2.2,
       Counts_V2.3,
       Counts_V2.4,
       Counts_V2.5,
       Counts_V2.6, 
       Counts_V1.1, 
       Counts_V1.2, 
       Counts_V1.3,
       Counts_V1.4,
       Counts_V1.5, 
       Counts_H1.1, 
       Counts_H1.2, 
       Counts_H1.3,
       Counts_H1.4, 
       Counts_H1.5
       )
    rm(Barcodes_H2_1,       Barcodes_H2_2,
                            Barcodes_H2_3,
                            Barcodes_H2_4,
                            Barcodes_H2_5,
                            Barcodes_V2_1,
                            Barcodes_V2_2,
                            Barcodes_V2_3,
                            Barcodes_V2_4,
                            Barcodes_V2_5,
                            Barcodes_V2_6,
                            Barcodes_V1_1,
                            Barcodes_V1_2,
                            Barcodes_V1_3,
                            Barcodes_V1_4,
                            Barcodes_V1_5,
                            Barcodes_H1_1,
                            Barcodes_H1_2,
                            Barcodes_H1_3,
                            Barcodes_H1_4,
                            Barcodes_H1_5)

    save.image(file = "Fig1d_joined_counts_my_work_space.RData")

# Joining all of the QC’d counts

    Counts_joined <- JoinedCounts_H2.1 %>% full_join(JoinedCounts_H2.2, by = "Genes")
    rm(JoinedCounts_H2.1, JoinedCounts_H2.2)
    Counts_joined <- Counts_joined %>% full_join(JoinedCounts_H2.3, by = "Genes")
    rm(JoinedCounts_H2.3)
    Counts_joined <- Counts_joined %>% full_join(JoinedCounts_H2.4, by = "Genes")
    rm(JoinedCounts_H2.4)
    Counts_joined <- Counts_joined %>% full_join(JoinedCounts_H2.5, by = "Genes")
    rm(JoinedCounts_H2.5)

    Counts_joined <- Counts_joined %>% full_join(JoinedCounts_V2.1, by = "Genes")
    rm(JoinedCounts_V2.1)
    Counts_joined <- Counts_joined %>% full_join(JoinedCounts_V2.2, by = "Genes")
    rm(JoinedCounts_V2.2)
    Counts_joined <- Counts_joined %>% full_join(JoinedCounts_V2.3, by = "Genes")
    rm(JoinedCounts_V2.3)
    Counts_joined <- Counts_joined %>% full_join(JoinedCounts_V2.4, by = "Genes")
    rm(JoinedCounts_V2.4)
    Counts_joined <- Counts_joined %>% full_join(JoinedCounts_V2.5, by = "Genes")
    rm(JoinedCounts_V2.5)
    Counts_joined <- Counts_joined %>% full_join(JoinedCounts_V2.6, by = "Genes")
    rm(JoinedCounts_V2.6)

    Counts_joined <- Counts_joined %>% full_join(JoinedCounts_V1.1, by = "Genes")
    rm(JoinedCounts_V1.1)
    Counts_joined <- Counts_joined %>% full_join(JoinedCounts_V1.2, by = "Genes")
    rm(JoinedCounts_V1.2)
    Counts_joined <- Counts_joined %>% full_join(JoinedCounts_V1.4, by = "Genes")
    rm(JoinedCounts_V1.4)
    Counts_joined <- Counts_joined %>% full_join(JoinedCounts_V1.3, by = "Genes")
    rm(JoinedCounts_V1.3)
    Counts_joined <- Counts_joined %>% full_join(JoinedCounts_V1.5, by = "Genes")
    rm(JoinedCounts_V1.5)

    Counts_joined <- Counts_joined %>% full_join(JoinedCounts_H1.1, by = "Genes")
    rm(JoinedCounts_H1.1)
    Counts_joined <- Counts_joined %>% full_join(JoinedCounts_H1.2, by = "Genes")
    rm(JoinedCounts_H1.2)
    Counts_joined <- Counts_joined %>% full_join(JoinedCounts_H1.3, by = "Genes")
    rm(JoinedCounts_H1.3)
    Counts_joined <- Counts_joined %>% full_join(JoinedCounts_H1.4, by = "Genes")
    rm(JoinedCounts_H1.4)
    Counts_joined <- Counts_joined %>% full_join(JoinedCounts_H1.5, by = "Genes")
    rm(JoinedCounts_H1.5)


    saveRDS(Counts_joined, file = "Fig1d_Counts_joined.rds")

# For the joined counts, replacing NAs

    Counts_joined <- Counts_joined %>% replace(., is.na(.), 0)

    saveRDS(Counts_joined, file = "Fig1d_noNA_Counts_joined.rds")

    Fig1d_Counts_joined <- Counts_joined

# Mapping the gene names (from count files) to chromosome coordinates to prepare for the gene order file

    #library(tidyverse)
    #library(data.table)
    GeneToENSMBL <- fread('https://data.broadinstitute.org/Trinity/CTAT/cnv/gencode_v19_gen_pos.complete.txt')
    GeneToENSMBL <- GeneToENSMBL %>% separate(V1, c("left","ENSMBLID"), sep = "\\|")

    names(GeneToENSMBL)[1] <- "Genes"
    names(GeneToENSMBL)[3] <- "chr"
    names(GeneToENSMBL)[4] <- "start"
    names(GeneToENSMBL)[5] <- "stop"

    #removing "."
    GeneToENSMBL <- GeneToENSMBL %>% separate(ENSMBLID, c("ENSMBLID", NA))

    Counts_joined <- Fig1d_Counts_joined
    Counts_joined <- Counts_joined %>%
                        separate(Genes, c("Genes", NA))

    Counts_joined <- Counts_joined %>% select(Genes)

    GenesForMapping <- GeneToENSMBL %>% select(ENSMBLID, chr, start, stop)
    GenesInSample <- Counts_joined %>% select(Genes)
    GenesInSamplevsOrdering <- inner_join(GenesInSample, GenesForMapping, by = c("Genes" = "ENSMBLID"))
      dedup_GenesInSamplevsOrdering <- GenesInSamplevsOrdering[!duplicated(GenesInSamplevsOrdering$Genes), ]
      dedup_GenesInSamplevsOrdering$chromorder <- gsub("chr","",dedup_GenesInSamplevsOrdering$chr)
      dedup_GenesInSamplevsOrdering$chromorder <- as.numeric(ifelse(dedup_GenesInSamplevsOrdering$chromorder == "X", 23,
                                                             ifelse(dedup_GenesInSamplevsOrdering$chromorder == "Y", 24,      dedup_GenesInSamplevsOrdering$chromorder)))
      dedup_GenesInSamplevsOrdering <- dedup_GenesInSamplevsOrdering[order(dedup_GenesInSamplevsOrdering$chromorder),]
      dedup_GenesInSamplevsOrdering <- dedup_GenesInSamplevsOrdering[,1:4]  

      
    MappingFileForInferCNV <- dedup_GenesInSamplevsOrdering

    saveRDS(MappingFileForInferCNV, file = "Fig1d_MappingFileForInferCNV.rds")  

# Selecting only counts from genes that are mapped in the gene order file, and then outputting all 3 requisite files for inferCNV::run

    CountmappedGenes <- select(MappingFileForInferCNV, Genes)

    Counts_joined <- Fig1d_Counts_joined
    Counts_joined <- Counts_joined %>%
                        separate(Genes, c("Genes", NA))

    Mapped_Counts_joined <- left_join(CountmappedGenes, Counts_joined)

    rm(Counts_joined)

    Mapped_Counts_joinedSliced <- Mapped_Counts_joined %>% slice(1L)
    Mapped_Counts_joinedSliced <- as.data.frame(t(Mapped_Counts_joinedSliced[, colnames(Mapped_Counts_joinedSliced)[c(1:length(Mapped_Counts_joinedSliced))]]))
    Mapped_Counts_joinedSliced <- Mapped_Counts_joinedSliced %>% rownames_to_column()
    Mapped_Counts_joinedSliced <- as.data.frame(Mapped_Counts_joinedSliced[2:(dim(Mapped_Counts_joinedSliced)[1]), 1])
    names(Mapped_Counts_joinedSliced)[1] <- "Barcode"

    JoinedBarcodes <- readRDS("Fig1d_Barcodes.rds")
    CorrectedBarcodes <- Mapped_Counts_joinedSliced
    CorrectedBarcodes$Histology <- paste0(substr(CorrectedBarcodes$Barcode, start = 1, stop = 4)
    )

    #Write GenesInSamplevsOrdering
    write.table(Mapped_Counts_joined, 
                "Fig1d_STOrganscale_Selected_Mapped_Counts.tsv",
                row.names = FALSE,
                sep = "\t")

    write.table(CorrectedBarcodes, 
                "Fig1d_STOrganscale_Selected_Mapped_Annotations.tsv", 
                quote = FALSE, 
                col.names = FALSE, 
                row.names = FALSE, 
                sep = "\t")

    write.table(MappingFileForInferCNV, 
                "Fig1d_STOrganscale_GeneOrderFile.tsv", 
                quote = FALSE, 
                col.names = FALSE, 
                row.names = FALSE, 
                sep = "\t")

# Creating the inferCNV object (prior to run)

    STOrganscale <- infercnv::CreateInfercnvObject(raw_counts_matrix="./Fig1d_STOrganscale_Selected_Mapped_Counts.tsv", 
                                                   gene_order_file="./Fig1d_STOrganscale_GeneOrderFile.tsv",
                                                   annotations_file="./Fig1d_STOrganscale_Selected_Mapped_Annotations.tsv",
                                                   delim="\t",
                                                   ref_group_names=NULL,
                                                                   chr_exclude = c("chrM"))

# Unsupervised Run - (Typically ran on cluster)

    STOrganscale = infercnv::run(STOrganscale,
                                 cutoff=0.1,
                                 out_dir="./Fig1_Step1/Outputs",
                                 cluster_by_groups=TRUE,
                                 num_threads = 20, 
                                 denoise=TRUE,
                                 HMM=TRUE,
                                 analysis_mode = "cells",
                                 HMM_report_by = "cell")