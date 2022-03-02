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

    library(SpatialInferCNV)

    ## Registered S3 method overwritten by 'spatstat.geom':
    ##   method     from
    ##   print.boxx cli

    ## Warning: replacing previous import 'phylogram::as.phylo' by 'ape::as.phylo' when
    ## loading 'SpatialInferCNV'

# Importing Data for Benigns

    #Import SCC, Patient 6, scRNAseq benigns that we subset out in step 1
    load("./Figure4c - Step 1 - Downloading scRNAseq Data/SCC_P6_Benigns.RData")

    head(SCC_P6_Benigns)

    SCC_P6_BenignReferences_Barcodes <- read.csv("./Figure4c_SCCP6_BenignReferenceSet.csv")
    names(SCC_P6_BenignReferences_Barcodes)[1] <- "Barcodes"
    SCC_P6_BenignReferences_Barcodes$Histology <- "PurestBenign_SCCPatient6"
    SCC_P6_Benigns <- SCC_P6_Benigns %>% rownames_to_column()
    names(SCC_P6_Benigns)[1] <- "Barcodes"

    SCC_P6_BenignReferences_Counts <- left_join(SCC_P6_BenignReferences_Barcodes, SCC_P6_Benigns, by = c("Barcodes" = "Barcodes"))
    rm(SCC_P6_Benigns)
    SCC_P6_BenignReferences_Counts <- SCC_P6_BenignReferences_Counts %>% select(-Histology)
    SCC_P6_BenignReferences_Counts <- SCC_P6_BenignReferences_Counts %>% column_to_rownames(var = "Barcodes")
    SCC_P6_BenignReferences_Counts <- as.data.frame(t(SCC_P6_BenignReferences_Counts))
    SCC_P6_BenignReferences_Counts <- SCC_P6_BenignReferences_Counts %>% rownames_to_column()
    names(SCC_P6_BenignReferences_Counts)[1] <- "Genes"

    saveRDS(SCC_P6_BenignReferences_Counts, file = "SCC_P6_BenignReferences_Counts.rds")
    saveRDS(SCC_P6_BenignReferences_Barcodes, file = "SCC_P6_BenignReferences_Barcodes.rds")

# Importing Data for Visium Data

    t28 <- readRDS("./Mendeley/SCC Patient/t28.Rds")

    SCC_P6_Visium_Counts <- as.data.frame(t28@assays$Spatial@counts)
    rm(t28)

    head(SCC_P6_Visium_Counts)

    SCC_P6_Visium_Counts <- as.data.frame(t(SCC_P6_Visium_Counts))
    SCC_P6_Visium_Counts <- rownames_to_column(SCC_P6_Visium_Counts)
    SCC_P6_Visium_Counts$section <- str_sub(SCC_P6_Visium_Counts$rowname, start= -1)
    table(SCC_P6_Visium_Counts$section)

    SCC_P6_Visium_Counts$barcode <- str_sub(SCC_P6_Visium_Counts$rowname, start = 1L, end = -3)

    SCC_P6_Visium_Counts <- SCC_P6_Visium_Counts %>% filter(section == 1)

    SCC_P6_Visium_Annotations <- SCC_P6_Visium_Counts %>% select(barcode, section)
    SCC_P6_Visium_Annotations$section <- "SCC_P6_Visium"
    names(SCC_P6_Visium_Annotations)[1] <- "Barcodes"
    names(SCC_P6_Visium_Annotations)[2] <- "Histology"

    saveRDS(SCC_P6_Visium_Annotations, file = "SCC_P6_Visium_Annotations.rds")

    SCC_P6_Visium_Counts <- column_to_rownames(SCC_P6_Visium_Counts, var = "barcode")
    SCC_P6_Visium_Counts <- SCC_P6_Visium_Counts %>% select(-rowname, -section)

    SCC_P6_Visium_Counts$Total <- rowSums(SCC_P6_Visium_Counts)
    SCC_P6_Visium_Counts <- SCC_P6_Visium_Counts %>% filter(Total >= 500)
    SCC_P6_Visium_Counts <- select(SCC_P6_Visium_Counts, -Total)
    SCC_P6_Visium_Counts <- as.data.frame(t(SCC_P6_Visium_Counts))
    SCC_P6_Visium_Counts <- SCC_P6_Visium_Counts[,colSums(is.na(SCC_P6_Visium_Counts))<nrow(SCC_P6_Visium_Counts)]
    SCC_P6_Visium_Counts <- tibble::rownames_to_column(SCC_P6_Visium_Counts, "Genes")

    saveRDS(SCC_P6_Visium_Counts, file = "SCC_P6_Visium_Counts.rds")

# Importing Data Visium and Benign Data

    SCC_P6_Visium_Annotations <- readRDS("./SCC_P6_Visium_Annotations.rds")
    SCC_P6_BenignReferences_Barcodes <- readRDS("./SCC_P6_BenignReferences_Barcodes.rds")

    Joined_Barcodes <- rbind(SCC_P6_Visium_Annotations, SCC_P6_BenignReferences_Barcodes)
    saveRDS(Joined_Barcodes, file = "SCC_P6_BenignRef_and_Visium_Annotations.rds")

# Importing Data Visium and Benign Data

    SCC_P6_BenignReferences_Counts <- readRDS("./SCC_P6_BenignReferences_Counts.rds")
    SCC_P6_Visium_Counts <- readRDS("./SCC_P6_Visium_Counts.rds")

    head(SCC_P6_BenignReferences_Counts)
    head(SCC_P6_Visium_Counts)

    SCC_P6_BenignRef_and_Visium_Counts <- SCC_P6_BenignReferences_Counts %>% full_join(SCC_P6_Visium_Counts, by = "Genes")
    SCC_P6_BenignRef_and_Visium_Counts <- SCC_P6_BenignRef_and_Visium_Counts %>% replace(., is.na(.), 0)

    saveRDS(SCC_P6_BenignRef_and_Visium_Counts, file = "SCC_P6_BenignRef_and_Visium_Counts.rds")

# Creating GeneToENSMBL dataframe

    GeneToENSMBL <- read.csv("./Mendeley/ProcessedFilesForFigures/Figure4/GeneToENSMBL.csv")

    #library(tidyverse)
    #library(data.table)
    #GeneToENSMBL <- fread('https://data.broadinstitute.org/Trinity/CTAT/cnv/gencode_v19_gen_pos.complete.txt')
    #GeneToENSMBL <- mydat %>% separate(V1, c("left","ENSMBLID"), sep = "\\|")

    #names(GeneToENSMBL)[1] <- "Genes"
    #names(GeneToENSMBL)[3] <- "chr"
    #names(GeneToENSMBL)[4] <- "start"
    #names(GeneToENSMBL)[5] <- "stop"

    #write.csv(GeneToENSMBL, "GeneToENSMBL.csv", row.names = FALSE)

# Mapping Gene Names to counts/barcodes, and then outputting the requisite files for infercnv::run, part 1

    #removing "."
    Counts_joined <- SCC_P6_BenignRef_and_Visium_Counts
    Counts_joined <- Counts_joined %>%
                        separate(Genes, c("Genes", NA))

    Counts_joined <- Counts_joined %>% select(Genes)

    GenesForMapping <- GeneToENSMBL %>% select(Genes, chr, start, stop)
    GenesInSample <- Counts_joined %>% select(Genes)
    GenesInSamplevsOrdering <- inner_join(GenesInSample, GenesForMapping, by = c("Genes" = "Genes"))
      dedup_GenesInSamplevsOrdering <- GenesInSamplevsOrdering[!duplicated(GenesInSamplevsOrdering$Genes), ]
      dedup_GenesInSamplevsOrdering$chromorder <- gsub("chr","",dedup_GenesInSamplevsOrdering$chr)
      dedup_GenesInSamplevsOrdering$chromorder <- as.numeric(ifelse(dedup_GenesInSamplevsOrdering$chromorder == "X", 23,
                                                             ifelse(dedup_GenesInSamplevsOrdering$chromorder == "Y", 24,      dedup_GenesInSamplevsOrdering$chromorder)))
      dedup_GenesInSamplevsOrdering <- dedup_GenesInSamplevsOrdering[order(dedup_GenesInSamplevsOrdering$chromorder),]
      dedup_GenesInSamplevsOrdering <- dedup_GenesInSamplevsOrdering[,1:4]  
      
    MappingFileForInferCNV <- dedup_GenesInSamplevsOrdering

    saveRDS(MappingFileForInferCNV, file = "MappingFileForSCC_P6_Visium_and_Bg.rds")  

# Outputting the requisite files for infercnv::run, part 2

    MappingFileForInferCNV <- readRDS("MappingFileForSCC_P6_Visium_and_Bg.rds")
    SCC_P6_BenignRef_and_Visium_Counts <- readRDS("SCC_P6_BenignRef_and_Visium_Counts.rds")

    CountmappedGenes <- select(MappingFileForInferCNV, Genes)

    Counts_joined <- SCC_P6_BenignRef_and_Visium_Counts
    Counts_joined <- Counts_joined %>%
                        separate(Genes, c("Genes", NA))

    Mapped_Counts_joined <- left_join(CountmappedGenes, Counts_joined)
    Mapped_Counts_joined <- Mapped_Counts_joined[!duplicated(Mapped_Counts_joined$Genes), ]
    Mapped_Counts_joinedSliced <- Mapped_Counts_joined %>% slice(1L)
    Mapped_Counts_joinedSliced <- as.data.frame(t(Mapped_Counts_joinedSliced[, colnames(Mapped_Counts_joinedSliced)[c(1:length(Mapped_Counts_joinedSliced))]]))
    Mapped_Counts_joinedSliced <- Mapped_Counts_joinedSliced %>% rownames_to_column()
    Mapped_Counts_joinedSliced <- as.data.frame(Mapped_Counts_joinedSliced[2:(dim(Mapped_Counts_joinedSliced)[1]), 1])
    names(Mapped_Counts_joinedSliced)[1] <- "Barcode"

    Mapped_Counts_joinedSliced$Histology <- ifelse(paste0(substr(Mapped_Counts_joinedSliced$Barcode, start = 1, stop = 4)) == "P6_N", "PurestBenign_SCCPatient6", "Visium")

    #P6_Normal
    #Write GenesInSamplevsOrdering
    write.table(Mapped_Counts_joined, 
                "SCC_P6_BenignRef_and_Visium_Mapped_Counts.tsv",
                row.names = FALSE,
                sep = "\t")

    write.table(MappingFileForInferCNV, 
                "SCC_P6_BenignRef_and_Visium_GeneOrderFile.tsv", 
                quote = FALSE, 
                col.names = FALSE, 
                row.names = FALSE, 
                sep = "\t")


    write.table(Mapped_Counts_joinedSliced, 
                "SCC_P6_BenignRef_and_Visium_Mapped_Annotations.tsv", 
                quote = FALSE, 
                col.names = FALSE, 
                row.names = FALSE, 
                sep = "\t")

# Creating the inferCNV object (prior to run)

    Visium_P6_Bg_infCNV <- infercnv::CreateInfercnvObject(raw_counts_matrix="./SCC_P6_BenignRef_and_Visium_Mapped_Counts.tsv", 
                                                   gene_order_file="./SCC_P6_BenignRef_and_Visium_GeneOrderFile.tsv",
                                                   annotations_file="./SCC_P6_BenignRef_and_Visium_Mapped_Annotations.tsv",
                                                   delim="\t",
                                                   ref_group_names="PurestBenign_SCCPatient6",
                                                   chr_exclude = c("chrM"))

# Unsupervised Run - (Typically ran on cluster)

    Visium_P6_Bg_infCNV = infercnv::run(Visium_P6_Bg_infCNV,
                                                  cutoff=0.1,
                                                out_dir="./Outputs", 
                                                  num_threads = 10,
                                                  cluster_by_groups=FALSE, 
                                                  denoise=TRUE,
                                                  HMM=FALSE)
