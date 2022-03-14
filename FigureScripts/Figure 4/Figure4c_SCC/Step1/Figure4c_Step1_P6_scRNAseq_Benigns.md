# Figure 4c - Step 1 - Selectiong of Benign references (from paired scRNAseq data)

In order to run the SCC Visium data with siCNV, we need a reference set.
We identifed a set of paired scRNA sequencing data, from benign skin
cells (from the same exact patient), provided by the authors as listed
below.

# Setup

Initializing packages.

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

# Creating a working directory

    dir.create("Figure4c_output")
    setwd("Figure4c_output")

# Downloading and formatting data, part 1

Warning, this step will take 10-60 min, even with a decent internet
connection.

    counturl <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE144236&format=file&file=GSE144236%5FcSCC%5Fcounts%2Etxt%2Egz"
    tmp <- tempfile()
    download.file(counturl,tmp)

    #Warning, this next step will take 10-60 minutes
    merge10pts_counts <- read.delim(gzfile(tmp))
    merge10pts_counts <- as.data.frame(t(merge10pts_counts))

    SCC_P6_Benigns <- merge10pts_counts %>%
                filter(Patient == 6) %>%
                filter(`Tissue: 0=Normal, 1=Tumor` == 0)

    SCC_P6_Benigns <- SCC_P6_Benigns %>%
          select(-Patient, -`Tissue: 0=Normal, 1=Tumor`)

    save(SCC_P6_Benigns, file = "SCC_P6_Benigns.RData")

# Downloading and formatting data, part 2

We then select only patient 6 data (corresponds to the specific patient
in our analyses).

    #Import SCC, Patient 6, scRNAseq benigns that we subset out above
    load("./SCC_P6_Benigns.RData")

    #Following code creates a barcode dataframe that we will need later
    P6_Benigns_forannotations <- SCC_P6_Benigns %>% rownames_to_column()

    Barcodes_P6 <- P6_Benigns_forannotations %>% 
                      select(rowname) %>%
                      mutate(Histology = "P6_Benigns")
    names(Barcodes_P6)[1] <- "Barcodes"

    #Next, we will prepare the gene order file required for infercnv:run
    SCC_P6_Benigns <- as.data.frame(t(SCC_P6_Benigns))
    SCC_P6_Benigns <- SCC_P6_Benigns %>% rownames_to_column()
    names(SCC_P6_Benigns)[1] <- "Genes"

# Creating GeneToENSMBL dataframe

The code below creates the GeneToENSMBL.csv file, but we have provided
this on our GitHub:

![](https://github.com/aerickso/SpatialInferCNV/blob/main/FigureScripts/Figure%204/Figure4c_SCC/GeneToENSMBL.csv).

    GeneToENSMBL <- read.csv("./GeneToENSMBL.csv")

    #library(tidyverse)
    #library(data.table)
    #GeneToENSMBL <- fread('https://data.broadinstitute.org/Trinity/CTAT/cnv/gencode_v19_gen_pos.complete.txt')
    #GeneToENSMBL <- mydat %>% separate(V1, c("left","ENSMBLID"), sep = "\\|")

    #names(GeneToENSMBL)[1] <- "Genes"
    #names(GeneToENSMBL)[3] <- "chr"
    #names(GeneToENSMBL)[4] <- "start"
    #names(GeneToENSMBL)[5] <- "stop"

    #write.csv(GeneToENSMBL, "GeneToENSMBL.csv", row.names = FALSE)

# Mapping Gene Names to counts/barcodes, and then outputting the requisite files for infercnv::run

We need to provide a gene ordering file to inferCNV, in the form of:
Gene Name / Chromosome Number / Start Loci / Stop Loci. As the files
provided by the authors are in “Gene Name”, and our chromosomal / loci
information are mapped to ENSMBLID’s, we need to map the Gene Names to
ENSMBLIDs.

    Counts_joined <- SCC_P6_Benigns
    Counts_joined <- Counts_joined %>%
                        separate(Genes, c("Genes", NA))

    Counts_joined <- Counts_joined %>% select(Genes)

    #Selecting Gene name, chromosome, start and stop locations
    GenesForMapping <- GeneToENSMBL %>% select(Genes, chr, start, stop)
    GenesInSample <- Counts_joined %>% select(Genes)

    #Next, reordering the entries from Chromsomes 1-22, followed by X and Y
    GenesInSamplevsOrdering <- inner_join(GenesInSample, GenesForMapping, by = c("Genes" = "Genes"))
      dedup_GenesInSamplevsOrdering <- GenesInSamplevsOrdering[!duplicated(GenesInSamplevsOrdering$Genes), ]
      dedup_GenesInSamplevsOrdering$chromorder <- gsub("chr","",dedup_GenesInSamplevsOrdering$chr)
      dedup_GenesInSamplevsOrdering$chromorder <- as.numeric(ifelse(dedup_GenesInSamplevsOrdering$chromorder == "X", 23,
                                                             ifelse(dedup_GenesInSamplevsOrdering$chromorder == "Y", 24,      dedup_GenesInSamplevsOrdering$chromorder)))
      dedup_GenesInSamplevsOrdering <- dedup_GenesInSamplevsOrdering[order(dedup_GenesInSamplevsOrdering$chromorder),]
      dedup_GenesInSamplevsOrdering <- dedup_GenesInSamplevsOrdering[,1:4]  

    MappingFileForInferCNV <- dedup_GenesInSamplevsOrdering

    #Selecting only genes that have location data
    CountmappedGenes <- select(MappingFileForInferCNV, Genes)
    Counts_joined <- SCC_P6_Benigns
    Counts_joined <- Counts_joined %>%
                        separate(Genes, c("Genes", NA))

    #Selecting only genes that have location and count data
    Mapped_Counts_joined <- left_join(CountmappedGenes, Counts_joined)
    #Removing duplicates
    Mapped_Counts_joined <- Mapped_Counts_joined[!duplicated(Mapped_Counts_joined$Genes), ]

# Outputting all files for inferCNV::run

    #Write GenesInSamplevsOrdering
    write.table(Mapped_Counts_joined, 
                "SCC_P6_Bg_Selected_Mapped_Counts.tsv",
                row.names = FALSE,
                sep = "\t")

    write.table(Barcodes_P6, 
                "SCC_P6_Bg_Selected_CorrectedBarcodes.tsv", 
                quote = FALSE, 
                col.names = FALSE, 
                row.names = FALSE, 
                sep = "\t")

    write.table(MappingFileForInferCNV, 
                "SCC_P6_Bg_MappingFileForInferCNV.tsv", 
                quote = FALSE, 
                col.names = FALSE, 
                row.names = FALSE, 
                sep = "\t")

# Creating the inferCNV object (prior to run)

Creating the object for infercnv::run.

    P6_Bg_infCNV <- infercnv::CreateInfercnvObject(raw_counts_matrix="./SCC_P6_Bg_Selected_Mapped_Counts.tsv", 
                                                   gene_order_file="./SCC_P6_Bg_MappingFileForInferCNV.tsv",
                                                   annotations_file="./SCC_P6_Bg_Selected_CorrectedBarcodes.tsv",
                                                   delim="\t",
                                                   ref_group_names=NULL,
                                                   chr_exclude = c("chrM"))

# Unsupervised Run - (Typically ran on cluster)

Running infercnv.

    P6_Bg_infCNV = infercnv::run(P6_Bg_infCNV,
                                                  cutoff=0.1,
                                                  out_dir="./Figure4c_Step1/Outputs", 
                                                  num_threads = 20,
                                                  cluster_by_groups=FALSE, 
                                                  denoise=TRUE,
                                                  HMM=FALSE)

InferCNV will output many files. We are primarily interested in the
final “infercnv.21\_denoised.png” file, as well as the text file
associated with the dendrogram associated with the hierarchical
clustering on the left hand side of the image
(infercnv.21\_denoised.observations\_dendrogram.txt).

![infercnv.21\_denoised.png](https://github.com/aerickso/SpatialInferCNV/blob/main/FigureScripts/Figure%204/Figure4c_SCC/Step2/infercnv.21_denoised.png)

# Importing dendrogram

Next, we want to import this dendrogram file, this was created just
above.

    library(ape)
    library(phylogram)
    SCC_P6_benigns_for_clustering <- read.dendrogram(file = "./Figure4c_Step1/Outputs/infercnv.21_denoised.observations_dendrogram.txt")

    SCC_P6_benigns_for_clustering_phylo <- as.phylo(SCC_P6_benigns_for_clustering)

# Visualizing dendrogram node numbers

Next, we want to visualize the numbers associated with the nodes of
interest (clones). We output a large image file that allows us to
manually inspect which nodes (cells) should be selected the purest
benign references. Here, we want the cells with the least signal
possible.

    my.subtrees = subtrees(SCC_P6_benigns_for_clustering_phylo)  # subtrees() to subset

    png("SCC_P6_benigns_for_clustering_phylo.png",width=10000,height=2500, res = 300)
    plot(SCC_P6_benigns_for_clustering_phylo,show.tip.label = FALSE)
    nodelabels(text=1:SCC_P6_benigns_for_clustering_phylo$Nnode,node=1:SCC_P6_benigns_for_clustering_phylo$Nnode+Ntip(SCC_P6_benigns_for_clustering_phylo))
    dev.off()

We provide the image output here:

![SCC\_P6\_benigns\_for\_clustering\_phylo.png](https://github.com/aerickso/SpatialInferCNV/tree/main/FigureScripts/Figure%204/Figure4c_SCC/Step1/SCC_P6_benigns_for_clustering_phylo.png)

# Clone selection

Next, view the output .png file, which provides a (albeit cluttered)
labeling of the dendrogram tree nodes. Manually select individual nodes
that correspond with a distinct subclonal grouping or signal, that will
be taken forward for re-clustering. This can be iteratively tweaked with
the next step + spatial visualization til optimal. We provide more
details
[here](https://github.com/aerickso/SpatialInferCNV/blob/main/FigureScripts/BenignRefs_ForFigs2and3/BenignRefs.md),
and provide the finalized selected SCC purest benign nodes here.

    #A - 4034
    #B - 3605   
    #B - 3360  
    #B - 2316
    #B - 724
    #C - 2

# Selecting purest benigns

Next, after identifying the numerical nodes that correspond to
dendrogram branches that correspond with a given set of molecular
signals, we then manually select these nodes in R, apply a label, then
join them all together for use in the next step.

    library(SpatialInferCNV)
    library(tidyverse)

    Node4034 <- SelectingSubTreeData(my.subtrees, 4034)
    Node2 <- SelectingSubTreeData(my.subtrees, 2)

    Merged <- rbind(Node4034, Node2)
    table(Merged$Node)

    Merged$Node <- ifelse(Merged$Node == "Node_4034", "PurestBenigns", "OtherBenigns")
    names(Merged)[2] <- "Histology"

    BenignRefs <- filter(Merged, Histology == "PurestBenigns") %>%
                                                select(Barcode, Histology)

    write.csv(BenignRefs, "Figure4c_SCCP6_BenignReferenceSet.csv", row.names = FALSE)
