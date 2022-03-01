# Code for generating inferCNV plot of pediatric brain tumor patient 1

# Set working directory, this directory should also include a folder containing all necessary files called “InferCNV\_pediatric\_patient\_1”

    #setwd("type_in_the_path_to_your_working_directory")

    # Load R packages
    library(STutility)

    ## Loading required package: Seurat

    ## Attaching SeuratObject

    ## Loading required package: ggplot2

    ## Registered S3 method overwritten by 'imager':
    ##   method      from
    ##   plot.imlist

    library(infercnv)

    ## Registered S3 method overwritten by 'ape':
    ##   method   from 
    ##   plot.mst spdep

    ## Registered S3 method overwritten by 'gplots':
    ##   method         from 
    ##   reorder.factor gdata

# Loading Data

    # Load infoTables: Sample from patient 3 included regions of stroma cells that will be excluded from the inferCNV analysis, which only contains spots from tumor regions.
    infoTable_pat_1_2 <- read.table("./inferCNV_pediatric_patient_1/infoTable_pat_1_2.csv", sep=";", header=T, stringsAsFactors = F)
    infoTable_pat_3 <- read.table("./inferCNV_pediatric_patient_1/infoTable_pat_3.csv", sep=";", header=T, stringsAsFactors = F)

    # Creat Seurat Objects
    se_pat_1_2 <- InputFromTable(infotable = infoTable_pat_1_2, 
                               min.gene.count = 100, 
                               min.gene.spots = 5,
                               min.spot.count = 500,
                               platform="Visium")


    se_pat_3 <- InputFromTable(infotable = infoTable_pat_3, 
                         min.gene.count = 100, 
                         min.gene.spots = 5,
                         min.spot.count = 500,
                         platform="Visium")

# Further Formatting

    # Add Pathology annotations to Meta.data in se_pat_3 object
    df <- read.csv(file = "./inferCNV_pediatric_patient_1/pathology_patient_3.csv")
    df$Barcode <- paste0(df$Barcode, "_1")
    rownames(df) <- df$Barcode
    se_pat_3$pathology <- df[rownames(se_pat_3[[]]), ]$Pathology

    # Check that pathology data was added to meta data of se_pat_3
    head(se_pat_3[[]])
    tail(se_pat_3[[]])
    table(se_pat_3$pathology)

    # Subsetting se_pat_3 to only contain spots with annotated tumor cells
    se_pat_3 <- SetIdent(se_pat_3, value = "pathology")
    se_pat_3 <- SubsetSTData(se_pat_3, idents = c("tumor cells"))

    # Check that only spots containing tumor cells are left
    table(se_pat_3$pathology)

    # Merge the se objects 
    se <- MergeSTData(se_pat_1_2, y = c(se_pat_3))

    # Check that the merge worked
    se
    head(se[[]])
    tail(se[[]])
    table(se$sample)

    # Set ident to sample
    se <- SetIdent(se, value = "sample")
    table(se$sample)

    # prepare a data.frame used as input for the inferCNV run
    se_sample <- as.data.frame(se$sample)
    head(se_sample)
    colnames(se_sample) <- c("sample")
    head(se_sample)
    se_sample <- cbind(Barcode = rownames(se_sample), se_sample)
    rownames(se_sample) <- NULL
    head(se_sample)
    tail(se_sample)

# Outputting Files

    # save the data.frame
    write.table(x = se_sample, file = "./inferCNV_annotions_se_pat_1_2_3.txt",sep = "\t", row.names = F, col.names = F)

    # extract 10x count data from used as input for the inferCNV run
    counts_matrix = GetAssayData(se, slot="counts")

# Create the infercnv object

    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                        annotations_file="./inferCNV_annotions_se_pat_1_2_3.txt",
                                        delim="\t",
                                        gene_order_file="./inferCNV_pediatric_patient_1/gencode.v25.annotation_gen_pos_v3.txt",
                                        ref_group_names=c("patient_2","patient_3"),
                                        chr_exclude=c("chrMT"))

# InferCNV run

    # perform infercnv operations to reveal cnv signal
    infercnv_obj = infercnv::run(infercnv_obj,
                                 cutoff=0.1,  
                                 out_dir="./inferCNV_pediatric_patient_1_output_dir",  # dir is auto-created for storing outputs
                                 cluster_by_groups=T,   # If observations are defined according to groups (ie. patients), each group will be clustered separately
                                 denoise=T,
                                 HMM=T)