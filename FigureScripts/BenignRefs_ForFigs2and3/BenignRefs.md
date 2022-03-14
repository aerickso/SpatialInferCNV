# Setup

Initiating libraries.

``` r
library(SpatialInferCNV)
library(devtools)
library(ape)
library(phylogram)
library(tidyverse)
```

# Download data

Download all the data from
[Mendeley](https://data.mendeley.com/v1/datasets/svw96g68dv/draft),
specifically all folders from: count_matrices/Patient
1/Visium_with_annotation/.

``` r
dir.create("Patient1_BenignRefs")
setwd("Patient1_BenignRefs")
```

# Selecting Benign Histogical Spot annotations

We then import the consensus pathology annotations and select benigns
only, and created an annotation dataframe.

``` r
H1_2_Cleaned <- ImportHistologicalAnnotations("H1_2", "./Patient1_BenignRefs/Visium_with_annotation/H1_2/H1_2_Final_Consensus_Annotations.csv")
H1_2_Benigns <- filter(H1_2_Cleaned, Histology == "Benign")
rm(H1_2_Cleaned)

H1_4_Cleaned <- ImportHistologicalAnnotations("H1_4", "./Patient1_BenignRefs/Visium_with_annotation/H1_4/H1_4_Final_Consensus_Annotations.csv")
H1_4_Benigns <- filter(H1_4_Cleaned, Histology == "Benign")
rm(H1_4_Cleaned)

H1_5_Cleaned <- ImportHistologicalAnnotations("H1_5", "./Patient1_BenignRefs/Visium_with_annotation/H1_5/H1_5_Final_Consensus_Annotations.csv")
H1_5_Benigns <- filter(H1_5_Cleaned, Histology == "Benign")
rm(H1_5_Cleaned)

H2_1_Cleaned <- ImportHistologicalAnnotations("H2_1", "./Patient1_BenignRefs/Visium_with_annotation/H2_1/H2_1_Final_Consensus_Annotations.csv")
H2_1_Benigns <- filter(H2_1_Cleaned, Histology == "Benign")
rm(H2_1_Cleaned)

H2_2_Cleaned <- ImportHistologicalAnnotations("H2_2","./Patient1_BenignRefs/Visium_with_annotation/H2_2/H2_2_Final_Consensus_Annotations.csv")
H2_2_Benigns <- filter(H2_2_Cleaned, Histology == "Benign")
rm(H2_2_Cleaned)

H2_5_Cleaned <- ImportHistologicalAnnotations("H2_5", "./Patient1_BenignRefs/Visium_with_annotation/H2_5/H2_5_Final_Consensus_Annotations.csv")
H2_5_Benigns <- filter(H2_5_Cleaned, Histology == "Benign")
rm(H2_5_Cleaned)

V1_2_Cleaned <- ImportHistologicalAnnotations("V1_2", "./Patient1_BenignRefs/Visium_with_annotation/V1_2/V1_2_Final_Consensus_Annotations.csv")
V1_2_Benigns <- filter(V1_2_Cleaned, Histology == "Benign" | Histology == "Benign*")
rm(V1_2_Cleaned)

AllBenigns <- rbind(H1_2_Benigns, H1_4_Benigns)
AllBenigns <- rbind(AllBenigns, H2_1_Benigns)
AllBenigns <- rbind(AllBenigns, H2_2_Benigns)
AllBenigns <- rbind(AllBenigns, H2_5_Benigns)
AllBenigns <- rbind(AllBenigns, V1_2_Benigns)

rm(H1_2_Benigns,
   H1_4_Benigns,
   H1_5_Benigns,
   H2_1_Benigns,
   H2_2_Benigns,
   H2_5_Benigns,
   V1_2_Benigns)

MergedAll <- AllBenigns
names(MergedAll)[2] <- "Histology"
rm(AllBenigns)
```

# Importing Count Data

This code chunk imports the .h5 files a default processed output from
[10x Genomics cell ranger pipeline
documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/molecule_info),
and appends a section label to the barcode.

We use the function ImportCountData(), which requires a section label,
and a path to the corresponding .h5 file. Again these are provided from
the Mendeley repository (as described above).

``` r
H2_1_ENSBMLID_Counts <- ImportCountData("H2_1", "./Patient1_BenignRefs/Visium_with_annotation/H2_1/filtered_feature_bc_matrix.h5")
H2_2_ENSBMLID_Counts <- ImportCountData("H2_2", "./Patient1_BenignRefs/Visium_with_annotation/H2_2/filtered_feature_bc_matrix.h5")
H1_2_ENSBMLID_Counts <- ImportCountData("H1_2", "./Patient1_BenignRefs/Visium_with_annotation/H1_2/filtered_feature_bc_matrix.h5")
H2_5_ENSBMLID_Counts <- ImportCountData("H2_5", "./Patient1_BenignRefs/Visium_with_annotation/H2_5/filtered_feature_bc_matrix.h5")
H1_4_ENSBMLID_Counts <- ImportCountData("H1_4", "./Patient1_BenignRefs/Visium_with_annotation/H1_4/filtered_feature_bc_matrix.h5")
V1_2_ENSBMLID_Counts <- ImportCountData("V1_2", "./Patient1_BenignRefs/Visium_with_annotation/V1_2/filtered_feature_bc_matrix.h5")
```

# QC, and Merging Count and Annotation Data

Next, we merge annotations with count data to get section wise count
matrices of only benign spots. This also applies a QC threshold (only
allowing spots with 500 UMIs or more to pass to the filtered
dataframes).

``` r
H2_1_Joined_Counts <- MergingCountAndAnnotationData("H2_1",MergedAll, H2_1_ENSBMLID_Counts)
H2_2_Joined_Counts <- MergingCountAndAnnotationData("H2_2",MergedAll, H2_2_ENSBMLID_Counts)
H1_2_Joined_Counts <- MergingCountAndAnnotationData("H1_2",MergedAll, H1_2_ENSBMLID_Counts)
H2_5_Joined_Counts <- MergingCountAndAnnotationData("H2_5",MergedAll, H2_5_ENSBMLID_Counts)
H1_4_Joined_Counts <- MergingCountAndAnnotationData("H1_4",MergedAll, H1_4_ENSBMLID_Counts)
V1_2_Joined_Counts <- MergingCountAndAnnotationData("V1_2",MergedAll, V1_2_ENSBMLID_Counts)

rm(H2_1_ENSBMLID_Counts, H2_2_ENSBMLID_Counts, H1_2_ENSBMLID_Counts, H2_5_ENSBMLID_Counts, H1_4_ENSBMLID_Counts, V1_2_ENSBMLID_Counts)
```

# Merging all count data into one object

We then merge all the sectionwise dataframes together, replace joined
NA’s with 0’s (inferCNV requires this), and output final count and
annotation .tsv files that are required for infercnv:run.

``` r
Counts_joined <- H2_1_Joined_Counts %>% full_join(H2_2_Joined_Counts, by = "Genes")
Counts_joined <- Counts_joined %>% full_join(H1_2_Joined_Counts, by = "Genes")
Counts_joined <- Counts_joined %>% full_join(H2_5_Joined_Counts, by = "Genes")
Counts_joined <- Counts_joined %>% full_join(H1_4_Joined_Counts, by = "Genes")
Counts_joined <- Counts_joined %>% full_join(V1_2_Joined_Counts, by = "Genes")

rm(H2_1_Joined_Counts ,H2_2_Joined_Counts, H1_2_Joined_Counts, H2_5_Joined_Counts, H1_4_Joined_Counts, V1_2_Joined_Counts)

Counts_joined <- Counts_joined %>% replace(., is.na(.), 0)
Counts_joined <- Counts_joined %>% column_to_rownames(., var = "Genes")

write.table(Counts_joined, "Organscale_Consensus_Benign_Counts.tsv", sep = "\t")

MergedAll_Final <- FinalAnnotations(MergedAll, Counts_joined)

write.table(MergedAll_Final, "Organscale_Consensus_Benign_Annotations.tsv", 
            sep = "\t",
            quote = FALSE, 
            col.names = FALSE, 
            row.names = FALSE)
```

# Confirming that the files are formatted correctly to create an inferCNV object

The siCNV_GeneOrderFile.tsv has been provided here:
<https://github.com/aerickso/SpatialInferCNV/tree/main/FigureScripts>.

``` r
AllBenigns_Consensus_Test_infCNV <- infercnv::CreateInfercnvObject(raw_counts_matrix="./Patient1_BenignRefs/Organscale_Consensus_Benign_Counts.tsv", 
                                               gene_order_file="./FigureScripts/siCNV_GeneOrderFile.tsv",
                                               annotations_file="./Patient1_BenignRefs/Organscale_Consensus_Benign_Annotations_04112020.tsv",
                                               delim="\t",
                                               ref_group_names=NULL)
```

# Running InferCNV (Unsupervised)

``` r
AllBenigns_Consensus_Test_infCNV = infercnv::run(AllBenigns_Consensus_Test_infCNV,
                                              cutoff=0.1,
                                              out_dir="./Patient1_BenignRefs/Outputs", 
                                              num_threads = 20,
                                              cluster_by_groups=FALSE, 
                                              denoise=TRUE,
                                              HMM=FALSE)
```

![](./Data/InferCNV%20Files/ExampleOutputs/BenignRefs/infercnv.21_denoised.png)

InferCNV will output many files. We are primarily interested in the
final “infercnv.21_denoised.png” file, as well as the text file
associated with the dendrogram associated with the hierarchical
clustering on the left hand side of the image
(infercnv.21_denoised.observations_dendrogram.txt).

# Importing dendrogram

Next, we want to import this dendrogram file fromo the above step:

``` r
Consensus_AllBenigns <- read.dendrogram(file = "./Patient1_BenignRefs/Outputs/infercnv.21_denoised.observations_dendrogram.txt")

Consensus_AllBenigns_phylo <- as.phylo(Consensus_AllBenigns)
```

# Visualizing dendrogram node numbers

``` r
my.subtrees = subtrees(Consensus_AllBenigns_phylo)  

png("Consensus_AllBenigns_phylo_Nodes.png",width=10000,height=2500, res = 300)
plot(Consensus_AllBenigns_phylo,show.tip.label = FALSE)
nodelabels(text=1:Consensus_AllBenigns_phylo$Nnode,node=1:Consensus_AllBenigns_phylo$Nnode+Ntip(Consensus_AllBenigns_phylo))
dev.off()
```

![](./Images/PurestBenigns.png)

![](./Data/Clone%20Selection/BenignRefs/Consensus_AllBenigns_phylo_Nodes_18112020.png)

![](./Images/NodeSelection_BenignRefs_Dendrogram_19032021.png)

# Node selection (Manual Task outside of R in an image editor)

Next, view the output .png file, which provides a (albeit cluttered)
labeling of the dendrogram tree nodes. Manually select individual nodes
that correspond with a distinct signal, in this case, nodes of visium
spots with little-to-no signal.

``` r
#3039 + 2560 
```

# Selecting clones in R

Next, after identifying the numerical nodes that correspond to
dendrogram branches that correspond with a given set of signals (aka,
clones), we then manually select these nodes in R, apply a label, then
join them all together and output as a .csv file for use as a
“Histologically Benign, inferCNV null” reference set to compare other
features of interest against.

``` r
Node3039 <- SelectingSubTreeData(my.subtrees, 3039) 
Node2560 <- SelectingSubTreeData(my.subtrees, 2560)

Merged <- rbind(Node3039, Node2560)

table(Merged$Node)

Merged$Node <- "Purest Benigns"
names(Merged)[2] <- "Histology"

write.csv(Merged, "Consensus_PurestBenigns.csv", row.names = FALSE)
```

The final file is provided at
[Mendeley](https://data.mendeley.com/v1/datasets/svw96g68dv/draft):
Count_matrices/Patient 1/Consensus Pathology.csv.
