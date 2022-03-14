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

    library(reshape2)

    ## 
    ## Attaching package: 'reshape2'

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     smiths

    library(grid)
    library(svglite)

# Creating Working Directory for Step 2

    dir.create("Fig1_Step2")
    setwd("Fig1_Step2")

# Importing Data - Part 1

Note that Fig1d\_STOrganscale\_Selected\_Mapped\_Annotations.tsv was
generated in [Step
1](https://github.com/aerickso/SpatialInferCNV/blob/main/FigureScripts/Figure%201/Step1_PreprocessingToSpotLevelHMMs/Fig1D_Step1_PreprocessingToSpotLevelHMMs.md).
This also imports the barcode files for the original 1k array data.

    AllBarcodes <- read.delim("./Fig1_Step1/Fig1d_STOrganscale_Selected_Mapped_Annotations.tsv", 
                              sep = "\t",
                              header = FALSE)
    names(AllBarcodes)[1] <- "Barcode"
    names(AllBarcodes)[2] <- "Histology"

    WGSize <-  3049315783

    L2_Barcodes <- read.table("https://github.com/SpatialTranscriptomicsResearch/st_pipeline/raw/master/ids/1000L2_barcodes.txt", sep = "\t")
    names(L2_Barcodes)[1] <- "Barcode"
    names(L2_Barcodes)[2] <- "X"
    names(L2_Barcodes)[3] <- "Y"
    L2_Barcodes$ModfiedBarcode <- paste0(L2_Barcodes$X, "x", L2_Barcodes$Y)
    L2_Barcodes <- L2_Barcodes %>% 
          select(ModfiedBarcode)
    names(L2_Barcodes)[1] <- "Barcode"

# Importing Data and formatting - Part 2

Note that 17\_HMM\_predHMMi6.hmm\_mode-cells.pred\_cnv\_genes.dat was
generated in [Step
1](https://github.com/aerickso/SpatialInferCNV/blob/main/FigureScripts/Figure%201/Step1_PreprocessingToSpotLevelHMMs/Fig1D_Step1_PreprocessingToSpotLevelHMMs.md).

We then use the first threshold to optimize signal-to-noise for spatial
visualization across all genes with an inferred CNV, all genes genes
with an inferred CNV (for this dataset, 35%).

    CNV_Genes_Organscale <- read.delim("./Fig1_Step1/Outputs/17_HMM_predHMMi6.hmm_mode-cells.pred_cnv_genes.dat")
    Counted <- CNV_Genes_Organscale %>% group_by(gene) %>% tally()
    MaxLength <- as.numeric(nrow(Counted))
    CountPercentageThreshold_35Perc <- round(.35 * MaxLength,0)
    CountedThresholded <- Counted %>%
                            filter(n > CountPercentageThreshold_35Perc)
    CNV_Genes_Filtered <- inner_join(CNV_Genes_Organscale, CountedThresholded) %>%
                            select(-n)

# Extracting Sectionwise Data

Note that 17\_HMM\_predHMMi6.hmm\_mode-cells.pred\_cnv\_genes.dat was
generated in [Step
1](https://github.com/aerickso/SpatialInferCNV/blob/main/FigureScripts/Figure%201/Step1_PreprocessingToSpotLevelHMMs/Fig1D_Step1_PreprocessingToSpotLevelHMMs.md).

We then use the second threshold to optimize signal-to-noise for spatial
visualization across data within a section itself (for this dataset,
45%).

    CNV_Genes_Filtered <- CNV_Genes_Filtered %>% mutate(section = substr(cell_group_name, 1, 4))

    H1_1_Sectionwise_CNVsGenes_Counted <- ExtractSectionWise("H1_1", CNV_Genes_Filtered, AllBarcodes, 0.45)
    H1_2_Sectionwise_CNVsGenes_Counted <- ExtractSectionWise("H1_2", CNV_Genes_Filtered, AllBarcodes, 0.45)
    H1_3_Sectionwise_CNVsGenes_Counted <- ExtractSectionWise("H1_3", CNV_Genes_Filtered, AllBarcodes, 0.45)
    H1_4_Sectionwise_CNVsGenes_Counted <- ExtractSectionWise("H1_4", CNV_Genes_Filtered, AllBarcodes, 0.45)
    H1_5_Sectionwise_CNVsGenes_Counted <- ExtractSectionWise("H1_5", CNV_Genes_Filtered, AllBarcodes, 0.45)
    H2_1_Sectionwise_CNVsGenes_Counted <- ExtractSectionWise("H2_1", CNV_Genes_Filtered, AllBarcodes, 0.45)
    H2_2_Sectionwise_CNVsGenes_Counted <- ExtractSectionWise("H2_2", CNV_Genes_Filtered, AllBarcodes, 0.45)
    H2_3_Sectionwise_CNVsGenes_Counted <- ExtractSectionWise("H2_3", CNV_Genes_Filtered, AllBarcodes, 0.45)
    H2_4_Sectionwise_CNVsGenes_Counted <- ExtractSectionWise("H2_4", CNV_Genes_Filtered, AllBarcodes, 0.45)
    H2_5_Sectionwise_CNVsGenes_Counted <- ExtractSectionWise("H2_5", CNV_Genes_Filtered, AllBarcodes, 0.45)
    V1_1_Sectionwise_CNVsGenes_Counted <- ExtractSectionWise("V1_1", CNV_Genes_Filtered, AllBarcodes, 0.45)
    V1_2_Sectionwise_CNVsGenes_Counted <- ExtractSectionWise("V1_2", CNV_Genes_Filtered, AllBarcodes, 0.45)
    V1_3_Sectionwise_CNVsGenes_Counted <- ExtractSectionWise("V1_3", CNV_Genes_Filtered, AllBarcodes, 0.45)
    V1_4_Sectionwise_CNVsGenes_Counted <- ExtractSectionWise("V1_4", CNV_Genes_Filtered, AllBarcodes, 0.45)
    V1_5_Sectionwise_CNVsGenes_Counted <- ExtractSectionWise("V1_5", CNV_Genes_Filtered, AllBarcodes, 0.45)
    V2_1_Sectionwise_CNVsGenes_Counted <- ExtractSectionWise("V2_1", CNV_Genes_Filtered, AllBarcodes, 0.45)
    V2_2_Sectionwise_CNVsGenes_Counted <- ExtractSectionWise("V2_2", CNV_Genes_Filtered, AllBarcodes, 0.45)
    V2_3_Sectionwise_CNVsGenes_Counted <- ExtractSectionWise("V2_3", CNV_Genes_Filtered, AllBarcodes, 0.45)
    V2_4_Sectionwise_CNVsGenes_Counted <- ExtractSectionWise("V2_4", CNV_Genes_Filtered, AllBarcodes, 0.45)
    V2_5_Sectionwise_CNVsGenes_Counted <- ExtractSectionWise("V2_5", CNV_Genes_Filtered, AllBarcodes, 0.45)
    V2_6_Sectionwise_CNVsGenes_Counted <- ExtractSectionWise("V2_6", CNV_Genes_Filtered, AllBarcodes, 0.45)

    H1_1Max <- as.numeric(max(H1_1_Sectionwise_CNVsGenes_Counted$PercentageGenomeAltered))
    H1_2Max <- as.numeric(max(H1_2_Sectionwise_CNVsGenes_Counted$PercentageGenomeAltered))
    H1_3Max <- as.numeric(max(H1_3_Sectionwise_CNVsGenes_Counted$PercentageGenomeAltered))
    H1_4Max <- as.numeric(max(H1_4_Sectionwise_CNVsGenes_Counted$PercentageGenomeAltered))
    H1_5Max <- as.numeric(max(H1_5_Sectionwise_CNVsGenes_Counted$PercentageGenomeAltered))
    H2_1Max <- as.numeric(max(H2_1_Sectionwise_CNVsGenes_Counted$PercentageGenomeAltered))
    H2_2Max <- as.numeric(max(H2_2_Sectionwise_CNVsGenes_Counted$PercentageGenomeAltered))
    H2_3Max <- as.numeric(max(H2_3_Sectionwise_CNVsGenes_Counted$PercentageGenomeAltered))
    H2_4Max <- as.numeric(max(H2_4_Sectionwise_CNVsGenes_Counted$PercentageGenomeAltered))
    H2_5Max <- as.numeric(max(H2_5_Sectionwise_CNVsGenes_Counted$PercentageGenomeAltered))
    V1_1Max <- as.numeric(max(V1_1_Sectionwise_CNVsGenes_Counted$PercentageGenomeAltered))
    V1_2Max <- as.numeric(max(V1_2_Sectionwise_CNVsGenes_Counted$PercentageGenomeAltered))
    V1_3Max <- as.numeric(max(V1_3_Sectionwise_CNVsGenes_Counted$PercentageGenomeAltered))
    V1_4Max <- as.numeric(max(V1_4_Sectionwise_CNVsGenes_Counted$PercentageGenomeAltered))
    V1_5Max <- as.numeric(max(V1_5_Sectionwise_CNVsGenes_Counted$PercentageGenomeAltered))
    V2_1Max <- as.numeric(max(V2_1_Sectionwise_CNVsGenes_Counted$PercentageGenomeAltered))
    V2_2Max <- as.numeric(max(V2_2_Sectionwise_CNVsGenes_Counted$PercentageGenomeAltered))
    V2_3Max <- as.numeric(max(V2_3_Sectionwise_CNVsGenes_Counted$PercentageGenomeAltered))
    V2_4Max <- as.numeric(max(V2_4_Sectionwise_CNVsGenes_Counted$PercentageGenomeAltered))
    V2_5Max <- as.numeric(max(V2_5_Sectionwise_CNVsGenes_Counted$PercentageGenomeAltered))
    V2_6Max <- as.numeric(max(V2_6_Sectionwise_CNVsGenes_Counted$PercentageGenomeAltered))

    MaxVal<- max(H1_1Max, 
                 H1_2Max,
                 H1_3Max, 
                 H1_4Max,
                 H1_5Max,
                 H2_1Max, 
                 H2_2Max,
                 H2_3Max, 
                 H2_4Max,
                 H2_5Max,
                 V1_1Max, 
                 V1_2Max,
                 V1_3Max, 
                 V1_4Max,
                 V1_5Max,
                 V2_1Max, 
                 V2_2Max,
                 V2_3Max, 
                 V2_4Max,
                 V2_5Max,
                 V2_6Max)

# Visualizing Spatial Outputs (Figure 1d)

Having obtained sectionwise dataframes, we then plot the spatial plots,
per section.

    dir.create("./Figure1D_sectionoutputs")
    setwd("./Figure1D_sectionoutputs")

    PGA_Matrix <- Output_PGA_Visualization_MatrixGreyNA("H1_1", H1_1_Sectionwise_CNVsGenes_Counted, L2_Barcodes)
    Plot_PGA_Visualization_Matrix("H1_1", PGA_Matrix, MaxVal)

    PGA_Matrix <- Output_PGA_Visualization_MatrixGreyNA("H1_2", H1_2_Sectionwise_CNVsGenes_Counted, L2_Barcodes)
    Plot_PGA_Visualization_Matrix("H1_2", PGA_Matrix, MaxVal)

    PGA_Matrix <- Output_PGA_Visualization_MatrixGreyNA("H1_3", H1_3_Sectionwise_CNVsGenes_Counted, L2_Barcodes)
    Plot_PGA_Visualization_Matrix("H1_3", PGA_Matrix, MaxVal)

    PGA_Matrix <- Output_PGA_Visualization_MatrixGreyNA("H1_4", H1_4_Sectionwise_CNVsGenes_Counted, L2_Barcodes)
    Plot_PGA_Visualization_Matrix("H1_4", PGA_Matrix, MaxVal)

    PGA_Matrix <- Output_PGA_Visualization_MatrixGreyNA("H1_5", H1_5_Sectionwise_CNVsGenes_Counted, L2_Barcodes)
    Plot_PGA_Visualization_Matrix("H1_5", PGA_Matrix, MaxVal)

    PGA_Matrix <- Output_PGA_Visualization_MatrixGreyNA("H2_1", H2_1_Sectionwise_CNVsGenes_Counted, L2_Barcodes)
    Plot_PGA_Visualization_Matrix("H2_1", PGA_Matrix, MaxVal)

    PGA_Matrix <- Output_PGA_Visualization_MatrixGreyNA("H2_2", H2_2_Sectionwise_CNVsGenes_Counted, L2_Barcodes)
    Plot_PGA_Visualization_Matrix("H2_2", PGA_Matrix, MaxVal)

    PGA_Matrix <- Output_PGA_Visualization_MatrixGreyNA("H2_3", H2_3_Sectionwise_CNVsGenes_Counted, L2_Barcodes)
    Plot_PGA_Visualization_Matrix("H2_3", PGA_Matrix, MaxVal)

    PGA_Matrix <- Output_PGA_Visualization_MatrixGreyNA("H2_4", H2_4_Sectionwise_CNVsGenes_Counted, L2_Barcodes)
    Plot_PGA_Visualization_Matrix("H2_4", PGA_Matrix, MaxVal)

    #Note: this is also used for the Figure 1G "Zoom in"
    PGA_Matrix <- Output_PGA_Visualization_MatrixGreyNA("H2_5", H2_5_Sectionwise_CNVsGenes_Counted, L2_Barcodes)
    Plot_PGA_Visualization_Matrix("H2_5", PGA_Matrix, MaxVal)

    PGA_Matrix <- Output_PGA_Visualization_MatrixGreyNA("V1_1", V1_1_Sectionwise_CNVsGenes_Counted, L2_Barcodes)
    Plot_PGA_Visualization_Matrix("V1_1", PGA_Matrix, MaxVal)

    PGA_Matrix <- Output_PGA_Visualization_MatrixGreyNA("V1_2", V1_2_Sectionwise_CNVsGenes_Counted, L2_Barcodes)
    Plot_PGA_Visualization_Matrix("V1_2", PGA_Matrix, MaxVal)

    PGA_Matrix <- Output_PGA_Visualization_MatrixGreyNA("V1_3", V1_3_Sectionwise_CNVsGenes_Counted, L2_Barcodes)
    Plot_PGA_Visualization_Matrix("V1_3", PGA_Matrix, MaxVal)

    PGA_Matrix <- Output_PGA_Visualization_MatrixGreyNA("V1_4", V1_4_Sectionwise_CNVsGenes_Counted, L2_Barcodes)
    Plot_PGA_Visualization_Matrix("V1_4", PGA_Matrix, MaxVal)

    PGA_Matrix <- Output_PGA_Visualization_MatrixGreyNA("V1_5", V1_5_Sectionwise_CNVsGenes_Counted, L2_Barcodes)
    Plot_PGA_Visualization_Matrix("V1_5", PGA_Matrix, MaxVal)

    PGA_Matrix <- Output_PGA_Visualization_MatrixGreyNA("V2_1", V2_1_Sectionwise_CNVsGenes_Counted, L2_Barcodes)
    Plot_PGA_Visualization_Matrix("V2_1", PGA_Matrix, MaxVal)

    PGA_Matrix <- Output_PGA_Visualization_MatrixGreyNA("V2_2", V2_2_Sectionwise_CNVsGenes_Counted, L2_Barcodes)
    Plot_PGA_Visualization_Matrix("V2_2", PGA_Matrix, MaxVal)

    PGA_Matrix <- Output_PGA_Visualization_MatrixGreyNA("V2_3", V2_3_Sectionwise_CNVsGenes_Counted, L2_Barcodes)
    Plot_PGA_Visualization_Matrix("V2_3", PGA_Matrix, MaxVal)

    PGA_Matrix <- Output_PGA_Visualization_MatrixGreyNA("V2_4", V2_4_Sectionwise_CNVsGenes_Counted, L2_Barcodes)
    Plot_PGA_Visualization_Matrix("V2_4", PGA_Matrix, MaxVal)

    PGA_Matrix <- Output_PGA_Visualization_MatrixGreyNA("V2_5", V2_5_Sectionwise_CNVsGenes_Counted, L2_Barcodes)
    Plot_PGA_Visualization_Matrix("V2_5", PGA_Matrix, MaxVal)

    PGA_Matrix <- Output_PGA_Visualization_MatrixGreyNA("V2_6", V2_6_Sectionwise_CNVsGenes_Counted, L2_Barcodes)
    Plot_PGA_Visualization_Matrix("V2_6", PGA_Matrix, MaxVal)

Here is an an example output for section H2\_5:

![example output for section
H2\_5](https://github.com/aerickso/SpatialInferCNV/blob/main/FigureScripts/Figure%201/Step2_FigureImages/H2_5_Revised_PGA_SpatialVisualization_2022-02-28.png).

# Visualizing BarPlot (Figure 1G)

Finally, we visualize the Barplot in Figure 1G.

    All_Organ35_andSectionWise45 <- rbind(H1_1_Sectionwise_CNVsGenes_Counted,
      H1_2_Sectionwise_CNVsGenes_Counted,
      H1_3_Sectionwise_CNVsGenes_Counted,
      H1_4_Sectionwise_CNVsGenes_Counted,
      H1_5_Sectionwise_CNVsGenes_Counted,
      H2_1_Sectionwise_CNVsGenes_Counted,
      H2_2_Sectionwise_CNVsGenes_Counted,
      H2_3_Sectionwise_CNVsGenes_Counted,
      H2_4_Sectionwise_CNVsGenes_Counted,
      H2_5_Sectionwise_CNVsGenes_Counted,
      V1_1_Sectionwise_CNVsGenes_Counted,
      V1_2_Sectionwise_CNVsGenes_Counted,
      V1_3_Sectionwise_CNVsGenes_Counted,
      V1_4_Sectionwise_CNVsGenes_Counted,
      V1_5_Sectionwise_CNVsGenes_Counted,
      V2_1_Sectionwise_CNVsGenes_Counted,
      V2_2_Sectionwise_CNVsGenes_Counted,
      V2_3_Sectionwise_CNVsGenes_Counted,
      V2_4_Sectionwise_CNVsGenes_Counted,
      V2_5_Sectionwise_CNVsGenes_Counted,
      V2_6_Sectionwise_CNVsGenes_Counted)

    All_Organ35_andSectionWise45 <- All_Organ35_andSectionWise45 %>% 
                mutate(section = substr(Barcode, 1, 4))

    ForVis <- left_join(All_Organ35_andSectionWise45, AllBarcodes)

    x <- c("V1_1",
           "V1_2", 
           "V1_3",
           "V1_4",
           "V1_5",
           "V2_1",
           "V2_2", 
           "V2_3",
           "V2_4",
           "V2_5",
           "H2_1", 
           "H2_2",
           "H2_3",
           "H2_4",
           "H2_5",
           "H1_1", 
           "H1_2",
           "H1_3",
           "H1_4",
           "H1_5",
           "V2_6")

    ForVis <- ForVis %>%
      mutate(section =  factor(section, levels = x)) %>%
      arrange(section)
    names(ForVis)[2] <- "GenesWithInferredCNV"

    ggplot(ForVis,aes(x = Barcode, fill=GenesWithInferredCNV,y=GenesWithInferredCNV)) +
      geom_bar(stat="identity") +
      scale_fill_gradient(low="blue",high="yellow") +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      facet_wrap( ~ section, strip.position = "bottom", scales = "free_x") +
        theme(panel.spacing = unit(0, "lines"),
            strip.background = element_blank(),
            strip.placement = "outside") +
        xlab("x-axis label")

    ggsave(paste0("siCNV_SectionBarPlot_Figure1G", ".png"),
      plot = last_plot() + labs(x=NULL, y=NULL),
      device = NULL,
      path = NULL,
      scale = 1,
      width = NA,
      height = NA,
      dpi = 300,
      limitsize = TRUE)

Here is the output:

![Barplot in Figure
1G](https://github.com/aerickso/SpatialInferCNV/blob/main/FigureScripts/Figure%201/Step2_FigureImages/siCNV_SectionBarPlot_Figure1G.png).
