# Seurat Spatial Import example

The data for [Erickson et
al](https://www.nature.com/articles/s41586-022-05023-2) can be found at
the following Mendeley link [(latest dataset version =
4)](https://data.mendeley.com/datasets/svw96g68dv/4). This version of
the dataset is not yet live as of 07.12.2022 but should be public
shortly.

The following code downloads the count matrix file, and the spaceranger
“spatial” folder files, and imports them into a
[Seurat](https://satijalab.org/seurat/index.html) object for further
analysis.

    #Install Seurat if not already installed
    #install.packages('Seurat')

    #Initialize the Seurat library
    library(Seurat)

    ## Warning: package 'Seurat' was built under R version 4.2.2

    ## Attaching SeuratObject

    #Downloading Patient 1 - H2_1 filtered_feature_bc_matrix.h5 file to working folder
    url = "https://data.mendeley.com/public-files/datasets/svw96g68dv/files/8b69170c-6c07-4e69-abf2-35fade0f5e2c/file_downloaded"
    download.file(url,'./filtered_feature_bc_matrix.h5', mode = 'wb')

    #Create subdirectory called "spatial"
    dir.create("spatial")

    ## Warning in dir.create("spatial"): 'spatial' already exists

    #07.12.2022 - This is manually downloaded for the user while waiting for Mendeley updates to be pushed
    #Downloading Patient 1 - H2_1 tissue_hires_image.png image file to spatial folder
    #url = "https://data.mendeley.com/api/datasets/svw96g68dv/draft/files/e1399690-dc45-43e5-ae39-7a065bf7d34e/file_downloaded"
    #download.file(url,'./spatial/H2_1_tissue_hires_image.png', mode = 'wb')

    #Downloading Patient 1 - H2_1 scalefactors_json.json file to spatial folder
    url = "https://data.mendeley.com/public-files/datasets/svw96g68dv/files/06eb7410-a6a3-4ea9-a364-c6a734a22169/file_downloaded"
    download.file(url,'./spatial/scalefactors_json.json', mode = 'wb')

    #Downloading Patient 1 - H2_1 tissue_positions_list.csv file to spatial folder
    url = "https://data.mendeley.com/public-files/datasets/svw96g68dv/files/e028d330-142b-4d8b-b32d-9114b5c48421/file_downloaded"
    download.file(url,'./spatial/tissue_positions_list.csv', mode = 'wb')

    #Creating an image file using the Seurat Read10X_image() function
    InputImage <- Read10X_Image(
      "./spatial",
      image.name = "H2_1_tissue_hires_image.png",
      filter.matrix = FALSE
    )

    #Loading a Seurat object using the Load10X_Spatial() function
    H2_1_Seurat <- Load10X_Spatial(
      ".",
      filename = "filtered_feature_bc_matrix.h5",
      assay = "Spatial",
      image = InputImage
    )

    #Confirming that the Seurat object was created and loaded.  
    summary(H2_1_Seurat)

    ## Length  Class   Mode 
    ##      1 Seurat     S4

    head(H2_1_Seurat)

    ##                       orig.ident nCount_Spatial nFeature_Spatial
    ## AAACAAGTATCTCCCA-1 SeuratProject           8758             2717
    ## AAACACCAATAACTGC-1 SeuratProject          13466             3889
    ## AAACAGCTTTCAGAAG-1 SeuratProject           9514             2511
    ## AAACAGGGTCTATATT-1 SeuratProject          15668             3601
    ## AAACAGTGTTCCTGGG-1 SeuratProject              0                0
    ## AAACATTTCCCGGATT-1 SeuratProject           5290             2211
    ## AAACCCGAACGAAATC-1 SeuratProject             27               26
    ## AAACCGGAAATGTTAA-1 SeuratProject              7                7
    ## AAACCGGGTAGGTACC-1 SeuratProject           9728             2781
    ## AAACCGTTCGTCCAGG-1 SeuratProject           3783             1660

    sessionInfo()

    ## R version 4.2.1 (2022-06-23 ucrt)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 22000)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.utf8 
    ## [2] LC_CTYPE=English_United States.utf8   
    ## [3] LC_MONETARY=English_United States.utf8
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.utf8    
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] SeuratObject_4.1.3 Seurat_4.3.0      
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] Rtsne_0.16             colorspace_2.0-3       deldir_1.0-6          
    ##   [4] ellipsis_0.3.2         ggridges_0.5.4         rstudioapi_0.14       
    ##   [7] spatstat.data_3.0-0    leiden_0.4.3           listenv_0.8.0         
    ##  [10] bit64_4.0.5            ggrepel_0.9.2          fansi_1.0.3           
    ##  [13] codetools_0.2-18       splines_4.2.1          knitr_1.40            
    ##  [16] polyclip_1.10-4        jsonlite_1.8.3         ica_1.0-3             
    ##  [19] cluster_2.1.3          png_0.1-7              uwot_0.1.14           
    ##  [22] shiny_1.7.3            sctransform_0.3.5      spatstat.sparse_3.0-0 
    ##  [25] compiler_4.2.1         httr_1.4.4             assertthat_0.2.1      
    ##  [28] Matrix_1.5-3           fastmap_1.1.0          lazyeval_0.2.2        
    ##  [31] cli_3.3.0              later_1.3.0            htmltools_0.5.2       
    ##  [34] tools_4.2.1            igraph_1.3.5           gtable_0.3.1          
    ##  [37] glue_1.6.2             RANN_2.6.1             reshape2_1.4.4        
    ##  [40] dplyr_1.0.10           Rcpp_1.0.9             scattermore_0.8       
    ##  [43] vctrs_0.5.1            nlme_3.1-157           spatstat.explore_3.0-5
    ##  [46] progressr_0.11.0       lmtest_0.9-40          spatstat.random_3.0-1 
    ##  [49] xfun_0.31              stringr_1.4.1          globals_0.16.2        
    ##  [52] mime_0.12              miniUI_0.1.1.1         lifecycle_1.0.3       
    ##  [55] irlba_2.3.5.1          goftest_1.2-3          future_1.29.0         
    ##  [58] MASS_7.3-57            zoo_1.8-11             scales_1.2.1          
    ##  [61] promises_1.2.0.1       spatstat.utils_3.0-1   parallel_4.2.1        
    ##  [64] RColorBrewer_1.1-3     yaml_2.3.5             reticulate_1.26       
    ##  [67] pbapply_1.6-0          gridExtra_2.3          ggplot2_3.4.0         
    ##  [70] stringi_1.7.8          rlang_1.0.6            pkgconfig_2.0.3       
    ##  [73] matrixStats_0.63.0     evaluate_0.18          lattice_0.20-45       
    ##  [76] ROCR_1.0-11            purrr_0.3.5            tensor_1.5            
    ##  [79] patchwork_1.1.2        htmlwidgets_1.5.4      bit_4.0.5             
    ##  [82] cowplot_1.1.1          tidyselect_1.2.0       parallelly_1.32.1     
    ##  [85] RcppAnnoy_0.0.20       plyr_1.8.8             magrittr_2.0.3        
    ##  [88] R6_2.5.1               generics_0.1.3         DBI_1.1.3             
    ##  [91] pillar_1.8.1           fitdistrplus_1.1-8     survival_3.3-1        
    ##  [94] abind_1.4-5            sp_1.5-1               tibble_3.1.8          
    ##  [97] future.apply_1.10.0    hdf5r_1.3.7            KernSmooth_2.23-20    
    ## [100] utf8_1.2.2             spatstat.geom_3.0-3    plotly_4.10.1         
    ## [103] rmarkdown_2.18         grid_4.2.1             data.table_1.14.6     
    ## [106] digest_0.6.29          xtable_1.8-4           tidyr_1.2.1           
    ## [109] httpuv_1.6.6           munsell_0.5.0          viridisLite_0.4.1
