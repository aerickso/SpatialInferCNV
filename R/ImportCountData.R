#' Importing Visium spatial transcriptomics count data from filtered_feature_bc_matrix.h5 file (output from SpaceRanger pipeline) and appending section name to barcodes
#'
#' ImportCountData()
#' 
#' @param SectionName A character string for section name.
#' @param InputCountFile A file path to a filtered_feature_bc_matrix.h5 file (output from 10X Genomics SpaceRanger pipeline)
#' @return A dataframe of barcodes with appended section names
#' @examples
#' ImportCountData("H2_1", "./filtered_feature_bc_matrix.h5")

ImportCountData <- function(SectionName, InputCountFile) {
  input <- Read10X_h5(InputCountFile, use.names = FALSE)
  input <- as.matrix(input)
  input <- as.data.frame(t(input))
  input <- rownames_to_column(input, "Barcode")
  input$Barcode <- paste0(SectionName, "_", input$Barcode)
  input$Barcode <- gsub("\\-", "\\.", input$Barcode)
  return(input)
}
